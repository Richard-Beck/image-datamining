#!/usr/bin/env python
"""Render FUCCI signal volume inside the largest segmented mask instance."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv
import tifffile
from scipy import ndimage as ndi
from skimage import filters, measure


SOFTWARE_RENDERER_MARKERS = (
    "llvmpipe",
    "softpipe",
    "software rasterizer",
    "mesa x11",
    "swrast",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render an opaque thresholded FUCCI volume from raw crop channels "
            "inside a translucent largest-instance mask surface."
        )
    )
    parser.add_argument("mask_tif", type=Path, help="Input 3D instance mask TIFF.")
    parser.add_argument("raw_crop_tif", type=Path, help="Raw crop TIFF in ZCYX order.")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults next to the mask run.",
    )
    parser.add_argument(
        "--fucci-channels",
        type=int,
        nargs=2,
        default=(0, 1),
        metavar=("CH_A", "CH_B"),
        help="Two raw ZCYX channels to normalize and combine by max.",
    )
    parser.add_argument(
        "--spacing",
        type=float,
        nargs=3,
        metavar=("Z", "Y", "X"),
        default=(1.0, 1.0, 1.0),
        help="Voxel spacing in z y x order.",
    )
    parser.add_argument(
        "--surface-opacity",
        type=float,
        default=0.12,
        help="Opacity for the segmented cell boundary.",
    )
    parser.add_argument(
        "--smooth-sigma",
        type=float,
        nargs=3,
        metavar=("Z", "Y", "X"),
        default=(0.6, 1.0, 1.0),
        help="Gaussian smoothing sigma for combined FUCCI signal.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Manual threshold on normalized/smoothed FUCCI signal. Defaults to Otsu inside the object.",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        default=(1200, 900),
        help="Render window size in pixels.",
    )
    parser.add_argument(
        "--allow-software-renderer",
        action="store_true",
        help="Do not fail when VTK reports a known software OpenGL renderer.",
    )
    return parser.parse_args()


def default_out_dir(mask_tif: Path) -> Path:
    return mask_tif.parent.parent / "render_largest_instance_fucci_volume" / mask_tif.stem


def load_mask(mask_tif: Path) -> np.ndarray:
    mask = tifffile.imread(mask_tif)
    if mask.ndim > 3:
        mask = np.squeeze(mask)
    if mask.ndim != 3:
        raise ValueError(f"Expected a 3D mask TIFF after squeeze, got shape {mask.shape}")
    return mask


def load_raw_crop(raw_crop_tif: Path) -> np.ndarray:
    raw = tifffile.imread(raw_crop_tif)
    if raw.ndim != 4:
        raise ValueError(f"Expected raw crop in ZCYX order, got shape {raw.shape}")
    return raw


def normalize01(values: np.ndarray) -> np.ndarray:
    values = values.astype(np.float32, copy=False)
    lo, hi = np.percentile(values, [1.0, 99.8])
    if hi <= lo:
        return np.zeros_like(values, dtype=np.float32)
    return np.clip((values - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)


def largest_instance(mask: np.ndarray) -> tuple[np.ndarray, int, int]:
    ids, counts = np.unique(mask, return_counts=True)
    foreground = ids != 0
    if not np.any(foreground):
        raise ValueError("Mask contains no nonzero instance IDs")
    ids = ids[foreground]
    counts = counts[foreground]
    largest_index = int(np.argmax(counts))
    largest_id = int(ids[largest_index])
    voxel_count = int(counts[largest_index])
    return mask == largest_id, largest_id, voxel_count


def crop_with_padding(mask: np.ndarray, padding: int = 2) -> tuple[tuple[slice, slice, slice], np.ndarray]:
    coords = np.argwhere(mask)
    lower = np.maximum(coords.min(axis=0) - padding, 0)
    upper = np.minimum(coords.max(axis=0) + padding + 1, mask.shape)
    crop_slices = tuple(slice(int(lo), int(hi)) for lo, hi in zip(lower, upper))
    return crop_slices, mask[crop_slices]


def mesh_from_component(
    component: np.ndarray,
    crop_slices: tuple[slice, slice, slice],
    spacing_zyx: tuple[float, float, float],
) -> pv.PolyData:
    padded = np.pad(component.astype(np.float32), 1, mode="constant")
    verts_zyx, faces, _normals, _values = measure.marching_cubes(
        padded,
        level=0.5,
        spacing=spacing_zyx,
    )
    crop_origin_zyx = np.array([sl.start for sl in crop_slices], dtype=np.float32)
    spacing = np.array(spacing_zyx, dtype=np.float32)
    verts_zyx += (crop_origin_zyx - 1.0) * spacing
    verts_xyz = verts_zyx[:, ::-1]

    pyvista_faces = np.empty((faces.shape[0], 4), dtype=np.int64)
    pyvista_faces[:, 0] = 3
    pyvista_faces[:, 1:] = faces
    return pv.PolyData(verts_xyz, pyvista_faces.ravel())


def build_fucci_volume(
    raw_zcyx: np.ndarray,
    component: np.ndarray,
    crop_slices: tuple[slice, slice, slice],
    channels: tuple[int, int],
    smooth_sigma: tuple[float, float, float],
    manual_threshold: float | None,
) -> tuple[np.ndarray, float, int]:
    if max(channels) >= raw_zcyx.shape[1] or min(channels) < 0:
        raise ValueError(f"Requested channels {channels}, but raw crop has {raw_zcyx.shape[1]} channels")

    ch_a = normalize01(raw_zcyx[:, channels[0]])
    ch_b = normalize01(raw_zcyx[:, channels[1]])
    combined = np.maximum(ch_a, ch_b)
    smoothed = ndi.gaussian_filter(combined, sigma=smooth_sigma).astype(np.float32)

    cropped_signal = smoothed[crop_slices]
    cropped_object = component[crop_slices]
    inside_values = cropped_signal[cropped_object]
    if inside_values.size == 0:
        raise ValueError("Largest object has no voxels after cropping")

    threshold = float(manual_threshold) if manual_threshold is not None else float(filters.threshold_otsu(inside_values))
    thresholded = np.zeros_like(cropped_signal, dtype=np.float32)
    keep = cropped_object & (cropped_signal >= threshold)
    thresholded[keep] = cropped_signal[keep]
    return thresholded, threshold, int(np.count_nonzero(keep))


def image_data_from_volume(
    volume_zyx: np.ndarray,
    crop_slices: tuple[slice, slice, slice],
    spacing_zyx: tuple[float, float, float],
    scalar_name: str = "values",
) -> pv.ImageData:
    spacing_xyz = tuple(spacing_zyx[::-1])
    origin_xyz = tuple(float(sl.start) * spacing for sl, spacing in zip(crop_slices[::-1], spacing_xyz))
    grid = pv.ImageData()
    grid.dimensions = volume_zyx.shape[::-1]
    grid.spacing = spacing_xyz
    grid.origin = origin_xyz
    grid.point_data[scalar_name] = np.ascontiguousarray(volume_zyx.transpose(2, 1, 0)).ravel(order="F")
    return grid


def parse_opengl_capabilities(capabilities: str) -> dict[str, str]:
    fields = {
        "OpenGL vendor string": "opengl_vendor",
        "OpenGL renderer string": "opengl_renderer",
        "OpenGL version string": "opengl_version",
    }
    parsed = {}
    for line in capabilities.splitlines():
        for label, key in fields.items():
            if line.startswith(label):
                parsed[key] = line.split(":", 1)[1].strip()
    return parsed


def assert_not_software_renderer(renderer_info: dict[str, str]) -> None:
    renderer = renderer_info.get("opengl_renderer", "")
    vendor = renderer_info.get("opengl_vendor", "")
    combined = f"{vendor} {renderer}".lower()
    if any(marker in combined for marker in SOFTWARE_RENDERER_MARKERS):
        raise RuntimeError(
            "VTK appears to be using a software OpenGL renderer: "
            f"vendor={vendor!r}, renderer={renderer!r}"
        )


def camera_positions(mesh: pv.PolyData) -> dict[str, tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]:
    center = np.array(mesh.center)
    bounds = np.array(mesh.bounds).reshape(3, 2)
    extent = bounds[:, 1] - bounds[:, 0]
    radius = float(np.linalg.norm(extent)) or 1.0
    views = {
        "front_near": ((0.0, -1.0, 0.25), 1.45),
        "front_far": ((0.0, -1.0, 0.25), 2.40),
        "right_oblique": ((1.0, -0.7, 0.35), 1.80),
        "left_high": ((-0.8, -0.6, 0.90), 1.95),
        "top_oblique": ((0.35, -0.35, 1.0), 1.65),
    }
    positions = {}
    for name, (direction, scale) in views.items():
        direction_arr = np.array(direction, dtype=np.float64)
        direction_arr /= np.linalg.norm(direction_arr)
        camera = center + direction_arr * radius * scale
        positions[name] = (tuple(camera), tuple(center), (0.0, 0.0, 1.0))
    return positions


def render_views(
    mesh: pv.PolyData,
    fucci_grid: pv.ImageData,
    out_dir: Path,
    surface_opacity: float,
    window_size: tuple[int, int],
) -> tuple[dict[str, str], list[Path], str]:
    renderer_info = {}
    output_paths = []
    volume_mapper = ""
    for name, camera_position in camera_positions(mesh).items():
        plotter = pv.Plotter(off_screen=True, window_size=window_size)
        plotter.set_background("white")
        plotter.add_mesh(
            mesh,
            color="#7b3294",
            opacity=surface_opacity,
            smooth_shading=True,
            specular=0.15,
            show_edges=False,
        )
        volume_actor = plotter.add_volume(
            fucci_grid,
            scalars="fucci",
            mapper="gpu",
            cmap="magma",
            opacity=[0.0, 1.0],
            shade=True,
        )
        plotter.add_axes()
        plotter.camera_position = camera_position
        output_path = out_dir / f"{name}.png"
        plotter.show(screenshot=str(output_path), auto_close=False)
        if not renderer_info:
            renderer_info = parse_opengl_capabilities(plotter.ren_win.ReportCapabilities())
            volume_mapper = volume_actor.GetMapper().GetClassName()
        plotter.close()
        output_paths.append(output_path)
    return renderer_info, output_paths, volume_mapper


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir if args.out_dir is not None else default_out_dir(args.mask_tif)
    out_dir.mkdir(parents=True, exist_ok=True)

    mask = load_mask(args.mask_tif)
    raw = load_raw_crop(args.raw_crop_tif)
    if raw.shape[0] != mask.shape[0] or raw.shape[2:] != mask.shape[1:]:
        raise ValueError(f"Raw shape {raw.shape} is incompatible with mask shape {mask.shape}")

    component, largest_id, object_voxels = largest_instance(mask)
    crop_slices, cropped_component = crop_with_padding(component)
    mesh = mesh_from_component(cropped_component, crop_slices, tuple(args.spacing))
    mesh = mesh.smooth(n_iter=20, relaxation_factor=0.08)
    mesh_path = out_dir / "largest_instance_boundary.ply"
    mesh.save(mesh_path)

    fucci_volume, threshold, fucci_voxels = build_fucci_volume(
        raw_zcyx=raw,
        component=component,
        crop_slices=crop_slices,
        channels=tuple(args.fucci_channels),
        smooth_sigma=tuple(args.smooth_sigma),
        manual_threshold=args.threshold,
    )
    fucci_tif_path = out_dir / "thresholded_fucci_inside_largest_instance.tif"
    tifffile.imwrite(fucci_tif_path, fucci_volume.astype(np.float32), metadata={"axes": "ZYX"})
    fucci_grid = image_data_from_volume(fucci_volume, crop_slices, tuple(args.spacing), scalar_name="fucci")

    renderer_info, output_paths, volume_mapper = render_views(
        mesh=mesh,
        fucci_grid=fucci_grid,
        out_dir=out_dir,
        surface_opacity=args.surface_opacity,
        window_size=tuple(args.window_size),
    )
    if not args.allow_software_renderer:
        assert_not_software_renderer(renderer_info)

    print(f"Mask: {args.mask_tif}")
    print(f"Raw crop: {args.raw_crop_tif}")
    print(f"Mask shape zyx: {mask.shape}, dtype: {mask.dtype}")
    print(f"Raw shape zcyx: {raw.shape}, dtype: {raw.dtype}")
    print(f"Largest instance ID: {largest_id}")
    print(f"Largest instance voxels: {object_voxels}")
    print(f"FUCCI channels: {tuple(args.fucci_channels)}")
    print(f"FUCCI threshold: {threshold:.6g}")
    print(f"Thresholded FUCCI voxels inside object: {fucci_voxels}")
    print(f"Crop slices zyx: {crop_slices}")
    print(f"Surface points: {mesh.n_points}, faces: {mesh.n_faces_strict}")
    print(f"Surface file: {mesh_path} ({mesh_path.stat().st_size} bytes)")
    print(f"Thresholded FUCCI TIFF: {fucci_tif_path} ({fucci_tif_path.stat().st_size} bytes)")
    for output_path in output_paths:
        print(f"Render file: {output_path} ({output_path.stat().st_size} bytes)")
    print(f"vtk_volume_mapper: {volume_mapper}")
    for key in sorted(renderer_info):
        print(f"{key}: {renderer_info[key]}")


if __name__ == "__main__":
    main()
