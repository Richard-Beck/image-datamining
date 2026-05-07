#!/usr/bin/env python
"""Render the largest instance ID in a mask TIFF as a translucent surface."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv
import tifffile
from skimage import measure


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
            "Extract the largest nonzero instance ID from a mask TIFF, create a "
            "marching-cubes surface, and render translucent views with VTK."
        )
    )
    parser.add_argument("mask_tif", type=Path, help="Input 3D mask TIFF.")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to a sibling render_largest_surface directory.",
    )
    parser.add_argument(
        "--spacing",
        type=float,
        nargs=3,
        metavar=("Z", "Y", "X"),
        default=(1.0, 1.0, 1.0),
        help="Voxel spacing in z y x order for marching cubes coordinates.",
    )
    parser.add_argument(
        "--opacity",
        type=float,
        default=0.35,
        help="Surface opacity from 0 to 1.",
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
    return mask_tif.parent.parent / "render_largest_instance_surface" / mask_tif.stem


def load_mask(mask_tif: Path) -> np.ndarray:
    mask = tifffile.imread(mask_tif)
    if mask.ndim > 3:
        mask = np.squeeze(mask)
    if mask.ndim != 3:
        raise ValueError(f"Expected a 3D TIFF after squeeze, got shape {mask.shape}")
    return mask


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


def crop_with_padding(mask: np.ndarray, padding: int = 1) -> tuple[np.ndarray, tuple[slice, slice, slice]]:
    coords = np.argwhere(mask)
    lower = np.maximum(coords.min(axis=0) - padding, 0)
    upper = np.minimum(coords.max(axis=0) + padding + 1, mask.shape)
    crop_slices = tuple(slice(int(lo), int(hi)) for lo, hi in zip(lower, upper))
    return mask[crop_slices], crop_slices


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
    radius = float(np.linalg.norm(extent))
    if radius == 0:
        radius = 1.0

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
    out_dir: Path,
    opacity: float,
    window_size: tuple[int, int],
) -> tuple[dict[str, str], list[Path]]:
    renderer_info = {}
    output_paths = []
    for name, camera_position in camera_positions(mesh).items():
        plotter = pv.Plotter(off_screen=True, window_size=window_size)
        plotter.set_background("white")
        plotter.add_mesh(
            mesh,
            color="#2b8cbe",
            opacity=opacity,
            smooth_shading=True,
            specular=0.35,
            specular_power=20,
            show_edges=False,
        )
        plotter.add_axes()
        plotter.camera_position = camera_position
        plotter.camera.zoom(1.0)
        output_path = out_dir / f"{name}.png"
        plotter.show(screenshot=str(output_path), auto_close=False)
        if not renderer_info:
            renderer_info = parse_opengl_capabilities(plotter.ren_win.ReportCapabilities())
        plotter.close()
        output_paths.append(output_path)
    return renderer_info, output_paths


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir if args.out_dir is not None else default_out_dir(args.mask_tif)
    out_dir.mkdir(parents=True, exist_ok=True)

    mask = load_mask(args.mask_tif)
    component, largest_id, voxel_count = largest_instance(mask)
    cropped_component, crop_slices = crop_with_padding(component)
    mesh = mesh_from_component(cropped_component, crop_slices, tuple(args.spacing))
    mesh = mesh.smooth(n_iter=20, relaxation_factor=0.08)
    mesh.save(out_dir / "largest_component_surface.ply")

    renderer_info, output_paths = render_views(
        mesh=mesh,
        out_dir=out_dir,
        opacity=args.opacity,
        window_size=tuple(args.window_size),
    )
    if not args.allow_software_renderer:
        assert_not_software_renderer(renderer_info)

    print(f"Input mask: {args.mask_tif}")
    print(f"Mask shape zyx: {mask.shape}, dtype: {mask.dtype}")
    print(f"Largest instance ID: {largest_id}")
    print(f"Largest instance voxels: {voxel_count}")
    print(f"Crop slices zyx: {crop_slices}")
    print(f"Surface points: {mesh.n_points}, faces: {mesh.n_faces_strict}")
    print(f"Surface file: {out_dir / 'largest_component_surface.ply'}")
    for output_path in output_paths:
        print(f"Render file: {output_path} ({output_path.stat().st_size} bytes)")
    for key in sorted(renderer_info):
        print(f"{key}: {renderer_info[key]}")


if __name__ == "__main__":
    main()
