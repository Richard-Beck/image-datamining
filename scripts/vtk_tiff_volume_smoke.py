#!/usr/bin/env python
"""Smoke test for TIFF-backed 3D volume rendering with PyVista/VTK.

The script creates a small synthetic 3D volume, writes it as a TIFF stack,
loads it back with tifffile, and renders it offscreen to a PNG. It avoids any
CellposeSAM/CPSAM execution.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pyvista as pv
import tifffile


SOFTWARE_RENDERER_MARKERS = (
    "llvmpipe",
    "softpipe",
    "software rasterizer",
    "mesa x11",
    "swrast",
)


def synthetic_volume(shape: tuple[int, int, int]) -> np.ndarray:
    """Create a uint16 volume with a few bright structures."""
    z_size, y_size, x_size = shape
    z, y, x = np.indices(shape, dtype=np.float32)

    volume = np.zeros(shape, dtype=np.float32)
    structures = [
        (0.35 * x_size, 0.45 * y_size, 0.40 * z_size, 0.18 * x_size, 1.0),
        (0.62 * x_size, 0.58 * y_size, 0.58 * z_size, 0.14 * x_size, 0.8),
        (0.50 * x_size, 0.28 * y_size, 0.70 * z_size, 0.10 * x_size, 0.7),
    ]

    for cx, cy, cz, radius, intensity in structures:
        distance = np.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2)
        volume += intensity * np.clip(1.0 - distance / radius, 0.0, 1.0)

    volume += 0.08 * (z / max(z_size - 1, 1))
    volume = np.clip(volume / volume.max(), 0.0, 1.0)
    return (volume * np.iinfo(np.uint16).max).astype(np.uint16)


def render_volume(volume: np.ndarray, output_png: Path) -> dict[str, str]:
    """Render a z-y-x NumPy volume to a PNG with VTK via PyVista."""
    grid = pv.ImageData()
    grid.dimensions = volume.shape[::-1]
    grid.spacing = (1.0, 1.0, 1.5)
    grid.origin = (0.0, 0.0, 0.0)
    grid.point_data["intensity"] = np.ascontiguousarray(volume.transpose(2, 1, 0)).ravel(
        order="F"
    )

    plotter = pv.Plotter(off_screen=True, window_size=(900, 700))
    actor = plotter.add_volume(
        grid,
        scalars="intensity",
        mapper="gpu",
        cmap="viridis",
        opacity="sigmoid",
        shade=True,
    )
    plotter.add_axes()
    plotter.camera_position = "iso"
    plotter.show(screenshot=str(output_png), auto_close=False)
    capabilities = plotter.ren_win.ReportCapabilities()
    renderer_info = parse_opengl_capabilities(capabilities)
    renderer_info["vtk_mapper"] = actor.GetMapper().GetClassName()
    plotter.close()
    return renderer_info


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create and render a synthetic TIFF volume using PyVista/VTK."
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("outputs/vtk_tiff_smoke"),
        help="Directory for the generated TIFF and PNG.",
    )
    parser.add_argument(
        "--shape",
        type=int,
        nargs=3,
        metavar=("Z", "Y", "X"),
        default=(48, 96, 96),
        help="Synthetic volume shape in z y x order.",
    )
    parser.add_argument(
        "--allow-software-renderer",
        action="store_true",
        help="Do not fail when VTK reports a known software OpenGL renderer.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tiff_path = args.out_dir / "synthetic_volume.tif"
    png_path = args.out_dir / "synthetic_volume_render.png"

    source = synthetic_volume(tuple(args.shape))
    tifffile.imwrite(tiff_path, source, photometric="minisblack")

    loaded = tifffile.imread(tiff_path)
    if loaded.shape != source.shape:
        raise RuntimeError(f"Loaded TIFF shape {loaded.shape} != expected {source.shape}")

    renderer_info = render_volume(loaded, png_path)
    if not args.allow_software_renderer:
        assert_not_software_renderer(renderer_info)

    print(f"Wrote TIFF: {tiff_path}")
    print(f"Wrote render: {png_path}")
    print(f"Render file bytes: {png_path.stat().st_size}")
    print(f"Volume shape zyx: {loaded.shape}, dtype: {loaded.dtype}")
    for key in sorted(renderer_info):
        print(f"{key}: {renderer_info[key]}")


if __name__ == "__main__":
    main()
