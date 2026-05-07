#!/usr/bin/env python
"""Open an interactive VTK window for FUCCI signal inside a mask boundary."""

from __future__ import annotations

import argparse
from pathlib import Path

import pyvista as pv

from render_instance_fucci_volume import (
    assert_not_software_renderer,
    build_fucci_volume,
    camera_positions,
    crop_with_padding,
    image_data_from_volume,
    largest_instance,
    load_mask,
    load_raw_crop,
    mesh_from_component,
    parse_opengl_capabilities,
)


DEFAULT_MASK = Path(
    "pipelines/fucci_3d_segmentation/runs/"
    "fucci_absz_sd2_product_cpsam_2d_fullstacks/stitched_simple/masks/"
    "MDAMB_240904_MDAMB_FUCCI_MDAMB453_80K_FoF3_fluorescent_nucleus_"
    "fucci_absz_sd2_product_stitched_simple_masks.tif"
)
DEFAULT_RAW = Path(
    "pipelines/fucci_3d_segmentation/runs/"
    "jackson_fucci_crop256_bf_fucci_fused_3d/raw_crops/"
    "MDAMB_240904_MDAMB_FUCCI_MDAMB453_80K_FoF3_fluorescent_nucleus_"
    "raw3ch_crop256_ZCYX.tif"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Open an interactive VTK render window for a FUCCI-in-cell overlay."
    )
    parser.add_argument("--mask-tif", type=Path, default=DEFAULT_MASK)
    parser.add_argument("--raw-crop-tif", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--fucci-channels", type=int, nargs=2, default=(0, 1))
    parser.add_argument("--spacing", type=float, nargs=3, default=(1.0, 1.0, 1.0), metavar=("Z", "Y", "X"))
    parser.add_argument("--smooth-sigma", type=float, nargs=3, default=(0.6, 1.0, 1.0), metavar=("Z", "Y", "X"))
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument("--surface-opacity", type=float, default=0.10)
    parser.add_argument("--window-size", type=int, nargs=2, default=(1400, 1000), metavar=("WIDTH", "HEIGHT"))
    parser.add_argument(
        "--allow-software-renderer",
        action="store_true",
        help="Do not fail when VTK reports a known software OpenGL renderer.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    mask = load_mask(args.mask_tif)
    raw = load_raw_crop(args.raw_crop_tif)
    if raw.shape[0] != mask.shape[0] or raw.shape[2:] != mask.shape[1:]:
        raise ValueError(f"Raw shape {raw.shape} is incompatible with mask shape {mask.shape}")

    component, largest_id, object_voxels = largest_instance(mask)
    crop_slices, cropped_component = crop_with_padding(component)
    mesh = mesh_from_component(cropped_component, crop_slices, tuple(args.spacing))
    mesh = mesh.smooth(n_iter=20, relaxation_factor=0.08)
    fucci_volume, threshold, fucci_voxels = build_fucci_volume(
        raw_zcyx=raw,
        component=component,
        crop_slices=crop_slices,
        channels=tuple(args.fucci_channels),
        smooth_sigma=tuple(args.smooth_sigma),
        manual_threshold=args.threshold,
    )
    fucci_grid = image_data_from_volume(fucci_volume, crop_slices, tuple(args.spacing), scalar_name="fucci")

    print(f"Mask: {args.mask_tif}")
    print(f"Raw crop: {args.raw_crop_tif}")
    print(f"Largest instance ID: {largest_id}")
    print(f"Largest instance voxels: {object_voxels}")
    print(f"FUCCI threshold: {threshold:.6g}")
    print(f"Thresholded FUCCI voxels inside object: {fucci_voxels}")

    plotter = pv.Plotter(off_screen=False, window_size=tuple(args.window_size))
    plotter.set_background("white")
    plotter.add_mesh(
        mesh,
        color="#7b3294",
        opacity=args.surface_opacity,
        smooth_shading=True,
        specular=0.15,
        show_edges=False,
        label="largest instance boundary",
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
    plotter.add_text("Largest instance boundary + thresholded FUCCI volume", font_size=10)
    plotter.camera_position = camera_positions(mesh)["right_oblique"]
    plotter.show(auto_close=False)

    renderer_info = parse_opengl_capabilities(plotter.ren_win.ReportCapabilities())
    if not args.allow_software_renderer:
        assert_not_software_renderer(renderer_info)
    print(f"vtk_volume_mapper: {volume_actor.GetMapper().GetClassName()}")
    for key in sorted(renderer_info):
        print(f"{key}: {renderer_info[key]}")
    plotter.close()


if __name__ == "__main__":
    main()
