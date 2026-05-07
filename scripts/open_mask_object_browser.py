#!/usr/bin/env python
"""Interactive VTK browser for stitched mask objects and FUCCI signal overlays."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyvista as pv

from render_instance_fucci_volume import (
    assert_not_software_renderer,
    build_fucci_volume,
    camera_positions,
    crop_with_padding,
    image_data_from_volume,
    load_mask,
    load_raw_crop,
    mesh_from_component,
    parse_opengl_capabilities,
)


DEFAULT_MASK = Path(
    "pipelines/fucci_3d_segmentation/runs/"
    "fucci_absz_sd2_product_cpsam_2d_fullstacks/stitched_simple/masks/"
    "MDAMB_240904_MDAMB_FUCCI_MDAMB231_30K_FoF3_fluorescent_nucleus_"
    "fucci_absz_sd2_product_stitched_simple_masks.tif"
)
DEFAULT_RAW = Path(
    "pipelines/fucci_3d_segmentation/runs/"
    "jackson_fucci_crop256_bf_fucci_fused_3d/raw_crops/"
    "MDAMB_240904_MDAMB_FUCCI_MDAMB231_30K_FoF3_fluorescent_nucleus_"
    "raw3ch_crop256_ZCYX.tif"
)

PALETTE = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#a6cee3",
    "#fb9a99",
    "#b2df8a",
    "#fdbf6f",
)

FUCCI_BASE_OPACITY = np.array([0.0, 0.05, 0.15, 0.30, 0.50], dtype=np.float32)
ABSZ_BASE_OPACITY = np.array([0.0, 0.10, 0.22, 0.40, 0.60], dtype=np.float32)


@dataclass(frozen=True)
class ObjectInfo:
    label: int
    voxels: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Open an interactive VTK window. The first view shows all mask objects "
            "above a voxel threshold as solid colored surfaces. Number keys "
            "highlight ranked objects, Enter switches the highlighted object to a "
            "translucent-boundary FUCCI volume view, and a/Escape restores the "
            "all-object view."
        )
    )
    parser.add_argument("--mask-tif", type=Path, default=DEFAULT_MASK)
    parser.add_argument("--raw-crop-tif", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--min-voxels", type=int, default=10_000)
    parser.add_argument("--fucci-channels", type=int, nargs=2, default=(0, 1))
    parser.add_argument("--absz-channel", type=int, default=2)
    parser.add_argument("--absz-min-z", type=float, default=2.5)
    parser.add_argument("--absz-percentile", type=float, default=98.0)
    parser.add_argument("--spacing", type=float, nargs=3, default=(1.0, 1.0, 1.0), metavar=("Z", "Y", "X"))
    parser.add_argument("--smooth-sigma", type=float, nargs=3, default=(0.6, 1.0, 1.0), metavar=("Z", "Y", "X"))
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument("--surface-opacity", type=float, default=0.10)
    parser.add_argument("--window-size", type=int, nargs=2, default=(1500, 1000), metavar=("WIDTH", "HEIGHT"))
    parser.add_argument(
        "--allow-software-renderer",
        action="store_true",
        help="Do not fail when VTK reports a known software OpenGL renderer.",
    )
    return parser.parse_args()


def object_infos(mask: np.ndarray, min_voxels: int) -> list[ObjectInfo]:
    ids, counts = np.unique(mask, return_counts=True)
    keep = (ids != 0) & (counts >= min_voxels)
    infos = [ObjectInfo(int(label), int(voxels)) for label, voxels in zip(ids[keep], counts[keep])]
    return sorted(infos, key=lambda item: item.voxels, reverse=True)


def build_absz_volume(
    raw_zcyx: np.ndarray,
    component: np.ndarray,
    crop_slices: tuple[slice, slice, slice],
    channel: int,
    min_abs_z: float,
    percentile: float,
) -> tuple[np.ndarray, float, int]:
    if channel < 0 or channel >= raw_zcyx.shape[1]:
        raise ValueError(f"Requested channel {channel}, but raw crop has {raw_zcyx.shape[1]} channels")

    raw_channel = raw_zcyx[:, channel].astype(np.float32, copy=False)
    inside_values = raw_channel[component]
    if inside_values.size == 0:
        raise ValueError("Selected object has no voxels")

    mean = float(inside_values.mean())
    std = float(inside_values.std())
    if std == 0.0:
        return np.zeros_like(raw_channel[crop_slices], dtype=np.float32), float("inf"), 0

    abs_z = np.abs((raw_channel - mean) / std).astype(np.float32)
    cropped_abs_z = abs_z[crop_slices]
    cropped_object = component[crop_slices]
    inside_abs_z = cropped_abs_z[cropped_object]
    threshold = max(float(min_abs_z), float(np.percentile(inside_abs_z, percentile)))

    thresholded = np.zeros_like(cropped_abs_z, dtype=np.float32)
    keep = cropped_object & (cropped_abs_z >= threshold)
    thresholded[keep] = cropped_abs_z[keep]
    return thresholded, threshold, int(np.count_nonzero(keep))


class ObjectBrowser:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.mask = load_mask(args.mask_tif)
        self.raw = load_raw_crop(args.raw_crop_tif)
        if self.raw.shape[0] != self.mask.shape[0] or self.raw.shape[2:] != self.mask.shape[1:]:
            raise ValueError(f"Raw shape {self.raw.shape} is incompatible with mask shape {self.mask.shape}")

        self.spacing = tuple(args.spacing)
        self.infos = object_infos(self.mask, args.min_voxels)
        if not self.infos:
            raise ValueError(f"No nonzero mask objects have at least {args.min_voxels} voxels")

        self.plotter = pv.Plotter(off_screen=False, window_size=tuple(args.window_size))
        self.plotter.set_background("white")
        self.mode = "all"
        self.last_camera_position = None
        self.selected_label: int | None = None
        self.highlight_actor = None
        self.current_label: int | None = None
        self.absz_min_z = float(args.absz_min_z)
        self.absz_percentile = float(args.absz_percentile)
        self.fucci_opacity_scale = 1.0
        self.absz_opacity_scale = 1.0
        self.surface_opacity = float(args.surface_opacity)

    def component_for_label(self, label: int) -> np.ndarray:
        return self.mask == label

    def surface_for_label(self, label: int, smooth: bool = True) -> pv.PolyData:
        component = self.component_for_label(label)
        crop_slices, cropped_component = crop_with_padding(component)
        mesh = mesh_from_component(cropped_component, crop_slices, self.spacing)
        if smooth:
            mesh = mesh.smooth(n_iter=20, relaxation_factor=0.08)
        return mesh

    def reset_scene(self) -> None:
        self.plotter.clear()
        self.plotter.set_background("white")
        self.plotter.add_axes()

    def fucci_opacity(self) -> list[float]:
        return np.clip(FUCCI_BASE_OPACITY * self.fucci_opacity_scale, 0.0, 1.0).tolist()

    def absz_opacity(self) -> list[float]:
        return np.clip(ABSZ_BASE_OPACITY * self.absz_opacity_scale, 0.0, 1.0).tolist()

    def rerender_current_detail(self) -> None:
        if self.mode == "detail" and self.current_label is not None:
            self.show_detail(self.current_label, recenter_camera=False)

    def camera_recentered_on_mesh(
        self,
        mesh: pv.PolyData,
        camera_position: tuple[
            tuple[float, float, float],
            tuple[float, float, float],
            tuple[float, float, float],
        ]
        | None,
    ) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
        new_focus = np.array(mesh.center, dtype=np.float64)
        if camera_position is None:
            base = camera_positions(mesh)["right_oblique"]
            camera = np.array(base[0], dtype=np.float64)
            focus = np.array(base[1], dtype=np.float64)
            view_up = base[2]
        else:
            camera = np.array(camera_position[0], dtype=np.float64)
            focus = np.array(camera_position[1], dtype=np.float64)
            view_up = camera_position[2]
        offset = camera - focus
        return tuple(new_focus + offset), tuple(new_focus), tuple(view_up)

    def show_all_objects(self) -> None:
        self.mode = "all"
        self.selected_label = None
        self.current_label = None
        self.highlight_actor = None
        self.reset_scene()

        for index, info in enumerate(self.infos):
            mesh = self.surface_for_label(info.label, smooth=True)
            self.plotter.add_mesh(
                mesh,
                color=PALETTE[index % len(PALETTE)],
                opacity=1.0,
                smooth_shading=True,
                specular=0.2,
                show_edges=False,
            )

        self.plotter.add_text(
            (
                f"Solid objects >= {self.args.min_voxels:,} voxels. "
                "Press 1-9 to highlight, Enter for detail, a/Escape to restore."
            ),
            font_size=10,
            name="mode_label",
        )
        if self.last_camera_position is not None:
            self.plotter.camera_position = self.last_camera_position
        else:
            self.plotter.view_isometric()
            self.plotter.reset_camera()
            self.plotter.camera.zoom(1.2)
        self.plotter.render()
        print(f"Showing {len(self.infos)} solid objects >= {self.args.min_voxels} voxels")

    def highlight_label(self, label: int, source: str) -> None:
        if self.mode != "all":
            print(f"{source}: detail mode active, press a/Escape to restore all objects first", flush=True)
            return
        if self.highlight_actor is not None:
            self.plotter.remove_actor(self.highlight_actor, render=False)

        mesh = self.surface_for_label(label, smooth=False)
        self.highlight_actor = self.plotter.add_mesh(
            mesh,
            color="red",
            style="wireframe",
            line_width=3,
            opacity=1.0,
            pickable=False,
        )
        self.selected_label = label
        self.plotter.add_text(
            f"Selected object {label}. Press Enter for FUCCI view, 1-9 to change, a/Escape to clear.",
            font_size=10,
            name="mode_label",
        )
        self.plotter.render()
        print(f"{source}: highlighted object {label}; press Enter to switch to FUCCI view", flush=True)

    def show_detail(self, label: int, recenter_camera: bool = True) -> None:
        self.mode = "detail"
        camera_position = self.plotter.camera_position
        self.last_camera_position = camera_position
        self.current_label = label
        self.selected_label = label
        self.highlight_actor = None
        self.reset_scene()

        component = self.component_for_label(label)
        crop_slices, cropped_component = crop_with_padding(component)
        mesh = mesh_from_component(cropped_component, crop_slices, self.spacing)
        mesh = mesh.smooth(n_iter=20, relaxation_factor=0.08)
        fucci_volume, threshold, fucci_voxels = build_fucci_volume(
            raw_zcyx=self.raw,
            component=component,
            crop_slices=crop_slices,
            channels=tuple(self.args.fucci_channels),
            smooth_sigma=tuple(self.args.smooth_sigma),
            manual_threshold=self.args.threshold,
        )
        fucci_grid = image_data_from_volume(fucci_volume, crop_slices, self.spacing, scalar_name="fucci")
        absz_volume, absz_threshold, absz_voxels = build_absz_volume(
            raw_zcyx=self.raw,
            component=component,
            crop_slices=crop_slices,
            channel=self.args.absz_channel,
            min_abs_z=self.absz_min_z,
            percentile=self.absz_percentile,
        )
        absz_grid = image_data_from_volume(absz_volume, crop_slices, self.spacing, scalar_name="absz")

        self.plotter.add_mesh(
            mesh,
            color="#7b3294",
            opacity=self.surface_opacity,
            smooth_shading=True,
            specular=0.15,
            show_edges=False,
        )
        volume_actor = self.plotter.add_volume(
            fucci_grid,
            scalars="fucci",
            mapper="gpu",
            cmap="Reds",
            opacity=self.fucci_opacity(),
            shade=False,
        )
        absz_actor = self.plotter.add_volume(
            absz_grid,
            scalars="absz",
            mapper="gpu",
            cmap="Blues",
            opacity=self.absz_opacity(),
            shade=False,
        )
        self.plotter.add_text(
            (
                f"Object {label}: FUCCI + ch{self.args.absz_channel:02d} abs-z. "
                f"z>={self.absz_min_z:.2f}, p{self.absz_percentile:.1f}. "
                "Use [] -/= u/U b/B s/S."
            ),
            font_size=10,
            name="mode_label",
        )
        if camera_position is not None and recenter_camera:
            self.plotter.camera_position = self.camera_recentered_on_mesh(mesh, camera_position)
        elif camera_position is not None:
            self.plotter.camera_position = camera_position
        else:
            self.plotter.camera_position = self.camera_recentered_on_mesh(mesh, None)
        self.plotter.render()
        print(
            f"Selected object {label}: object_voxels={int(np.count_nonzero(component))}, "
            f"fucci_threshold={threshold:.6g}, fucci_voxels={fucci_voxels}, "
            f"absz_threshold={absz_threshold:.6g}, absz_voxels={absz_voxels}, "
            f"fucci_opacity_scale={self.fucci_opacity_scale:.2f}, "
            f"absz_opacity_scale={self.absz_opacity_scale:.2f}, "
            f"surface_opacity={self.surface_opacity:.2f}, "
            f"vtk_volume_mapper={volume_actor.GetMapper().GetClassName()}, "
            f"absz_mapper={absz_actor.GetMapper().GetClassName()}"
        )

    def select_label(self, label: int, source: str) -> None:
        self.highlight_label(label, source)

    def select_rank(self, rank_index: int, source: str) -> None:
        if self.mode == "detail":
            print(f"{source}: detail mode active, restoring all objects", flush=True)
            self.show_all_objects()
            return
        if rank_index >= len(self.infos):
            print(f"{source}: no object at rank {rank_index + 1}", flush=True)
            return
        self.select_label(self.infos[rank_index].label, source)

    def restore_all_from_key(self) -> None:
        if self.mode == "detail":
            print("restore: restoring all objects", flush=True)
            self.show_all_objects()
            return
        if self.highlight_actor is not None:
            self.plotter.remove_actor(self.highlight_actor, render=False)
            self.highlight_actor = None
        self.selected_label = None
        self.plotter.add_text(
            f"Solid objects >= {self.args.min_voxels:,} voxels. Press 1-9 to highlight, Enter for detail.",
            font_size=10,
            name="mode_label",
        )
        self.plotter.render()
        print("restore: cleared selection", flush=True)

    def commit_selected_label(self) -> None:
        if self.mode == "detail":
            print("enter: already in FUCCI detail view", flush=True)
            return
        if self.selected_label is None:
            print("enter: no highlighted object; press 1-9 first", flush=True)
            return
        print(f"enter: switching object {self.selected_label} to FUCCI view", flush=True)
        self.show_detail(self.selected_label)

    def adjust_absz_percentile(self, delta: float) -> None:
        self.absz_percentile = float(np.clip(self.absz_percentile + delta, 50.0, 99.9))
        print(f"abs-z percentile -> {self.absz_percentile:.1f}", flush=True)
        self.rerender_current_detail()

    def adjust_absz_min_z(self, delta: float) -> None:
        self.absz_min_z = max(0.0, self.absz_min_z + delta)
        print(f"abs-z min z -> {self.absz_min_z:.2f}", flush=True)
        self.rerender_current_detail()

    def adjust_fucci_opacity(self, factor: float) -> None:
        self.fucci_opacity_scale = float(np.clip(self.fucci_opacity_scale * factor, 0.05, 4.0))
        print(f"FUCCI opacity scale -> {self.fucci_opacity_scale:.2f}", flush=True)
        self.rerender_current_detail()

    def adjust_absz_opacity(self, factor: float) -> None:
        self.absz_opacity_scale = float(np.clip(self.absz_opacity_scale * factor, 0.05, 4.0))
        print(f"abs-z opacity scale -> {self.absz_opacity_scale:.2f}", flush=True)
        self.rerender_current_detail()

    def adjust_surface_opacity(self, factor: float) -> None:
        self.surface_opacity = float(np.clip(self.surface_opacity * factor, 0.01, 1.0))
        print(f"surface opacity -> {self.surface_opacity:.2f}", flush=True)
        self.rerender_current_detail()

    def open(self) -> None:
        print(f"Mask: {self.args.mask_tif}")
        print(f"Raw crop: {self.args.raw_crop_tif}")
        print(f"Mask shape zyx: {self.mask.shape}, dtype: {self.mask.dtype}")
        print(f"Raw shape zcyx: {self.raw.shape}, dtype: {self.raw.dtype}")
        print(f"Objects shown: {len(self.infos)} with min_voxels={self.args.min_voxels}")
        print("Top objects:", ", ".join(f"{i + 1}={info.label}:{info.voxels}" for i, info in enumerate(self.infos[:10])))
        print("Controls:")
        print("  1-9: highlight object by voxel-volume rank")
        print("  Enter: switch highlighted object to detail view")
        print("  a or Escape: clear selection or restore all-object view")
        print("  [ / ]: decrease/increase abs-z percentile threshold")
        print("  - / =: decrease/increase abs-z minimum z threshold")
        print("  u / U: decrease/increase FUCCI opacity")
        print("  b / B: decrease/increase abs-z opacity")
        print("  s / S: decrease/increase boundary opacity")
        print(
            f"Abs-z layer: channel={self.args.absz_channel}, "
            f"threshold=max({self.absz_min_z}, p{self.absz_percentile} inside object)"
        )

        self.show_all_objects()
        self.plotter.add_key_event("a", self.restore_all_from_key)
        self.plotter.add_key_event("Escape", self.restore_all_from_key)
        self.plotter.add_key_event("Return", self.commit_selected_label)
        self.plotter.add_key_event("Enter", self.commit_selected_label)
        self.plotter.add_key_event("KP_Enter", self.commit_selected_label)
        for rank in range(1, 10):
            self.plotter.add_key_event(str(rank), lambda rank=rank: self.select_rank(rank - 1, source=f"key-{rank}"))
        self.plotter.add_key_event("[", lambda: self.adjust_absz_percentile(-1.0))
        self.plotter.add_key_event("]", lambda: self.adjust_absz_percentile(1.0))
        self.plotter.add_key_event("minus", lambda: self.adjust_absz_min_z(-0.25))
        self.plotter.add_key_event("equal", lambda: self.adjust_absz_min_z(0.25))
        self.plotter.add_key_event("-", lambda: self.adjust_absz_min_z(-0.25))
        self.plotter.add_key_event("=", lambda: self.adjust_absz_min_z(0.25))
        self.plotter.add_key_event("u", lambda: self.adjust_fucci_opacity(0.8))
        self.plotter.add_key_event("U", lambda: self.adjust_fucci_opacity(1.25))
        self.plotter.add_key_event("b", lambda: self.adjust_absz_opacity(0.8))
        self.plotter.add_key_event("B", lambda: self.adjust_absz_opacity(1.25))
        self.plotter.add_key_event("s", lambda: self.adjust_surface_opacity(0.8))
        self.plotter.add_key_event("S", lambda: self.adjust_surface_opacity(1.25))
        self.plotter.show(interactive=False, auto_close=False)
        self.plotter.iren.initialize()
        self.plotter.iren.start()

        renderer_info = parse_opengl_capabilities(self.plotter.ren_win.ReportCapabilities())
        if not self.args.allow_software_renderer:
            assert_not_software_renderer(renderer_info)
        for key in sorted(renderer_info):
            print(f"{key}: {renderer_info[key]}")
        self.plotter.close()


def main() -> None:
    args = parse_args()
    ObjectBrowser(args).open()


if __name__ == "__main__":
    main()
