#!/usr/bin/env python3
"""Track one K00 CPSAM mask stack with btrack."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import time
from pathlib import Path

import btrack
import numpy as np
import tifffile as tf
from btrack.config import TrackerConfig
from btrack.io.utils import localizations_to_objects
from btrack.models import MotionModel
from skimage.measure import regionprops


DEFAULT_OUT_DIR = Path("analyses/K00_GemcitabineExposure_033023/btrack_full_stacks")
SITE_RE = re.compile(r"^(?P<site_id>.+?)(?:_cpsam_masks)?\.tiff?$", re.IGNORECASE)

RAW_FEATURES = [
    "log_area",
    "log_major_axis_length",
    "log_minor_axis_length",
    "log_equivalent_diameter",
    "eccentricity",
    "solidity",
    "extent",
    "circularity",
    "orientation_cos2",
    "orientation_sin2",
]
TRACKING_FEATURES = [f"{name}_scaled" for name in RAW_FEATURES]
SIZE_TRACKING_FEATURES = ["log_area_scaled", "log_equivalent_diameter_scaled"]


def unique_fieldnames(names: list[str]) -> list[str]:
    seen: set[str] = set()
    unique: list[str] = []
    for name in names:
        if name not in seen:
            unique.append(name)
            seen.add(name)
    return unique


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run btrack on one CPSAM labeled mask stack. Defaults use "
            "size-only visual matching with a bounded spatial search radius."
        )
    )
    parser.add_argument("--mask-tiff", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--max-search-radius", type=float, default=40.0)
    parser.add_argument("--max-lost", type=int, default=1)
    parser.add_argument("--prob-not-assign", type=float, default=0.1)
    parser.add_argument("--accuracy", type=float, default=2.0)
    parser.add_argument("--step-size", type=int, default=100)
    parser.add_argument(
        "--motion-model",
        choices=("random-walk", "constant-velocity"),
        default="random-walk",
        help="Use random-walk for erratic cells; constant-velocity is retained for comparison.",
    )
    parser.add_argument(
        "--tracking-mode",
        choices=("motion", "motion-visual", "visual"),
        default="visual",
        help="Which btrack update terms to use. Use motion-visual only after motion-only looks sane.",
    )
    parser.add_argument(
        "--tracking-feature-set",
        choices=("size", "shape"),
        default="size",
        help="Feature set used when tracking-mode includes visual updates.",
    )
    parser.add_argument("--start-frame", type=int, default=0, help="Zero-based frame index to start tracking from.")
    parser.add_argument("--max-frames", type=int, default=None, help="Track only N mask frames after --start-frame.")
    parser.add_argument("--num-workers", type=int, default=1, help="Reserved for symmetry with Slurm settings.")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def site_id_from_mask(path: Path) -> str:
    match = SITE_RE.match(path.name)
    return match.group("site_id") if match else path.stem


def read_mask_stack(path: Path) -> np.ndarray:
    with tf.TiffFile(path) as tif:
        masks = tif.asarray()
        axes = tif.series[0].axes if tif.series else ""
    if masks.ndim == 2:
        masks = masks[None, :, :]
    if masks.ndim != 3:
        raise ValueError(f"Expected a TYX mask stack, got shape {masks.shape} from {path}")
    if axes and axes[-2:] != "YX":
        raise ValueError(f"Expected YX as the final TIFF axes, got axes={axes!r} from {path}")
    return masks


def equivalent_diameter(prop) -> float:
    if hasattr(prop, "equivalent_diameter_area"):
        return float(prop.equivalent_diameter_area)
    return float(prop.equivalent_diameter)


def object_rows_from_masks(masks: np.ndarray) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for frame, frame_mask in enumerate(masks):
        for prop in regionprops(frame_mask):
            area = float(prop.area)
            major = float(prop.major_axis_length)
            minor = float(prop.minor_axis_length)
            equiv = equivalent_diameter(prop)
            perimeter = float(prop.perimeter)
            circularity = (4.0 * math.pi * area / (perimeter * perimeter)) if perimeter > 0 else 0.0
            orientation = float(prop.orientation)

            rows.append(
                {
                    "t": int(frame),
                    "x": float(prop.centroid[1]),
                    "y": float(prop.centroid[0]),
                    "z": 0.0,
                    "label": int(prop.label),
                    "mask_label": int(prop.label),
                    "area": area,
                    "major_axis_length": major,
                    "minor_axis_length": minor,
                    "equivalent_diameter": equiv,
                    "eccentricity": float(prop.eccentricity),
                    "solidity": float(prop.solidity),
                    "extent": float(prop.extent),
                    "circularity": float(np.clip(circularity, 0.0, 1.0)),
                    "orientation": orientation,
                    "log_area": math.log1p(area),
                    "log_major_axis_length": math.log1p(major),
                    "log_minor_axis_length": math.log1p(minor),
                    "log_equivalent_diameter": math.log1p(equiv),
                    "orientation_cos2": 0.5 * (math.cos(2.0 * orientation) + 1.0),
                    "orientation_sin2": 0.5 * (math.sin(2.0 * orientation) + 1.0),
                }
            )
    return rows


def add_scaled_features(rows: list[dict[str, float | int]]) -> None:
    if not rows:
        return
    for feature in RAW_FEATURES:
        values = np.asarray([float(row[feature]) for row in rows], dtype=np.float64)
        lo, hi = np.nanpercentile(values, [1.0, 99.0])
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            scaled = np.full(values.shape, 0.5, dtype=np.float64)
        else:
            scaled = np.clip((values - lo) / (hi - lo), 0.0, 1.0)
        for row, value in zip(rows, scaled):
            row[f"{feature}_scaled"] = float(value)


def make_constant_velocity_motion_model(max_lost: int, prob_not_assign: float, accuracy: float) -> MotionModel:
    # K00 frame spacing is long enough for real motion, but the dense detections
    # make aggressive velocity uncertainty produce many bad cross-links.
    return MotionModel(
        name="k00_2d_constant_velocity",
        measurements=3,
        states=6,
        A=np.asarray(
            [
                [1, 0, 0, 1, 0, 0],
                [0, 1, 0, 0, 1, 0],
                [0, 0, 1, 0, 0, 1],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            dtype=float,
        ),
        H=np.asarray([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0]], dtype=float),
        P=np.diag([9.0, 9.0, 1.0, 16.0, 16.0, 1.0]),
        G=np.asarray([[1.5, 1.5, 0.05, 3.0, 3.0, 0.05]], dtype=float),
        R=np.diag([3.0, 3.0, 1.0]),
        dt=1.0,
        accuracy=accuracy,
        max_lost=max_lost,
        prob_not_assign=prob_not_assign,
    )


def make_random_walk_motion_model(max_lost: int, prob_not_assign: float, accuracy: float) -> MotionModel:
    return MotionModel(
        name="k00_2d_random_walk",
        measurements=3,
        states=3,
        A=np.eye(3, dtype=float),
        H=np.eye(3, dtype=float),
        P=np.diag([25.0, 25.0, 1.0]),
        R=np.diag([4.0, 4.0, 1.0]),
        Q=np.diag([36.0, 36.0, 1.0]),
        dt=1.0,
        accuracy=accuracy,
        max_lost=max_lost,
        prob_not_assign=prob_not_assign,
    )


def make_motion_model(kind: str, max_lost: int, prob_not_assign: float, accuracy: float) -> MotionModel:
    if kind == "random-walk":
        return make_random_walk_motion_model(max_lost, prob_not_assign, accuracy)
    if kind == "constant-velocity":
        return make_constant_velocity_motion_model(max_lost, prob_not_assign, accuracy)
    raise ValueError(f"Unsupported motion model: {kind}")


def tracking_updates_from_mode(mode: str) -> list[str]:
    if mode == "motion":
        return ["motion"]
    if mode == "motion-visual":
        return ["motion", "visual"]
    if mode == "visual":
        return ["visual"]
    raise ValueError(f"Unsupported tracking mode: {mode}")


def tracking_features_from_args(mode: str, feature_set: str) -> list[str]:
    if "visual" not in tracking_updates_from_mode(mode):
        return []
    if feature_set == "size":
        return SIZE_TRACKING_FEATURES
    if feature_set == "shape":
        return TRACKING_FEATURES
    raise ValueError(f"Unsupported tracking feature set: {feature_set}")


def write_objects_csv(path: Path, rows: list[dict[str, float | int]]) -> None:
    fieldnames = unique_fieldnames([
        "t",
        "x",
        "y",
        "z",
        "mask_label",
        "area",
        "major_axis_length",
        "minor_axis_length",
        "equivalent_diameter",
        "eccentricity",
        "solidity",
        "extent",
        "circularity",
        "orientation",
        *RAW_FEATURES,
        *TRACKING_FEATURES,
    ])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fieldnames} for row in rows)


def write_tracks_csv(path: Path, tracks: list) -> tuple[int, int]:
    fieldnames = unique_fieldnames([
        "track_id",
        "parent",
        "root",
        "generation",
        "t",
        "x",
        "y",
        "z",
        "dummy",
        "mask_label",
        "area",
        "major_axis_length",
        "minor_axis_length",
        "equivalent_diameter",
        "eccentricity",
        "solidity",
        "extent",
        "circularity",
        "orientation",
        *RAW_FEATURES,
        *TRACKING_FEATURES,
    ])
    n_rows = 0
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for track in tracks:
            props = track.properties
            parent = "" if track.parent is None else int(track.parent)
            root = "" if track.root is None else int(track.root)
            for i in range(len(track)):
                row = {
                    "track_id": int(track.ID),
                    "parent": parent,
                    "root": root,
                    "generation": int(track.generation),
                    "t": int(track.t[i]),
                    "x": float(track.x[i]),
                    "y": float(track.y[i]),
                    "z": float(track.z[i]),
                    "dummy": bool(track.dummy[i]),
                }
                for key in fieldnames:
                    if key in row:
                        continue
                    values = props.get(key)
                    row[key] = "" if values is None else values[i]
                writer.writerow(row)
                n_rows += 1
    return len(tracks), n_rows


def main() -> None:
    args = parse_args()
    mask_tiff = args.mask_tiff.resolve()
    site_id = site_id_from_mask(mask_tiff)

    tracks_dir = args.out_dir / "tracks"
    objects_dir = args.out_dir / "objects"
    manifest_dir = args.out_dir / "manifest_rows"
    tracks_dir.mkdir(parents=True, exist_ok=True)
    objects_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir.mkdir(parents=True, exist_ok=True)

    tracks_csv = tracks_dir / f"{site_id}_btrack_tracks.csv"
    objects_csv = objects_dir / f"{site_id}_btrack_objects.csv"
    manifest_json = manifest_dir / f"{site_id}_btrack_manifest.json"
    if tracks_csv.exists() and manifest_json.exists() and not args.overwrite:
        print(f"Existing outputs found for {site_id}; use --overwrite to rerun.")
        return

    started = time.time()
    masks = read_mask_stack(mask_tiff)
    if args.start_frame < 0:
        raise ValueError("--start-frame must be non-negative")
    if args.start_frame:
        masks = masks[args.start_frame :]
    if args.max_frames is not None:
        if args.max_frames < 1:
            raise ValueError("--max-frames must be positive")
        masks = masks[: args.max_frames]
    rows = object_rows_from_masks(masks)
    add_scaled_features(rows)
    write_objects_csv(objects_csv, rows)

    if rows:
        objects = localizations_to_objects({key: np.asarray([row[key] for row in rows]) for key in rows[0]})
    else:
        objects = []

    n_tracks = 0
    n_track_rows = 0
    tracking_updates = tracking_updates_from_mode(args.tracking_mode)
    tracking_features = tracking_features_from_args(args.tracking_mode, args.tracking_feature_set)
    if objects:
        cfg = TrackerConfig(
            name=f"k00_cpsam_{args.tracking_mode}_tracking",
            motion_model=make_motion_model(args.motion_model, args.max_lost, args.prob_not_assign, args.accuracy),
            max_search_radius=args.max_search_radius,
            features=tracking_features,
            tracking_updates=tracking_updates,
            enable_optimisation=False,
            verbose=args.verbose,
            volume=((0, int(masks.shape[2])), (0, int(masks.shape[1])), (-1, 1)),
        )
        with btrack.BayesianTracker(verbose=args.verbose) as tracker:
            tracker.configure(cfg)
            tracker.append(objects)
            tracker.track(step_size=args.step_size, tracking_updates=tracking_updates)
            tracks = tracker.tracks
        n_tracks, n_track_rows = write_tracks_csv(tracks_csv, tracks)
    else:
        write_tracks_csv(tracks_csv, [])

    frame_counts = [int(np.unique(frame[frame > 0]).size) for frame in masks]
    summary = {
        "site_id": site_id,
        "mask_tiff": str(mask_tiff),
        "tracks_csv": str(tracks_csv.resolve()),
        "objects_csv": str(objects_csv.resolve()),
        "n_frames": int(masks.shape[0]),
        "height": int(masks.shape[1]),
        "width": int(masks.shape[2]),
        "total_objects": int(len(rows)),
        "n_tracks": int(n_tracks),
        "n_track_rows": int(n_track_rows),
        "frame_object_counts": frame_counts,
        "tracking_mode": args.tracking_mode,
        "tracking_feature_set": args.tracking_feature_set,
        "motion_model": args.motion_model,
        "tracking_updates": tracking_updates,
        "tracking_features": tracking_features,
        "max_search_radius": args.max_search_radius,
        "max_lost": args.max_lost,
        "prob_not_assign": args.prob_not_assign,
        "accuracy": args.accuracy,
        "start_frame": args.start_frame,
        "max_frames": args.max_frames,
        "global_optimization": False,
        "elapsed_seconds": round(time.time() - started, 3),
        "status": "ok",
    }
    manifest_json.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {tracks_csv}")
    print(f"Wrote {objects_csv}")
    print(f"Wrote {manifest_json}")


if __name__ == "__main__":
    main()
