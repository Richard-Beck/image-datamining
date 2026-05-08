#!/usr/bin/env python3
"""Compute per-object nearest-object distances for one K00 CPSAM mask stack."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import tifffile as tf

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from image_datamining.mask_distances import nearest_mask_neighbors


DEFAULT_OUT_DIR = Path("analyses/K00_GemcitabineExposure_033023/cpsam_full_stacks/nearest_distances")
SITE_RE = re.compile(r"^(?P<site_id>.+?)(?:_cpsam_masks)?\.tiff?$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute exact per-object nearest-object mask-pixel-center distances "
            "for every frame in one CPSAM labeled mask stack."
        )
    )
    parser.add_argument("--mask-tiff", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--connectivity", type=int, choices=(4, 8), default=8)
    parser.add_argument("--start-frame", type=int, default=0, help="Zero-based first frame to process.")
    parser.add_argument("--max-frames", type=int, default=None, help="Process at most N frames after --start-frame.")
    parser.add_argument("--overwrite", action="store_true")
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


def frame_rows(site_id: str, frame: int, mask: np.ndarray, connectivity: int) -> list[dict[str, object]]:
    graph = nearest_mask_neighbors(mask, connectivity=connectivity)
    if graph.labels.size == 0:
        return []

    counts = np.bincount(mask.ravel())
    rows: list[dict[str, object]] = []
    n_objects = int(graph.labels.size)
    for label in graph.labels:
        label = int(label)
        dist = graph.nearest_other_dist[label]
        has_neighbor = bool(np.isfinite(dist))
        center_distance = float(dist) if has_neighbor else ""
        rows.append(
            {
                "site_id": site_id,
                "frame": int(frame),
                "cpsam_label": label,
                "cpsam_area_px": int(counts[label]) if label < counts.size else 0,
                "nearest_cpsam_label": int(graph.nearest_other_label[label]) if has_neighbor else "",
                "nearest_mask_center_distance_px": center_distance,
                "nearest_empty_gap_px": max(float(dist) - 1.0, 0.0) if has_neighbor else "",
                "nearest_y0": int(graph.source_y[label]) if has_neighbor else "",
                "nearest_x0": int(graph.source_x[label]) if has_neighbor else "",
                "nearest_y1": int(graph.neighbor_y[label]) if has_neighbor else "",
                "nearest_x1": int(graph.neighbor_x[label]) if has_neighbor else "",
                "n_objects_frame": n_objects,
                "connectivity": int(connectivity),
                "distance_definition": "mask_pixel_center_to_mask_pixel_center",
            }
        )
    return rows


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "site_id",
        "frame",
        "cpsam_label",
        "cpsam_area_px",
        "nearest_cpsam_label",
        "nearest_mask_center_distance_px",
        "nearest_empty_gap_px",
        "nearest_y0",
        "nearest_x0",
        "nearest_y1",
        "nearest_x1",
        "n_objects_frame",
        "connectivity",
        "distance_definition",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    mask_tiff = args.mask_tiff.resolve()
    site_id = site_id_from_mask(mask_tiff)
    out_dir = args.out_dir
    table_dir = out_dir / "object_tables"
    summary_dir = out_dir / "summaries"
    out_tsv = table_dir / f"{site_id}_nearest_object_distances.tsv"
    out_summary = summary_dir / f"{site_id}_nearest_object_distances_summary.json"

    if out_tsv.exists() and out_summary.exists() and not args.overwrite:
        print(f"Existing outputs found for {site_id}; use --overwrite to rerun.")
        return

    started = time.time()
    masks = read_mask_stack(mask_tiff)
    start_frame = max(args.start_frame, 0)
    end_frame = masks.shape[0] if args.max_frames is None else min(masks.shape[0], start_frame + args.max_frames)
    if start_frame >= end_frame:
        raise ValueError(f"No frames selected from stack with {masks.shape[0]} frames")

    print(f"Mask TIFF: {mask_tiff}")
    print(f"Site ID: {site_id}")
    print(f"Stack shape TYX: {masks.shape}")
    print(f"Frames: {start_frame}..{end_frame - 1}")
    print(f"Output TSV: {out_tsv}")

    all_rows: list[dict[str, object]] = []
    frame_seconds: list[float] = []
    for frame in range(start_frame, end_frame):
        frame_started = time.time()
        rows = frame_rows(site_id, frame, masks[frame], args.connectivity)
        frame_seconds.append(time.time() - frame_started)
        all_rows.extend(rows)
        print(
            f"frame={frame} objects={len(rows)} seconds={frame_seconds[-1]:.3f}",
            flush=True,
        )

    write_rows(out_tsv, all_rows)
    elapsed = time.time() - started
    distances = [
        float(row["nearest_mask_center_distance_px"])
        for row in all_rows
        if row["nearest_mask_center_distance_px"] != ""
    ]
    summary = {
        "site_id": site_id,
        "mask_tiff": str(mask_tiff),
        "output_tsv": str(out_tsv),
        "n_frames_total": int(masks.shape[0]),
        "start_frame": int(start_frame),
        "end_frame_exclusive": int(end_frame),
        "height": int(masks.shape[1]),
        "width": int(masks.shape[2]),
        "n_object_rows": int(len(all_rows)),
        "n_missing_neighbor": int(len(all_rows) - len(distances)),
        "connectivity": int(args.connectivity),
        "elapsed_seconds": round(elapsed, 3),
        "median_frame_seconds": round(float(np.median(frame_seconds)), 3) if frame_seconds else None,
        "max_frame_seconds": round(float(np.max(frame_seconds)), 3) if frame_seconds else None,
        "nearest_distance_px_median": round(float(np.median(distances)), 6) if distances else None,
        "nearest_distance_px_p90": round(float(np.quantile(distances, 0.90)), 6) if distances else None,
        "nearest_distance_px_p99": round(float(np.quantile(distances, 0.99)), 6) if distances else None,
        "nearest_distance_px_max": round(float(np.max(distances)), 6) if distances else None,
    }
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    out_summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
