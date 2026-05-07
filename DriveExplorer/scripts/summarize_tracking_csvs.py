#!/usr/bin/env python3
"""Summarize existing tracking CSV exports for basic QC."""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


DEFAULT_THRESHOLDS = (3, 5, 10, 25)
KEY_COLUMNS = ("trackId", "Object_Center_0", "Object_Center_1", "Object_Area_0")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize tracking CSV files.")
    parser.add_argument(
        "root",
        nargs="?",
        default="/share/andor_lab/Jackson/FUCCI_TimepointAnalysis/QI_Core_Analysis_Results/Tracking_Data/Tiffs",
        help="Root containing per-cell-line csv folders.",
    )
    parser.add_argument(
        "--output",
        default="data_inventory/fucci_tracking_csv_qc.tsv",
        help="Output TSV path.",
    )
    parser.add_argument(
        "--thresholds",
        default=",".join(str(value) for value in DEFAULT_THRESHOLDS),
        help="Comma-separated track-length thresholds in frames.",
    )
    return parser.parse_args()


def parse_thresholds(raw_thresholds: str) -> list[int]:
    thresholds = sorted({int(value) for value in raw_thresholds.split(",") if value.strip()})
    if not thresholds:
        raise ValueError("At least one threshold is required.")
    return thresholds


def csv_paths(root: Path) -> list[Path]:
    return sorted(root.glob("*/csv/*.csv"))


def infer_cell_line(path: Path, root: Path) -> str:
    relative = path.relative_to(root)
    return relative.parts[0]


def infer_field_id(path: Path) -> str:
    name = path.name
    marker = "_CSV-Table_"
    if marker in name:
        return name.split(marker, 1)[0]
    return path.stem


def median_int(values: list[int]) -> str:
    if not values:
        return ""
    return f"{statistics.median(values):.1f}"


def summarize_csv(path: Path, root: Path, thresholds: list[int]) -> dict[str, Any]:
    row_count = 0
    frames: set[int] = set()
    track_frames: dict[str, set[int]] = defaultdict(set)
    lineage_ids: set[str] = set()
    parent_track_links: set[tuple[str, str]] = set()
    merger_labels: set[str] = set()
    objects_per_frame: Counter[int] = Counter()
    missing_counts: Counter[str] = Counter()
    columns_seen: list[str] = []

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        columns_seen = reader.fieldnames or []
        for row in reader:
            row_count += 1

            frame_raw = row.get("frame", "")
            try:
                frame = int(frame_raw)
                frames.add(frame)
                objects_per_frame[frame] += 1
            except ValueError:
                frame = None

            track_id = row.get("trackId", "")
            if track_id and frame is not None:
                track_frames[track_id].add(frame)

            lineage_id = row.get("lineageId", "")
            if lineage_id:
                lineage_ids.add(lineage_id)

            parent_track_id = row.get("parentTrackId", "")
            if track_id and parent_track_id and parent_track_id not in {"0", "-1"}:
                parent_track_links.add((track_id, parent_track_id))

            merger_label_id = row.get("mergerLabelId", "")
            if merger_label_id and merger_label_id not in {"0", "-1"}:
                merger_labels.add(merger_label_id)

            for column in KEY_COLUMNS:
                if column not in row or row[column] in {"", "nan", "NaN"}:
                    missing_counts[column] += 1

    track_lengths = [len(frame_set) for frame_set in track_frames.values()]
    object_counts = list(objects_per_frame.values())
    frame_min = min(frames) if frames else ""
    frame_max = max(frames) if frames else ""

    summary: dict[str, Any] = {
        "cell_line": infer_cell_line(path, root),
        "field_id": infer_field_id(path),
        "csv_path": str(path),
        "n_rows": row_count,
        "n_columns": len(columns_seen),
        "n_frames": len(frames),
        "min_frame": frame_min,
        "max_frame": frame_max,
        "n_tracks": len(track_frames),
        "median_track_length_frames": median_int(track_lengths),
        "max_track_length_frames": max(track_lengths) if track_lengths else "",
        "median_objects_per_frame": median_int(object_counts),
        "min_objects_per_frame": min(object_counts) if object_counts else "",
        "max_objects_per_frame": max(object_counts) if object_counts else "",
        "n_lineages": len(lineage_ids),
        "n_parent_track_links": len(parent_track_links),
        "n_merger_labels": len(merger_labels),
    }

    for threshold in thresholds:
        summary[f"n_tracks_ge_{threshold}_frames"] = sum(1 for length in track_lengths if length >= threshold)

    for column in KEY_COLUMNS:
        summary[f"missing_{column}"] = missing_counts[column]

    return summary


def write_rows(output: Path, rows: list[dict[str, Any]], thresholds: list[int]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "cell_line",
        "field_id",
        "csv_path",
        "n_rows",
        "n_columns",
        "n_frames",
        "min_frame",
        "max_frame",
        "n_tracks",
        "median_track_length_frames",
        "max_track_length_frames",
        *[f"n_tracks_ge_{threshold}_frames" for threshold in thresholds],
        "median_objects_per_frame",
        "min_objects_per_frame",
        "max_objects_per_frame",
        "n_lineages",
        "n_parent_track_links",
        "n_merger_labels",
        *[f"missing_{column}" for column in KEY_COLUMNS],
    ]

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    root = Path(args.root)
    thresholds = parse_thresholds(args.thresholds)
    paths = csv_paths(root)
    if not paths:
        raise SystemExit(f"No CSV files found under {root}")

    rows = [summarize_csv(path, root, thresholds) for path in paths]
    rows.sort(key=lambda row: (row["cell_line"], row["field_id"]))
    write_rows(Path(args.output), rows, thresholds)
    print(f"wrote {len(rows)} CSV summaries to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
