#!/usr/bin/env python3
"""Collapse z-stack candidate directories into likely experiment-level groups."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


FIELD_OR_SITE_PATTERNS = (
    re.compile(r"fof\d+", re.IGNORECASE),
    re.compile(r"r\d{1,2}c\d{1,2}f\d{1,3}(?:p\d{1,4})?", re.IGNORECASE),
    re.compile(r"sk\d{1,3}", re.IGNORECASE),
)
METADATA_OR_DERIVED_TOKENS = {
    "digitalprojections",
    "metadata",
    "output",
    "processedframes",
}
IMAGE_CONTAINER_TOKENS = {
    "images",
    "imagesorganized",
    "raw",
    "rawimages",
    "raw_images",
}
VOLUME_EXPERIMENT_TOKENS = {
    "3dbrightfieldimages",
    "sum159volumetricanalysis",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize z-stack candidate directories by likely experiment root."
    )
    parser.add_argument(
        "input_tsv",
        nargs="?",
        default="data_inventory/z_stack_candidates.tsv",
        help="Z-stack candidate TSV from profile_acquisition_candidates.py.",
    )
    parser.add_argument(
        "--output",
        default="data_inventory/z_stack_experiment_groups.tsv",
        help="Output experiment-level TSV.",
    )
    return parser.parse_args()


def normalized(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def strip_known_nonexperiment_suffixes(parts: list[str]) -> list[str]:
    while parts and normalized(parts[-1]) in METADATA_OR_DERIVED_TOKENS:
        parts = parts[:-1]
    return parts


def matches_any_pattern(text: str, patterns: tuple[re.Pattern[str], ...]) -> bool:
    token = normalized(text)
    return any(pattern.search(token) for pattern in patterns)


def infer_experiment_root(path: str) -> str:
    parts = list(Path(path).parts)
    parts = strip_known_nonexperiment_suffixes(parts)
    lower_parts = [normalized(part) for part in parts]

    measurement_indices = [
        index
        for index, part in enumerate(lower_parts)
        if "measurement" in part or re.search(r"t\d{4}\d{2}\d{2}", part)
    ]
    if measurement_indices:
        return str(Path(*parts[: measurement_indices[-1] + 1]))

    volume_indices = [
        index
        for index, part in enumerate(lower_parts)
        if part in VOLUME_EXPERIMENT_TOKENS
    ]
    if volume_indices:
        index = volume_indices[-1]
        if index + 1 < len(parts):
            return str(Path(*parts[: index + 2]))
        return str(Path(*parts[: index + 1]))

    image_container_indices = [
        index
        for index, part in enumerate(lower_parts)
        if part in IMAGE_CONTAINER_TOKENS
    ]
    if image_container_indices:
        index = image_container_indices[-1]
        if index + 1 >= len(parts):
            return str(Path(*parts[: index + 1]))
        if matches_any_pattern(parts[-1], FIELD_OR_SITE_PATTERNS):
            return str(Path(*parts[: index + 1]))

    if parts and matches_any_pattern(parts[-1], FIELD_OR_SITE_PATTERNS):
        return str(Path(*parts[:-1]))

    return path


def int_or_zero(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def unique_join(values: list[str], limit: int = 16) -> str:
    seen: list[str] = []
    for value in values:
        if not value:
            continue
        for item in [part.strip() for part in value.split(";")]:
            if item and item not in seen:
                seen.append(item)
    rendered = seen[:limit]
    if len(seen) > limit:
        rendered.append(f"...(+{len(seen) - limit})")
    return "; ".join(rendered)


def unique_items(values: list[str]) -> list[str]:
    seen: list[str] = []
    for value in values:
        if not value:
            continue
        for item in [part.strip() for part in value.split(";")]:
            if item and item not in seen:
                seen.append(item)
    return seen


def suffix_after_root(root: str, path: str) -> str:
    try:
        relpath = Path(path).relative_to(root)
    except ValueError:
        return ""
    rendered = str(relpath)
    return "" if rendered == "." else rendered


def summarize_group(root: str, rows: list[dict[str, str]]) -> dict[str, Any]:
    dir_scores = [int_or_zero(row["score"]) for row in rows]
    image_counts = [int_or_zero(row["image_like_count"]) for row in rows]
    z_counts = [int_or_zero(row["distinct_z"]) for row in rows]
    plane_counts = [int_or_zero(row["distinct_planes"]) for row in rows]
    channel_counts = [
        int_or_zero(row.get("distinct_channels", ""))
        or len([part for part in row["channels"].split(";") if part.strip()])
        for row in rows
    ]
    group_channels = unique_items([row["channels"] for row in rows])
    multi_channel_rows = [row for row in rows if int_or_zero(row.get("has_multiple_channels", "")) >= 1]
    raw_or_derived = [row["raw_or_derived"] for row in rows]
    statuses = [row["status"] for row in rows]

    best_rows = sorted(
        rows,
        key=lambda row: (
            -int_or_zero(row["score"]),
            -max(int_or_zero(row["distinct_z"]), int_or_zero(row["distinct_planes"])),
            row["raw_or_derived"] == "derived",
            -int_or_zero(row["image_like_count"]),
            row["path"],
        ),
    )
    root_norm = normalized(root)

    score = 0
    score += max(dir_scores or [0])
    score += min(len(rows), 12)
    score += min(max(z_counts or [0]), 10)
    score += min(max(plane_counts or [0]), 10)
    if any(value == "raw_hint" for value in raw_or_derived):
        score += 3
    if any(value == "derived" for value in raw_or_derived):
        score -= 3
    if "volumetricanalysis" in root_norm or "3dbrightfield" in root_norm:
        score += 4
    if len(group_channels) >= 2:
        score += 2
    if any(status == "timeout_blackhole" for status in statuses):
        score -= 1

    return {
        "group_score": score,
        "experiment_root": root,
        "candidate_dirs": len(rows),
        "total_image_like_count": sum(image_counts),
        "best_dir_score": max(dir_scores or [0]),
        "max_distinct_z_sampled": max(z_counts or [0]),
        "max_distinct_planes_sampled": max(plane_counts or [0]),
        "max_channels_sampled": max(channel_counts or [0]),
        "group_distinct_channels": len(group_channels),
        "has_multiple_channels": "1" if len(group_channels) >= 2 else "0",
        "multi_channel_candidate_dirs": len(multi_channel_rows),
        "cell_line_hits": unique_join([row["cell_line_hits"] for row in rows]),
        "channels": unique_join([row["channels"] for row in rows]),
        "channel_z_or_plane_coverage": unique_join(
            [row.get("channel_z_or_plane_coverage", "") for row in rows]
        ),
        "raw_or_derived_values": unique_join(raw_or_derived),
        "statuses": unique_join(statuses),
        "condition_or_field_suffixes": unique_join(
            [suffix_after_root(root, row["path"]) for row in best_rows], limit=20
        ),
        "best_candidate_dir": best_rows[0]["path"],
        "example_dirs": "; ".join(row["path"] for row in best_rows[:5]),
        "example_filenames": unique_join([row["examples"] for row in best_rows[:2]], limit=12),
    }


def main() -> int:
    args = parse_args()
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)

    with Path(args.input_tsv).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            groups[infer_experiment_root(row["path"])].append(row)

    summaries = [summarize_group(root, rows) for root, rows in groups.items()]
    summaries.sort(
        key=lambda row: (
            -int(row["group_score"]),
            -int(row["candidate_dirs"]),
            -int(row["total_image_like_count"]),
            row["experiment_root"],
        )
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "group_score",
        "experiment_root",
        "candidate_dirs",
        "total_image_like_count",
        "best_dir_score",
        "max_distinct_z_sampled",
        "max_distinct_planes_sampled",
        "max_channels_sampled",
        "group_distinct_channels",
        "has_multiple_channels",
        "multi_channel_candidate_dirs",
        "cell_line_hits",
        "channels",
        "channel_z_or_plane_coverage",
        "raw_or_derived_values",
        "statuses",
        "condition_or_field_suffixes",
        "best_candidate_dir",
        "example_dirs",
        "example_filenames",
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    print(f"wrote {len(summaries)} z-stack experiment groups to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
