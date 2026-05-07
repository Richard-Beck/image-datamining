#!/usr/bin/env python3
"""Collapse focused tracking candidate directories into experiment-level groups."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


FOCUS_TERMS = ("chemotaxis", "migration", "motility", "tracking")
MODALITY_TOKENS = (
    "brightfield",
    "dead",
    "deadbottom",
    "deadtop",
    "fluorescence",
    "green",
    "live",
    "livebottom",
    "livetop",
    "phase",
    "phasebottom",
    "phasecontrast",
    "phasetop",
    "red",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize focused tracking candidate directories by likely experiment root."
    )
    parser.add_argument(
        "input_tsv",
        nargs="?",
        default="data_inventory/focused_tracking_candidates.tsv",
        help="Focused tracking candidate TSV from profile_acquisition_candidates.py.",
    )
    parser.add_argument(
        "--output",
        default="data_inventory/tracking_experiment_groups.tsv",
        help="Output experiment-level TSV.",
    )
    return parser.parse_args()


def normalized(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def infer_experiment_root(path: str) -> str:
    parts = Path(path).parts
    lower_parts = [normalized(part) for part in parts]

    focus_indices = [
        index
        for index, part in enumerate(lower_parts)
        if any(term in part for term in FOCUS_TERMS)
    ]
    if not focus_indices:
        return str(Path(path).parent)

    last_focus_index = focus_indices[-1]
    if last_focus_index + 1 >= len(parts):
        return path

    next_part = lower_parts[last_focus_index + 1]
    if re.fullmatch(r"col\d+", next_part) or next_part in MODALITY_TOKENS:
        return str(Path(*parts[: last_focus_index + 1]))

    return str(Path(*parts[: last_focus_index + 2]))


def int_or_zero(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def maybe_int(value: str) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


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


def modality_for_path(path: str) -> str:
    parts = [part for part in Path(path).parts if part not in {"/", "share", "andor_lab"}]
    tokens = []
    for part in parts[-3:]:
        token = normalized(part)
        if token in MODALITY_TOKENS or any(modality in token for modality in MODALITY_TOKENS):
            tokens.append(part)
    return "; ".join(tokens[-2:])


def summarize_group(root: str, rows: list[dict[str, str]]) -> dict[str, Any]:
    min_spacings = [value for value in (maybe_int(row["min_time_spacing_minutes"]) for row in rows) if value is not None]
    median_spacings = [value for value in (maybe_int(row["median_time_spacing_minutes"]) for row in rows) if value is not None]
    timepoint_counts = [int_or_zero(row["distinct_timepoints"]) for row in rows]
    image_counts = [int_or_zero(row["image_like_count"]) for row in rows]
    dir_scores = [int_or_zero(row["score"]) for row in rows]
    raw_or_derived = [row["raw_or_derived"] for row in rows]
    statuses = [row["status"] for row in rows]

    modalities = [modality_for_path(row["path"]) for row in rows]
    best_rows = sorted(
        rows,
        key=lambda row: (
            maybe_int(row["min_time_spacing_minutes"]) is None,
            maybe_int(row["min_time_spacing_minutes"]) or 10**9,
            -int_or_zero(row["distinct_timepoints"]),
            -int_or_zero(row["image_like_count"]),
            row["path"],
        ),
    )

    score = 0
    if min_spacings:
        score += max(0, 6 - min(min_spacings) // 60)
    score += min(max(timepoint_counts or [0]), 10)
    score += min(len(rows), 10)
    if any("phase" in normalized(value) or "brightfield" in normalized(value) for value in modalities):
        score += 3
    if any("chemotaxis" in row["tracking_focus_matches"] for row in rows):
        score += 2
    if any(status == "timeout_blackhole" for status in statuses):
        score -= 1

    return {
        "group_score": score,
        "experiment_root": root,
        "candidate_dirs": len(rows),
        "total_image_like_count": sum(image_counts),
        "best_dir_score": max(dir_scores or [0]),
        "min_time_spacing_minutes": "" if not min_spacings else min(min_spacings),
        "best_median_time_spacing_minutes": "" if not median_spacings else min(median_spacings),
        "max_distinct_timepoints_sampled": max(timepoint_counts or [0]),
        "modalities_or_views": unique_join(modalities),
        "cell_line_hits": unique_join([row["cell_line_hits"] for row in rows]),
        "tracking_focus_matches": unique_join([row["tracking_focus_matches"] for row in rows]),
        "raw_or_derived_values": unique_join(raw_or_derived),
        "statuses": unique_join(statuses),
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
            row["min_time_spacing_minutes"] == "",
            int(row["min_time_spacing_minutes"] or 10**9),
            -int(row["candidate_dirs"]),
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
        "min_time_spacing_minutes",
        "best_median_time_spacing_minutes",
        "max_distinct_timepoints_sampled",
        "modalities_or_views",
        "cell_line_hits",
        "tracking_focus_matches",
        "raw_or_derived_values",
        "statuses",
        "best_candidate_dir",
        "example_dirs",
        "example_filenames",
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)

    print(f"wrote {len(summaries)} experiment groups to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
