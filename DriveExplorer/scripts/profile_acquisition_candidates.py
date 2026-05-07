#!/usr/bin/env python3
"""Rank inventory directories for z-stack volume and tracking follow-up.

This script only uses directory paths and sampled filenames from the existing
tree-walk JSONL. It does not read image files, image headers, or pixels.
"""

from __future__ import annotations

import argparse
import csv
import difflib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any


IMAGE_EXTENSIONS = {
    ".bmp",
    ".czi",
    ".gif",
    ".jpg",
    ".jpeg",
    ".lif",
    ".png",
    ".svs",
    ".tif",
    ".tiff",
}

RAW_HINTS = {
    "brightfield",
    "fluorescence",
    "images",
    "phase",
    "phasecontrast",
    "raw",
    "raw_images",
}

DERIVED_HINTS = {
    "bounding_boxes",
    "cell_segmentation_results",
    "cellpose",
    "cp_masks",
    "flattened",
    "flattenedimages",
    "log_transforms",
    "masks",
    "objectimages",
    "processed",
    "processed_images",
    "results",
    "segmentation",
}

TRACKING_FOCUS_TERMS = ("chemotaxis", "migration", "motility", "tracking")
EXCLUDED_PATH_TERMS = ("myco",)

CELL_LINE_PATTERN = re.compile(
    r"(?<![A-Za-z0-9])(?:"
    r"SNU[-_]?\d+[A-Za-z0-9]*|"
    r"SUM[-_]?\d+[A-Za-z0-9]*|"
    r"NCI[-_]?[A-Z]{1,3}\d+[A-Za-z0-9]*|"
    r"NUGC[-_]?\d+|"
    r"MCF[-_]?\d+|"
    r"MDA[-_]?[A-Z]{2}[-_]?\d+|"
    r"HEK[-_]?\d+|"
    r"U2OS|A549|HeLa|HCT[-_]?\d+"
    r")(?![A-Za-z0-9])",
    re.IGNORECASE,
)

EXPLICIT_Z_PATTERN = re.compile(r"(?<![A-Za-z0-9])z(?P<z>\d{1,4})(?![A-Za-z0-9])", re.IGNORECASE)
ZEISS_TIME_PATTERN = re.compile(r"(?<![A-Za-z0-9])t(?P<t>\d{1,5})(?![A-Za-z0-9])", re.IGNORECASE)
ELAPSED_PATTERN = re.compile(r"(?P<days>\d{1,3})d(?P<hours>\d{1,2})h(?P<minutes>\d{1,2})m", re.IGNORECASE)
INC_DATE_PATTERN = re.compile(
    r"(?P<year>\d{4})y(?P<month>\d{1,2})m(?P<day>\d{1,2})d_"
    r"(?P<hour>\d{1,2})h(?P<minute>\d{1,2})m",
    re.IGNORECASE,
)
OPERAPHENIX_PATTERN = re.compile(
    r"r(?P<row>\d{1,2})c(?P<column>\d{1,2})f(?P<field>\d{1,3})p(?P<plane>\d{1,4})"
    r"(?:-ch(?P<channel>\d{1,2}))?",
    re.IGNORECASE,
)
WELL_SITE_PATTERN = re.compile(r"(?<![A-Za-z0-9])(?P<well>[A-H]\d{1,2})_(?P<site>\d{1,3})(?![A-Za-z0-9])")
CHANNEL_PATTERN = re.compile(r"(?<![A-Za-z0-9])ch(?P<channel>\d{1,2})(?![A-Za-z0-9])", re.IGNORECASE)
OME_SIZE_PATTERN = re.compile(r'\b(Size[ZTCXY])="([^"]+)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank directories for z-stack volume and longitudinal tracking follow-up."
    )
    parser.add_argument("input_jsonl", help="Input JSONL from walk_directory_tree.py.")
    parser.add_argument(
        "--z-output",
        default="data_inventory/z_stack_candidates.tsv",
        help="Output TSV for z-stack candidates.",
    )
    parser.add_argument(
        "--tracking-output",
        default="data_inventory/tracking_candidates.tsv",
        help="Output TSV for tracking candidates.",
    )
    parser.add_argument(
        "--focused-tracking-output",
        default="data_inventory/focused_tracking_candidates.tsv",
        help="Output TSV for tracking candidates with path matches to tracking-related terms.",
    )
    parser.add_argument(
        "--min-score",
        type=int,
        default=4,
        help="Only emit rows with score at least this value.",
    )
    parser.add_argument(
        "--example-limit",
        type=int,
        default=8,
        help="Maximum example filenames stored per row.",
    )
    parser.add_argument(
        "--metadata-stack-scan-limit",
        type=int,
        default=20,
        help=(
            "Maximum sampled direct-child TIFF files per directory to inspect for stack metadata. "
            "The scan stops early after the first stack hit. Set to 0 to disable TIFF metadata reads."
        ),
    )
    return parser.parse_args()


def normalize_token(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def image_like_count(record: dict[str, Any]) -> int:
    ext_counts = record.get("extension_counts") or record.get("sampleExtensions") or {}
    return sum(int(ext_counts.get(ext, 0)) for ext in IMAGE_EXTENSIONS)


def status_label(path: str, names: list[str]) -> str:
    haystack = " ".join([path, *names])
    normalized = {normalize_token(part) for part in re.split(r"[/\s;]+", haystack) if part}
    if normalized & DERIVED_HINTS:
        return "derived"
    if normalized & RAW_HINTS:
        return "raw_hint"
    return "unknown"


def collect_path_context(path: str) -> tuple[list[str], list[str]]:
    cell_lines = sorted({match.group(0).replace("_", "-") for match in CELL_LINE_PATTERN.finditer(path)})
    parts = [part for part in Path(path).parts if part not in {"/", "share", "andor_lab"}]
    context = parts[-5:]
    return cell_lines, context


def path_focus_match(path: str) -> tuple[str, str, str]:
    normalized_path = path.lower()
    tokens = [token for token in re.split(r"[^a-z0-9]+", normalized_path) if token]
    matches: list[str] = []
    best_term = ""
    best_token = ""
    best_ratio = 0.0

    for term in TRACKING_FOCUS_TERMS:
        if term in normalized_path:
            matches.append(f"{term}:exact")
            if best_ratio < 1.0:
                best_term = term
                best_token = term
                best_ratio = 1.0
            continue

        close_tokens = difflib.get_close_matches(term, tokens, n=1, cutoff=0.78)
        if close_tokens:
            token = close_tokens[0]
            ratio = difflib.SequenceMatcher(a=term, b=token).ratio()
            matches.append(f"{term}:soft:{token}:{ratio:.2f}")
            if ratio > best_ratio:
                best_term = term
                best_token = token
                best_ratio = ratio

    return "; ".join(matches), f"{best_ratio:.2f}" if best_ratio else "", f"{best_term}:{best_token}" if best_term else ""


def is_excluded_path(path: str) -> bool:
    normalized_path = path.lower()
    return any(term in normalized_path for term in EXCLUDED_PATH_TERMS)


def profile_names(names: list[str]) -> dict[str, Any]:
    z_values: set[str] = set()
    zeiss_t_values: set[str] = set()
    elapsed_minutes: set[int] = set()
    inc_datetimes: set[str] = set()
    absolute_minutes: set[int] = set()
    plane_values: set[str] = set()
    wells: set[str] = set()
    fields: set[str] = set()
    channels: set[str] = set()
    channel_z_values: dict[str, set[str]] = defaultdict(set)
    channel_plane_values: dict[str, set[str]] = defaultdict(set)
    z_examples: list[str] = []
    t_examples: list[str] = []
    image_examples: list[str] = []
    extension_counts: Counter[str] = Counter()

    for name in names:
        extension = Path(name).suffix.lower()
        extension_counts[extension or "[no extension]"] += 1
        is_image = extension in IMAGE_EXTENSIONS
        if is_image:
            image_examples.append(name)

        z_match = EXPLICIT_Z_PATTERN.search(name)
        if z_match:
            z_values.add(z_match.group("z"))
            z_examples.append(name)

        time_match = ZEISS_TIME_PATTERN.search(name)
        if time_match:
            zeiss_t_values.add(time_match.group("t"))
            t_examples.append(name)

        elapsed_match = ELAPSED_PATTERN.search(name)
        if elapsed_match:
            minutes = (
                int(elapsed_match.group("days")) * 24 * 60
                + int(elapsed_match.group("hours")) * 60
                + int(elapsed_match.group("minutes"))
            )
            elapsed_minutes.add(minutes)
            t_examples.append(name)

        inc_match = INC_DATE_PATTERN.search(name)
        if inc_match:
            acquired_at = "{year}-{month:0>2}-{day:0>2} {hour:0>2}:{minute:0>2}".format(
                **inc_match.groupdict()
            )
            inc_datetimes.add(acquired_at)
            absolute_minutes.add(
                int(datetime.strptime(acquired_at, "%Y-%m-%d %H:%M").timestamp() // 60)
            )
            t_examples.append(name)

        op_match = OPERAPHENIX_PATTERN.search(name)
        if op_match:
            fields.add(op_match.group("field"))
            plane_values.add(op_match.group("plane"))
            if op_match.group("channel"):
                channel = op_match.group("channel")
                channels.add(channel)
                channel_plane_values[channel].add(op_match.group("plane"))

        well_match = WELL_SITE_PATTERN.search(name)
        if well_match:
            wells.add(well_match.group("well").upper())
            fields.add(well_match.group("site"))

        channel_matches = list(CHANNEL_PATTERN.finditer(name))
        for channel_match in channel_matches:
            channel = channel_match.group("channel")
            channels.add(channel)
            if z_match:
                channel_z_values[channel].add(z_match.group("z"))

    timepoint_count = len(zeiss_t_values) + len(elapsed_minutes) + len(inc_datetimes)
    return {
        "extensions": extension_counts,
        "image_examples": image_examples,
        "z_values": z_values,
        "zeiss_t_values": zeiss_t_values,
        "elapsed_minutes": elapsed_minutes,
        "inc_datetimes": inc_datetimes,
        "absolute_minutes": absolute_minutes,
        "plane_values": plane_values,
        "wells": wells,
        "fields": fields,
        "channels": channels,
        "channel_z_values": channel_z_values,
        "channel_plane_values": channel_plane_values,
        "timepoint_count": timepoint_count,
        "z_examples": z_examples,
        "t_examples": t_examples,
    }


def inspect_tiff_stacks(record: dict[str, Any], names: list[str], limit: int) -> dict[str, Any]:
    """Read TIFF headers for a small sample of direct-child files.

    This intentionally reads metadata only, not image pixels. It catches
    OME-TIFF and generic multi-page TIFF stacks that do not encode z/t in
    filenames.
    """
    empty = {
        "stack_tiff_count": 0,
        "stack_tiff_examples": [],
        "stack_tiff_kind": set(),
        "stack_tiff_ome_sizes": set(),
        "stack_tiff_shapes": set(),
        "stack_tiff_pages": set(),
        "stack_tiff_errors": [],
    }
    if limit <= 0:
        return empty

    try:
        import tifffile
    except ImportError:
        return empty

    directory = Path(record["path"])
    inspected = 0
    stack_examples: list[str] = []
    kind_values: set[str] = set()
    ome_size_values: set[str] = set()
    shape_values: set[str] = set()
    page_values: set[int] = set()
    errors: list[str] = []

    for name in names:
        if inspected >= limit:
            break
        if Path(name).suffix.lower() not in {".tif", ".tiff"}:
            continue

        file_path = directory / name
        if not file_path.is_file():
            continue

        inspected += 1
        try:
            with tifffile.TiffFile(file_path) as tif:
                page0 = tif.pages[0]
                description = page0.description or ""
                ome_sizes = {key: int(value) for key, value in OME_SIZE_PATTERN.findall(description)}
                page0_shape = tuple(page0.shape)
                page0_dtype = str(page0.dtype)
                page_count = 0
                stack_kind = ""

                if ome_sizes:
                    size_z = ome_sizes.get("SizeZ", 1)
                    size_t = ome_sizes.get("SizeT", 1)
                    size_c = ome_sizes.get("SizeC", 1)
                    if max(size_z, size_t) >= 2:
                        stack_kind = "ome_z_or_t_stack"
                    elif size_c >= 2:
                        stack_kind = "ome_multichannel_single_plane"
                else:
                    page_count = len(tif.pages)
                    if page_count >= 10:
                        stack_kind = "ambiguous_multipage_tiff"
        except Exception as error:
            if len(errors) < 3:
                errors.append(f"{name}:{type(error).__name__}")
            continue

        if not stack_kind:
            continue

        stack_examples.append(name)
        kind_values.add(stack_kind)
        if ome_sizes:
            ome_size_values.add(
                ",".join(f"{key}={ome_sizes[key]}" for key in ("SizeC", "SizeT", "SizeZ", "SizeY", "SizeX") if key in ome_sizes)
            )
        shape_values.add(f"page0:{'x'.join(str(value) for value in page0_shape)}:{page0_dtype}")
        if page_count:
            page_values.add(page_count)
        break

    return {
        "stack_tiff_count": len(stack_examples),
        "stack_tiff_examples": stack_examples,
        "stack_tiff_kind": kind_values,
        "stack_tiff_ome_sizes": ome_size_values,
        "stack_tiff_shapes": shape_values,
        "stack_tiff_pages": page_values,
        "stack_tiff_errors": errors,
    }


def spacing_summary(profile: dict[str, Any]) -> tuple[str, str, str]:
    if len(profile["elapsed_minutes"]) >= 2:
        values = sorted(profile["elapsed_minutes"])
        unit = "elapsed_minutes"
    elif len(profile["absolute_minutes"]) >= 2:
        values = sorted(profile["absolute_minutes"])
        unit = "absolute_minutes"
    else:
        return "", "", ""

    intervals = [after - before for before, after in zip(values, values[1:]) if after > before]
    if not intervals:
        return "", "", unit

    positive_intervals = sorted(intervals)
    median = positive_intervals[len(positive_intervals) // 2]
    return str(positive_intervals[0]), str(median), unit


def score_z(
    profile: dict[str, Any],
    raw_status: str,
    image_count: int,
    cell_lines: list[str],
    stack_profile: dict[str, Any],
) -> int:
    score = 0
    if image_count:
        score += 1
    if len(profile["z_values"]) >= 2:
        score += 5
    elif len(profile["z_values"]) == 1:
        score += 2
    if len(profile["plane_values"]) >= 5:
        score += 2
    elif len(profile["plane_values"]) >= 2:
        score += 1
    if profile["channels"]:
        score += 1
    if len(profile["channels"]) >= 2:
        score += 1
    if profile["wells"] or profile["fields"]:
        score += 1
    if cell_lines:
        score += 1
    if stack_profile["stack_tiff_count"]:
        score += 4
    if raw_status == "derived":
        score -= 3
    return score


def score_tracking(
    profile: dict[str, Any],
    raw_status: str,
    image_count: int,
    cell_lines: list[str],
    stack_profile: dict[str, Any],
) -> int:
    score = 0
    if image_count:
        score += 1
    if profile["timepoint_count"] >= 5:
        score += 5
    elif profile["timepoint_count"] >= 2:
        score += 3
    elif profile["timepoint_count"] == 1:
        score += 1
    if len(profile["wells"]) >= 2:
        score += 2
    elif profile["wells"]:
        score += 1
    if len(profile["fields"]) >= 2:
        score += 1
    if profile["channels"]:
        score += 1
    if cell_lines:
        score += 1
    if stack_profile["stack_tiff_count"]:
        score += 2
    if raw_status == "derived":
        score -= 2
    return score


def compact_values(values: set[Any], limit: int = 12) -> str:
    ordered = sorted(values)
    rendered = [str(value) for value in ordered[:limit]]
    if len(ordered) > limit:
        rendered.append(f"...(+{len(ordered) - limit})")
    return "; ".join(rendered)


def compact_channel_coverage(profile: dict[str, Any], limit: int = 12) -> str:
    channels = sorted(profile["channels"])
    rendered: list[str] = []
    for channel in channels[:limit]:
        z_count = len(profile["channel_z_values"].get(channel, set()))
        plane_count = len(profile["channel_plane_values"].get(channel, set()))
        if z_count and plane_count:
            rendered.append(f"{channel}:z{z_count}/p{plane_count}")
        elif z_count:
            rendered.append(f"{channel}:z{z_count}")
        elif plane_count:
            rendered.append(f"{channel}:p{plane_count}")
        else:
            rendered.append(f"{channel}:unquantified")
    if len(channels) > limit:
        rendered.append(f"...(+{len(channels) - limit})")
    return "; ".join(rendered)


def make_row(
    record: dict[str, Any],
    profile: dict[str, Any],
    stack_profile: dict[str, Any],
    score: int,
    raw_status: str,
    cell_lines: list[str],
    path_context: list[str],
    focus_matches: str,
    focus_score: str,
    focus_best_match: str,
    example_key: str,
    example_limit: int,
) -> dict[str, Any]:
    examples = profile[example_key] or profile["image_examples"]
    min_spacing, median_spacing, spacing_unit = spacing_summary(profile)
    return {
        "score": score,
        "path": record["path"],
        "status": record.get("status", ""),
        "raw_or_derived": raw_status,
        "Nentries": "" if record.get("Nentries") is None else record.get("Nentries"),
        "Nfiles": "" if record.get("Nfiles") is None else record.get("Nfiles"),
        "image_like_count": image_like_count(record),
        "sample_image_count": len(profile["image_examples"]),
        "stack_tiff_count_sampled": stack_profile["stack_tiff_count"],
        "stack_tiff_kind": compact_values(stack_profile["stack_tiff_kind"]),
        "stack_tiff_ome_sizes": compact_values(stack_profile["stack_tiff_ome_sizes"]),
        "stack_tiff_pages": compact_values(stack_profile["stack_tiff_pages"]),
        "stack_tiff_shapes": compact_values(stack_profile["stack_tiff_shapes"]),
        "distinct_z": len(profile["z_values"]),
        "distinct_planes": len(profile["plane_values"]),
        "distinct_timepoints": profile["timepoint_count"],
        "min_time_spacing_minutes": min_spacing,
        "median_time_spacing_minutes": median_spacing,
        "time_spacing_unit": spacing_unit,
        "distinct_wells": len(profile["wells"]),
        "distinct_fields_or_sites": len(profile["fields"]),
        "channels": compact_values(profile["channels"]),
        "distinct_channels": len(profile["channels"]),
        "has_multiple_channels": "1" if len(profile["channels"]) >= 2 else "0",
        "channel_z_or_plane_coverage": compact_channel_coverage(profile),
        "cell_line_hits": "; ".join(cell_lines),
        "tracking_focus_score": focus_score,
        "tracking_focus_best_match": focus_best_match,
        "tracking_focus_matches": focus_matches,
        "path_context": " / ".join(path_context),
        "z_values": compact_values(profile["z_values"]),
        "time_values": compact_values(
            set(profile["zeiss_t_values"]) | set(profile["elapsed_minutes"]) | set(profile["inc_datetimes"])
        ),
        "stack_tiff_examples": "; ".join(stack_profile["stack_tiff_examples"][:example_limit]),
        "stack_tiff_errors": "; ".join(stack_profile["stack_tiff_errors"]),
        "examples": "; ".join(examples[:example_limit]),
    }


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "score",
        "path",
        "status",
        "raw_or_derived",
        "Nentries",
        "Nfiles",
        "image_like_count",
        "sample_image_count",
        "stack_tiff_count_sampled",
        "stack_tiff_kind",
        "stack_tiff_ome_sizes",
        "stack_tiff_pages",
        "stack_tiff_shapes",
        "distinct_z",
        "distinct_planes",
        "distinct_timepoints",
        "min_time_spacing_minutes",
        "median_time_spacing_minutes",
        "time_spacing_unit",
        "distinct_wells",
        "distinct_fields_or_sites",
        "channels",
        "distinct_channels",
        "has_multiple_channels",
        "channel_z_or_plane_coverage",
        "cell_line_hits",
        "tracking_focus_score",
        "tracking_focus_best_match",
        "tracking_focus_matches",
        "path_context",
        "z_values",
        "time_values",
        "stack_tiff_examples",
        "stack_tiff_errors",
        "examples",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    z_rows: list[dict[str, Any]] = []
    tracking_rows: list[dict[str, Any]] = []
    focused_tracking_rows: list[dict[str, Any]] = []

    with Path(args.input_jsonl).open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            names = record.get("sampleNames") or record.get("example_names") or []
            if not names:
                continue
            if is_excluded_path(record["path"]):
                continue

            profile = profile_names(names)
            stack_profile = inspect_tiff_stacks(record, names, args.metadata_stack_scan_limit)
            count = image_like_count(record)
            raw_status = status_label(record["path"], names)
            cell_lines, path_context = collect_path_context(record["path"])
            focus_matches, focus_score, focus_best_match = path_focus_match(record["path"])

            z_score = score_z(profile, raw_status, count, cell_lines, stack_profile)
            if z_score >= args.min_score:
                z_rows.append(
                    make_row(
                        record,
                        profile,
                        stack_profile,
                        z_score,
                        raw_status,
                        cell_lines,
                        path_context,
                        focus_matches,
                        focus_score,
                        focus_best_match,
                        "z_examples",
                        args.example_limit,
                    )
                )

            tracking_score = score_tracking(profile, raw_status, count, cell_lines, stack_profile)
            if tracking_score >= args.min_score:
                tracking_rows.append(
                    make_row(
                        record,
                        profile,
                        stack_profile,
                        tracking_score,
                        raw_status,
                        cell_lines,
                        path_context,
                        focus_matches,
                        focus_score,
                        focus_best_match,
                        "t_examples",
                        args.example_limit,
                    )
                )
                if focus_matches:
                    focused_tracking_rows.append(tracking_rows[-1])

    z_rows.sort(key=lambda row: (-int(row["score"]), -int(row["distinct_z"]), -int(row["image_like_count"]), row["path"]))
    tracking_rows.sort(
        key=lambda row: (
            row["min_time_spacing_minutes"] == "",
            int(row["min_time_spacing_minutes"] or 10**9),
            -int(row["distinct_timepoints"]),
            -int(row["score"]),
            -int(row["image_like_count"]),
            row["path"],
        )
    )
    focused_tracking_rows.sort(
        key=lambda row: (
            -float(row["tracking_focus_score"] or 0),
            row["min_time_spacing_minutes"] == "",
            int(row["min_time_spacing_minutes"] or 10**9),
            -int(row["distinct_timepoints"]),
            -int(row["score"]),
            -int(row["image_like_count"]),
            row["path"],
        )
    )

    write_rows(Path(args.z_output), z_rows)
    write_rows(Path(args.tracking_output), tracking_rows)
    write_rows(Path(args.focused_tracking_output), focused_tracking_rows)
    print(f"wrote {len(z_rows)} z-stack candidate rows to {args.z_output}")
    print(f"wrote {len(tracking_rows)} tracking candidate rows to {args.tracking_output}")
    print(f"wrote {len(focused_tracking_rows)} focused tracking candidate rows to {args.focused_tracking_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
