#!/usr/bin/env python3
"""Rank directory summaries by count of image-like filenames."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank name-based directory summaries by image-like entry counts."
    )
    parser.add_argument("input_jsonl", help="Input JSONL from walk_directory_tree.py or inventory_dirs.py.")
    parser.add_argument("--output", required=True, help="Output TSV path.")
    parser.add_argument(
        "--min-matches",
        type=int,
        default=1,
        help="Only include directories with at least this many image-like entries.",
    )
    parser.add_argument(
        "--example-limit",
        type=int,
        default=10,
        help="Maximum matching example filenames to include per row.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = []

    with Path(args.input_jsonl).open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            ext_counts = record.get("extension_counts") or record.get("name_extension_counts", {})
            matches = sum(ext_counts.get(ext, 0) for ext in IMAGE_EXTENSIONS)
            if matches < args.min_matches:
                continue

            sample_names = record.get("sampleNames") or record.get("example_names", [])
            examples = [
                name
                for name in sample_names
                if Path(name).suffix.lower() in IMAGE_EXTENSIONS
            ][: args.example_limit]

            rows.append(
                {
                    "path": record["path"],
                    "Nmatches": matches,
                    "exampleFilenames": "; ".join(examples),
                }
            )

    rows.sort(key=lambda row: (-row["Nmatches"], row["path"]))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["path", "Nmatches", "exampleFilenames"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows)} rows to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
