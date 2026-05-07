#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from _bootstrap import add_src_to_path

add_src_to_path()

from image_datamining.fucci_3d_segmentation.manifest import discover_z_stack_fields, write_manifest


DEFAULT_ROOTS = [
    Path("/share/andor_lab/Jackson/FUCCI/NCI-N87/fucci_ncin87"),
]


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a manifest of FUCCI split-channel z-stack fields.")
    parser.add_argument("--root", action="append", type=Path, help="Root containing field directories. May be repeated.")
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Search below each root for field directories instead of only checking direct children.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("pipelines/fucci_3d_segmentation/manifests/jackson_fucci_zstacks.tsv"),
    )
    args = parser.parse_args()

    roots = args.root if args.root else DEFAULT_ROOTS
    records = discover_z_stack_fields(roots, recursive=args.recursive)
    write_manifest(records, args.output)
    print(f"wrote {len(records)} fields to {args.output}")


if __name__ == "__main__":
    main()
