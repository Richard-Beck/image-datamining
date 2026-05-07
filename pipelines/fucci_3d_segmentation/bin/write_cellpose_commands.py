#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from _bootstrap import add_src_to_path

add_src_to_path()

from image_datamining.fucci_3d_segmentation.commands import write_cellpose_command_manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Write a CellposeSAM command manifest without running segmentation.")
    parser.add_argument("--prepared-manifest", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--mode", choices=["3d", "2d-per-z"], required=True)
    parser.add_argument("--anisotropy", type=float, default=1.44)
    parser.add_argument("--min-size", type=int, default=2000)
    parser.add_argument("--cpu", action="store_true", help="Generate commands without --gpu.")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    output = args.output or args.run_dir / f"cellpose_commands_{args.mode}.tsv"
    write_cellpose_command_manifest(
        prepared_manifest=args.prepared_manifest,
        run_dir=args.run_dir,
        output=output,
        mode=args.mode,
        anisotropy=args.anisotropy,
        min_size=args.min_size,
        gpu=not args.cpu,
    )
    print(f"wrote command manifest to {output}")


if __name__ == "__main__":
    main()

