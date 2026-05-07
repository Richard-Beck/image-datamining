#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from _bootstrap import add_src_to_path

add_src_to_path()

from image_datamining.fucci_3d_segmentation.io import read_zcyx_from_manifest_row, write_zcyx_tiff
from image_datamining.fucci_3d_segmentation.manifest import read_manifest
from image_datamining.fucci_3d_segmentation.preprocess import build_bf_fucci_fused_input, centered_xy_crop


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare cropped/full ZCYX inputs for segmentation experiments.")
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--crop-size", type=int, default=0, help="Centered XY crop size. Use 0 for full XY.")
    parser.add_argument("--limit", type=int, default=0, help="Optional number of manifest rows to prepare.")
    parser.add_argument("--invert-bf", action="store_true")
    parser.add_argument("--input-kind", choices=["raw_zcyx", "bf_fucci_fused"], default="bf_fucci_fused")
    args = parser.parse_args()

    run_dir = args.run_dir
    prepared_dir = run_dir / "prepared_inputs"
    prepared_dir.mkdir(parents=True, exist_ok=True)

    rows = read_manifest(args.manifest)
    if args.limit:
        rows = rows[: args.limit]

    config = {
        "manifest": str(args.manifest),
        "crop_size": args.crop_size,
        "invert_bf": args.invert_bf,
        "input_kind": args.input_kind,
    }
    (run_dir / "config.prepare.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")

    output_manifest = run_dir / "prepared_inputs.tsv"
    with output_manifest.open("w", newline="") as handle:
        fieldnames = ["sample_id", "cell_line", "source_dir", "prepared_input", "input_kind", "crop_size", "invert_bf"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in rows:
            zcyx = read_zcyx_from_manifest_row(row)
            zcyx = centered_xy_crop(zcyx, args.crop_size)
            if args.input_kind == "bf_fucci_fused":
                out = build_bf_fucci_fused_input(zcyx, invert_bf=args.invert_bf)
                dtype_tag = "float32"
            else:
                out = zcyx
                dtype_tag = str(out.dtype)

            sample_id = row["sample_id"]
            prepared_input = prepared_dir / f"{sample_id}_{args.input_kind}_crop{args.crop_size}_{dtype_tag}_ZCYX.tif"
            write_zcyx_tiff(prepared_input, out)
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "cell_line": row["cell_line"],
                    "source_dir": row["source_dir"],
                    "prepared_input": str(prepared_input),
                    "input_kind": args.input_kind,
                    "crop_size": args.crop_size,
                    "invert_bf": int(args.invert_bf),
                }
            )
            print(f"wrote {prepared_input}")

    print(f"wrote prepared manifest to {output_manifest}")


if __name__ == "__main__":
    main()

