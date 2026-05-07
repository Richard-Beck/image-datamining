#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
import traceback
from pathlib import Path

from _bootstrap import add_src_to_path

add_src_to_path()

from image_datamining.fucci_3d_segmentation.manifest import (
    discover_z_stack_fields,
    read_manifest,
    write_manifest,
)


DEFAULT_ROOTS = [Path("/share/andor_lab/Jackson/FUCCI")]


def prepare_row(
    row: dict[str, str],
    input_kind: str,
    crop_size: int,
    invert_bf: bool,
    focus_blur_sigma_xy: float,
    focus_z_anisotropy: float,
):
    from image_datamining.fucci_3d_segmentation.io import read_zcyx_from_manifest_row
    from image_datamining.fucci_3d_segmentation.preprocess import (
        build_bf_fucci_fused_input,
        build_fucci_bf_absz_product_input,
        centered_xy_crop,
    )

    zcyx = read_zcyx_from_manifest_row(row)
    raw_crop = centered_xy_crop(zcyx, crop_size)
    if input_kind == "bf_fucci_fused":
        cellpose_input = build_bf_fucci_fused_input(raw_crop, invert_bf=invert_bf)
        return raw_crop, cellpose_input, "float32"
    if input_kind == "fucci_bf_absz_product":
        cellpose_input = build_fucci_bf_absz_product_input(
            raw_crop,
            blur_sigma_xy=focus_blur_sigma_xy,
            z_anisotropy=focus_z_anisotropy,
        )
        return raw_crop, cellpose_input, "float32"
    return raw_crop, raw_crop, str(raw_crop.dtype)


def append_result(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "cell_line",
        "source_dir",
        "status",
        "mode",
        "mask_path",
        "raw_crop",
        "cellpose_input",
        "elapsed_seconds",
        "error",
    ]
    exists = path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stream FUCCI z-stack preparation and CellposeSAM segmentation one field at a time."
    )
    parser.add_argument("--root", action="append", type=Path, help="Root to discover recursively. May be repeated.")
    parser.add_argument("--manifest", type=Path, help="Existing input manifest. If omitted, roots are discovered.")
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--mode", choices=["3d", "2d-per-z"], default="3d")
    parser.add_argument("--crop-size", type=int, default=256)
    parser.add_argument(
        "--input-kind",
        choices=["raw_zcyx", "bf_fucci_fused", "fucci_bf_absz_product"],
        default="bf_fucci_fused",
    )
    parser.add_argument("--invert-bf", action="store_true")
    parser.add_argument(
        "--focus-blur-sigma-xy",
        type=float,
        default=20,
        help="XY Gaussian sigma for fucci_bf_absz_product focus channel.",
    )
    parser.add_argument(
        "--focus-z-anisotropy",
        type=float,
        default=1.44,
        help="Z spacing divided by XY spacing for fucci_bf_absz_product focus channel.",
    )
    parser.add_argument("--anisotropy", type=float, default=1.44)
    parser.add_argument("--min-size", type=int, default=2000)
    parser.add_argument("--cpu", action="store_true", help="Run CellposeSAM without GPU.")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--resume", action="store_true", help="Skip rows whose mask output already exists.")
    parser.add_argument(
        "--discard-prepared-inputs",
        action="store_true",
        help="Do not keep raw crop and Cellpose input TIFFs after segmentation.",
    )
    parser.add_argument("--plan-only", action="store_true", help="Write manifest/config and planned outputs, then exit.")
    args = parser.parse_args()

    run_dir = args.run_dir
    masks_dir = run_dir / "masks"
    raw_crop_dir = run_dir / "raw_crops"
    cellpose_input_dir = run_dir / "cellpose_inputs"
    logs_dir = run_dir / "logs"
    masks_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    keep_prepared_inputs = not args.discard_prepared_inputs
    if keep_prepared_inputs:
        raw_crop_dir.mkdir(parents=True, exist_ok=True)
        cellpose_input_dir.mkdir(parents=True, exist_ok=True)

    if args.manifest:
        manifest_path = args.manifest
        rows = read_manifest(manifest_path)
    else:
        roots = args.root if args.root else DEFAULT_ROOTS
        manifest_path = run_dir / "input_manifest.tsv"
        records = discover_z_stack_fields(roots, recursive=True)
        write_manifest(records, manifest_path)
        rows = read_manifest(manifest_path)

    if args.limit:
        rows = rows[: args.limit]

    config = {
        "manifest": str(manifest_path),
        "roots": [str(root) for root in (args.root or DEFAULT_ROOTS)],
        "mode": args.mode,
        "crop_size": args.crop_size,
        "input_kind": args.input_kind,
        "invert_bf": args.invert_bf,
        "focus_blur_sigma_xy": args.focus_blur_sigma_xy,
        "focus_z_anisotropy": args.focus_z_anisotropy,
        "anisotropy": args.anisotropy,
        "min_size": args.min_size,
        "gpu": not args.cpu,
        "keep_prepared_inputs": keep_prepared_inputs,
    }
    (run_dir / "config.batch_segment.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")

    plan_path = run_dir / "batch_plan.tsv"
    with plan_path.open("w", newline="") as handle:
        fieldnames = [
            "sample_id",
            "cell_line",
            "source_dir",
            "mode",
            "raw_crop",
            "cellpose_input",
            "mask_path",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            sample_id = row["sample_id"]
            raw_crop = ""
            cellpose_input = ""
            if keep_prepared_inputs:
                raw_crop = str(raw_crop_dir / f"{sample_id}_raw3ch_crop{args.crop_size}_ZCYX.tif")
                cellpose_input = str(
                    cellpose_input_dir / f"{sample_id}_{args.input_kind}_crop{args.crop_size}_ZCYX.tif"
                )
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "cell_line": row["cell_line"],
                    "source_dir": row["source_dir"],
                    "mode": args.mode,
                    "raw_crop": raw_crop,
                    "cellpose_input": cellpose_input,
                    "mask_path": str(masks_dir / f"{sample_id}_{args.mode}_cpsam_masks.tif"),
                }
            )

    print(f"planned {len(rows)} fields in {plan_path}")
    if args.plan_only:
        return

    results_path = run_dir / "batch_results.tsv"
    for i, row in enumerate(rows, start=1):
        from image_datamining.fucci_3d_segmentation.io import write_mask_tiff, write_zcyx_tiff
        from image_datamining.fucci_3d_segmentation.segment import run_cellpose_cpsam

        sample_id = row["sample_id"]
        mask_path = masks_dir / f"{sample_id}_{args.mode}_cpsam_masks.tif"
        raw_crop_path = ""
        cellpose_input_path = ""
        if keep_prepared_inputs:
            raw_crop_path = str(raw_crop_dir / f"{sample_id}_raw3ch_crop{args.crop_size}_ZCYX.tif")
            cellpose_input_path = str(
                cellpose_input_dir / f"{sample_id}_{args.input_kind}_crop{args.crop_size}_ZCYX.tif"
            )

        if args.resume and mask_path.exists():
            append_result(
                results_path,
                {
                    "sample_id": sample_id,
                    "cell_line": row["cell_line"],
                    "source_dir": row["source_dir"],
                    "status": "skipped_existing",
                    "mode": args.mode,
                    "mask_path": str(mask_path),
                    "raw_crop": raw_crop_path,
                    "cellpose_input": cellpose_input_path,
                    "elapsed_seconds": 0,
                    "error": "",
                },
            )
            print(f"[{i}/{len(rows)}] skipped existing {mask_path}")
            continue

        start = time.time()
        try:
            print(f"[{i}/{len(rows)}] preparing {sample_id}", flush=True)
            raw_crop, cellpose_input, _dtype_tag = prepare_row(
                row,
                args.input_kind,
                args.crop_size,
                args.invert_bf,
                args.focus_blur_sigma_xy,
                args.focus_z_anisotropy,
            )
            if keep_prepared_inputs:
                write_zcyx_tiff(Path(raw_crop_path), raw_crop)
                write_zcyx_tiff(Path(cellpose_input_path), cellpose_input)

            print(f"[{i}/{len(rows)}] segmenting {sample_id}", flush=True)
            masks = run_cellpose_cpsam(
                cellpose_input,
                mode=args.mode,
                anisotropy=args.anisotropy,
                min_size=args.min_size,
                gpu=not args.cpu,
            )
            write_mask_tiff(mask_path, masks)
            elapsed = round(time.time() - start, 3)
            append_result(
                results_path,
                {
                    "sample_id": sample_id,
                    "cell_line": row["cell_line"],
                    "source_dir": row["source_dir"],
                    "status": "ok",
                    "mode": args.mode,
                    "mask_path": str(mask_path),
                    "raw_crop": raw_crop_path,
                    "cellpose_input": cellpose_input_path,
                    "elapsed_seconds": elapsed,
                    "error": "",
                },
            )
            print(f"[{i}/{len(rows)}] wrote {mask_path} in {elapsed}s", flush=True)
        except Exception as exc:
            error_path = logs_dir / f"{sample_id}_{args.mode}_error.txt"
            error_path.write_text(traceback.format_exc())
            elapsed = round(time.time() - start, 3)
            append_result(
                results_path,
                {
                    "sample_id": sample_id,
                    "cell_line": row["cell_line"],
                    "source_dir": row["source_dir"],
                    "status": "error",
                    "mode": args.mode,
                    "mask_path": str(mask_path),
                    "raw_crop": raw_crop_path,
                    "cellpose_input": cellpose_input_path,
                    "elapsed_seconds": elapsed,
                    "error": f"{type(exc).__name__}: {exc}",
                },
            )
            print(f"[{i}/{len(rows)}] ERROR {sample_id}: {exc}", flush=True)

    print(f"wrote results to {results_path}")


if __name__ == "__main__":
    main()
