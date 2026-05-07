#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from _bootstrap import add_src_to_path

add_src_to_path()

from image_datamining.fucci_3d_segmentation.manifest import (
    discover_z_stack_fields,
    read_manifest,
    write_manifest,
)


DEFAULT_ROOTS = [Path("/share/andor_lab/Jackson/FUCCI")]

DEFAULT_VARIANTS = [
    {
        "name": "sd_w9",
        "method": "local_sd",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 1.2,
    },
    {
        "name": "sd_w17",
        "method": "local_sd",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 17,
        "log_sigma": 1.2,
    },
    {
        "name": "ten_hp8",
        "method": "tenengrad",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 1.2,
    },
    {
        "name": "ten_hp16",
        "method": "tenengrad",
        "background_sigma": 40,
        "highpass_sigma": 16,
        "local_std_window": 9,
        "log_sigma": 1.2,
    },
    {
        "name": "log1.2",
        "method": "log",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 1.2,
    },
    {
        "name": "log2.0",
        "method": "log",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 2.0,
    },
    {
        "name": "combo_base",
        "method": "composite",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 1.2,
    },
    {
        "name": "combo_w17",
        "method": "composite",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 17,
        "log_sigma": 1.2,
    },
]

LOG_ONLY_VARIANTS = [
    {
        "name": "log0.7",
        "method": "log",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 0.7,
    },
    {
        "name": "log1.0",
        "method": "log",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 1.0,
    },
    {
        "name": "log1.5",
        "method": "log",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 1.5,
    },
    {
        "name": "log2.0",
        "method": "log",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 2.0,
    },
    {
        "name": "log3.0",
        "method": "log",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 3.0,
    },
    {
        "name": "log4.0",
        "method": "log",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 4.0,
    },
    {
        "name": "log6.0",
        "method": "log",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 6.0,
    },
    {
        "name": "log8.0",
        "method": "log",
        "background_sigma": 24,
        "highpass_sigma": 8,
        "local_std_window": 9,
        "log_sigma": 8.0,
    },
]


def parse_float_list(text: str) -> list[float]:
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("Expected at least one comma-separated number")
    return values


def parse_int_list(text: str) -> list[int]:
    values = [int(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("Expected at least one comma-separated integer")
    return values


def read_channel_stack(paths_text: str):
    import tifffile as tf

    from image_datamining.fucci_3d_segmentation.io import read_channel_file_list

    return tf.imread([str(path) for path in read_channel_file_list(paths_text)])


def append_summary(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "cell_line",
        "source_dir",
        "bf_channel",
        "crop_size",
        "central_z",
        "central_grid_png",
        "projection_grid_png",
        "variant_names",
        "raw_bf_crop_tif",
        "focus_confidence_tif",
        "focus_score_tif",
        "extra_outputs_json",
    ]
    exists = path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create small focus-channel previews from FUCCI brightfield z-stacks.")
    parser.add_argument("--root", action="append", type=Path, help="Root to discover recursively. May be repeated.")
    parser.add_argument("--manifest", type=Path, help="Existing z-stack manifest. If omitted, roots are discovered.")
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--crop-size", type=int, default=256)
    parser.add_argument("--limit", type=int, default=6)
    parser.add_argument("--bf-channel", choices=["ch00", "ch01", "ch02"], default="ch01")
    parser.add_argument("--background-sigma", type=float, default=24)
    parser.add_argument("--highpass-sigma", type=float, default=8)
    parser.add_argument("--local-std-window", type=int, default=9)
    parser.add_argument("--log-sigma", type=float, default=1.2)
    parser.add_argument("--method", choices=["local_sd", "tenengrad", "log", "composite"], default="composite")
    parser.add_argument("--sd-window", type=int, default=12, help="Local SD window for --preset sd-blur.")
    parser.add_argument("--combo-sd-window", type=int, default=10, help="Local SD window for --preset sd-log-combo.")
    parser.add_argument("--combo-log-sigma", type=float, default=2.0, help="LoG sigma for --preset sd-log-combo.")
    parser.add_argument("--blur-sigma-xy", type=float, default=20, help="XY Gaussian sigma for --preset sd-blur.")
    parser.add_argument(
        "--grid-method",
        choices=["local_sd", "tenengrad", "log", "composite"],
        default="log",
        help="Detector method for --preset detector-blur-grid.",
    )
    parser.add_argument(
        "--detector-scales",
        type=parse_float_list,
        default=parse_float_list("1.0,2.0,4.0"),
        help="Comma-separated detector scales for --preset detector-blur-grid. For log this is LoG sigma; for local_sd this is window size.",
    )
    parser.add_argument(
        "--blur-sigmas-xy",
        type=parse_float_list,
        default=parse_float_list("5,10,20,30,40"),
        help="Comma-separated XY blur sigmas for --preset detector-blur-grid.",
    )
    parser.add_argument(
        "--z-anisotropy",
        type=float,
        default=1.44,
        help="Z spacing divided by XY spacing; blur presets use sigma_z = sigma_xy / z_anisotropy.",
    )
    parser.add_argument(
        "--preset",
        choices=["mixed", "log-scales", "sd-blur", "detector-blur-grid", "sd-log-combo", "abs-zscore-smooth"],
        default="mixed",
        help="Default grid preset. 'detector-blur-grid' makes detector-scale x blur-scale product grids.",
    )
    parser.add_argument("--save-stacks", action="store_true", help="Also write raw BF, focus confidence, and score TIFF stacks.")
    parser.add_argument(
        "--single-setting",
        action="store_true",
        help="Only compute the explicitly supplied focus setting instead of the default 8-setting grid.",
    )
    args = parser.parse_args()

    run_dir = args.run_dir
    png_dir = run_dir / "png"
    stack_dir = run_dir / "stacks"
    png_dir.mkdir(parents=True, exist_ok=True)
    stack_dir.mkdir(parents=True, exist_ok=True)

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
        "crop_size": args.crop_size,
        "limit": args.limit,
        "bf_channel": args.bf_channel,
        "background_sigma": args.background_sigma,
        "highpass_sigma": args.highpass_sigma,
        "local_std_window": args.local_std_window,
        "log_sigma": args.log_sigma,
        "method": args.method,
        "sd_window": args.sd_window,
        "combo_sd_window": args.combo_sd_window,
        "combo_log_sigma": args.combo_log_sigma,
        "blur_sigma_xy": args.blur_sigma_xy,
        "grid_method": args.grid_method,
        "detector_scales": args.detector_scales,
        "blur_sigmas_xy": args.blur_sigmas_xy,
        "z_anisotropy": args.z_anisotropy,
        "preset": args.preset,
        "save_stacks": args.save_stacks,
        "single_setting": args.single_setting,
        "default_variants": DEFAULT_VARIANTS,
        "log_only_variants": LOG_ONLY_VARIANTS,
    }
    (run_dir / "config.focus_preview.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")

    summary_path = run_dir / "focus_preview_outputs.tsv"
    for i, row in enumerate(rows, start=1):
        import tifffile as tf

        from image_datamining.fucci_3d_segmentation.preprocess import centered_xy_crop
        from image_datamining.fucci_focus_channel.focus import (
            build_abs_zscore_smooth,
            build_component_focus_confidence,
            build_sd_blur_product,
            build_sd_log_combo,
            otsu_threshold,
        )
        from image_datamining.fucci_focus_channel.viz import (
            save_labeled_grid_png,
        )

        sample_id = row["sample_id"]
        print(f"[{i}/{len(rows)}] {sample_id}", flush=True)
        bf = read_channel_stack(row[f"{args.bf_channel}_files"])
        bf_crop = centered_xy_crop(bf, args.crop_size)

        central_z = int(bf_crop.shape[0] // 2)
        raw_tif_text = ""
        confidence_tif_text = ""
        score_tif_text = ""
        extra_outputs: dict[str, str | float] = {}

        if args.preset == "abs-zscore-smooth" and not args.single_setting:
            result = build_abs_zscore_smooth(
                bf_crop,
                blur_sigma_xy=args.blur_sigma_xy,
                z_anisotropy=args.z_anisotropy,
            )
            threshold = otsu_threshold(result["smooth"])
            mask = (result["smooth"] >= threshold).astype("uint8")
            labels = [
                f"{args.bf_channel} raw",
                "abs global z",
                f"smooth b{args.blur_sigma_xy:g}",
                "otsu(smooth)",
            ]
            central_images = [
                bf_crop[central_z],
                result["abs_zscore"][central_z],
                result["smooth"][central_z],
                mask[central_z],
            ]
            projection_images = [
                bf_crop.sum(axis=0),
                result["abs_zscore"].sum(axis=0),
                result["smooth"].sum(axis=0),
                mask.sum(axis=0),
            ]
            central_grid_png = png_dir / f"{sample_id}_z{central_z:02d}_abs_zscore_smooth_grid.png"
            projection_grid_png = png_dir / f"{sample_id}_zsum_abs_zscore_smooth_grid.png"
            save_labeled_grid_png(central_grid_png, central_images, labels, ncols=2)
            save_labeled_grid_png(projection_grid_png, projection_images, labels, ncols=2)

            raw_tif = stack_dir / f"{sample_id}_raw_bf_{args.bf_channel}_crop{args.crop_size}_ZYX.tif"
            abs_tif = stack_dir / f"{sample_id}_abs_global_zscore_ZYX.tif"
            smooth_tif = stack_dir / f"{sample_id}_abs_global_zscore_smooth_b{args.blur_sigma_xy:g}_ZYX.tif"
            mask_tif = stack_dir / f"{sample_id}_abs_global_zscore_smooth_otsu_mask_ZYX.tif"
            tf.imwrite(raw_tif, bf_crop)
            tf.imwrite(abs_tif, result["abs_zscore"].astype("float32"))
            tf.imwrite(smooth_tif, result["smooth"].astype("float32"))
            tf.imwrite(mask_tif, mask.astype("uint8"))
            raw_tif_text = str(raw_tif)
            confidence_tif_text = str(smooth_tif)
            score_tif_text = str(abs_tif)
            extra_outputs = {
                "abs_zscore_tif": str(abs_tif),
                "smooth_tif": str(smooth_tif),
                "otsu_mask_tif": str(mask_tif),
                "otsu_threshold": float(threshold),
            }

            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "cell_line": row["cell_line"],
                    "source_dir": row["source_dir"],
                    "bf_channel": args.bf_channel,
                    "crop_size": args.crop_size,
                    "central_z": central_z,
                    "central_grid_png": str(central_grid_png),
                    "projection_grid_png": str(projection_grid_png),
                    "variant_names": ";".join(labels),
                    "raw_bf_crop_tif": raw_tif_text,
                    "focus_confidence_tif": confidence_tif_text,
                    "focus_score_tif": score_tif_text,
                    "extra_outputs_json": json.dumps(extra_outputs, sort_keys=True),
                },
            )
            continue

        if args.preset == "sd-log-combo" and not args.single_setting:
            combo = build_sd_log_combo(
                bf_crop,
                sd_window=args.combo_sd_window,
                log_sigma=args.combo_log_sigma,
                blur_sigma_xy=args.blur_sigma_xy,
                z_anisotropy=args.z_anisotropy,
                background_sigma=args.background_sigma,
                highpass_sigma=args.highpass_sigma,
            )
            labels = [
                f"{args.bf_channel} raw",
                f"sd_w{args.combo_sd_window}",
                f"log{args.combo_log_sigma:g}",
                "sqrt(sd*log)",
                f"blur b{args.blur_sigma_xy:g}",
                "combo*blur",
                "otsu(blur)",
                "product*otsu",
            ]
            central_images = [
                bf_crop[central_z],
                combo["sd_confidence"][central_z],
                combo["log_confidence"][central_z],
                combo["combined"][central_z],
                combo["blurred"][central_z],
                combo["product"][central_z],
                combo["otsu_mask"][central_z],
                combo["otsu_product"][central_z],
            ]
            projection_images = [
                bf_crop.sum(axis=0),
                combo["sd_confidence"].sum(axis=0),
                combo["log_confidence"].sum(axis=0),
                combo["combined"].sum(axis=0),
                combo["blurred"].sum(axis=0),
                combo["product"].sum(axis=0),
                combo["otsu_mask"].sum(axis=0),
                combo["otsu_product"].sum(axis=0),
            ]
            central_grid_png = png_dir / f"{sample_id}_z{central_z:02d}_sd_log_combo_grid.png"
            projection_grid_png = png_dir / f"{sample_id}_zsum_sd_log_combo_grid.png"
            save_labeled_grid_png(central_grid_png, central_images, labels, ncols=4)
            save_labeled_grid_png(projection_grid_png, projection_images, labels, ncols=4)

            raw_tif = stack_dir / f"{sample_id}_raw_bf_{args.bf_channel}_crop{args.crop_size}_ZYX.tif"
            sd_tif = stack_dir / f"{sample_id}_sd_w{args.combo_sd_window}_confidence_ZYX.tif"
            log_tif = stack_dir / f"{sample_id}_log{args.combo_log_sigma:g}_confidence_ZYX.tif"
            combined_tif = stack_dir / f"{sample_id}_sd_log_combined_ZYX.tif"
            blurred_tif = stack_dir / f"{sample_id}_sd_log_blur_b{args.blur_sigma_xy:g}_ZYX.tif"
            product_tif = stack_dir / f"{sample_id}_sd_log_product_ZYX.tif"
            otsu_mask_tif = stack_dir / f"{sample_id}_sd_log_otsu_blur_mask_ZYX.tif"
            otsu_product_tif = stack_dir / f"{sample_id}_sd_log_otsu_product_ZYX.tif"
            tf.imwrite(raw_tif, bf_crop)
            tf.imwrite(sd_tif, combo["sd_confidence"].astype("float32"))
            tf.imwrite(log_tif, combo["log_confidence"].astype("float32"))
            tf.imwrite(combined_tif, combo["combined"].astype("float32"))
            tf.imwrite(blurred_tif, combo["blurred"].astype("float32"))
            tf.imwrite(product_tif, combo["product"].astype("float32"))
            tf.imwrite(otsu_mask_tif, combo["otsu_mask"].astype("uint8"))
            tf.imwrite(otsu_product_tif, combo["otsu_product"].astype("float32"))
            raw_tif_text = str(raw_tif)
            confidence_tif_text = str(product_tif)
            score_tif_text = str(combined_tif)
            extra_outputs = {
                "sd_confidence_tif": str(sd_tif),
                "log_confidence_tif": str(log_tif),
                "combined_tif": str(combined_tif),
                "blurred_tif": str(blurred_tif),
                "product_tif": str(product_tif),
                "otsu_mask_tif": str(otsu_mask_tif),
                "otsu_product_tif": str(otsu_product_tif),
                "otsu_threshold": float(combo["otsu_threshold"]),
            }

            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "cell_line": row["cell_line"],
                    "source_dir": row["source_dir"],
                    "bf_channel": args.bf_channel,
                    "crop_size": args.crop_size,
                    "central_z": central_z,
                    "central_grid_png": str(central_grid_png),
                    "projection_grid_png": str(projection_grid_png),
                    "variant_names": ";".join(labels),
                    "raw_bf_crop_tif": raw_tif_text,
                    "focus_confidence_tif": confidence_tif_text,
                    "focus_score_tif": score_tif_text,
                    "extra_outputs_json": json.dumps(extra_outputs, sort_keys=True),
                },
            )
            continue

        if args.preset == "sd-blur" and not args.single_setting:
            blur, sd_confidence, product, score = build_sd_blur_product(
                bf_crop,
                local_std_window=args.sd_window,
                blur_sigma_xy=args.blur_sigma_xy,
                z_anisotropy=args.z_anisotropy,
                background_sigma=args.background_sigma,
                highpass_sigma=args.highpass_sigma,
            )
            labels = [
                f"{args.bf_channel} raw",
                f"blur sd sigma{args.blur_sigma_xy:g}",
                f"sd_w{args.sd_window}",
                "sd * blur",
            ]
            central_images = [bf_crop[central_z], blur[central_z], sd_confidence[central_z], product[central_z]]
            projection_images = [bf_crop.sum(axis=0), blur.sum(axis=0), sd_confidence.sum(axis=0), product.sum(axis=0)]
            central_grid_png = png_dir / f"{sample_id}_z{central_z:02d}_sd_blur_grid_2x2.png"
            projection_grid_png = png_dir / f"{sample_id}_zsum_sd_blur_grid_2x2.png"
            save_labeled_grid_png(central_grid_png, central_images, labels, ncols=2)
            save_labeled_grid_png(projection_grid_png, projection_images, labels, ncols=2)

            if args.save_stacks:
                raw_tif = stack_dir / f"{sample_id}_raw_bf_{args.bf_channel}_crop{args.crop_size}_ZYX.tif"
                confidence_tif = stack_dir / f"{sample_id}_sd_blur_product_crop{args.crop_size}_ZYX.tif"
                score_tif = stack_dir / f"{sample_id}_sd_score_crop{args.crop_size}_ZYX.tif"
                tf.imwrite(raw_tif, bf_crop)
                tf.imwrite(confidence_tif, product.astype("float32"))
                tf.imwrite(score_tif, score.astype("float32"))
                raw_tif_text = str(raw_tif)
                confidence_tif_text = str(confidence_tif)
                score_tif_text = str(score_tif)

            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "cell_line": row["cell_line"],
                    "source_dir": row["source_dir"],
                    "bf_channel": args.bf_channel,
                    "crop_size": args.crop_size,
                    "central_z": central_z,
                    "central_grid_png": str(central_grid_png),
                    "projection_grid_png": str(projection_grid_png),
                    "variant_names": ";".join(labels),
                    "raw_bf_crop_tif": raw_tif_text,
                    "focus_confidence_tif": confidence_tif_text,
                    "focus_score_tif": score_tif_text,
                    "extra_outputs_json": json.dumps(extra_outputs, sort_keys=True),
                },
            )
            continue

        if args.preset == "detector-blur-grid" and not args.single_setting:
            import numpy as np
            from scipy import ndimage as ndi

            central_images = []
            projection_images = []
            labels = []
            base_confidence = None
            base_score = None

            for detector_scale in args.detector_scales:
                for blur_sigma_xy in args.blur_sigmas_xy:
                    local_std_window = args.local_std_window
                    log_sigma = args.log_sigma
                    label_scale = detector_scale
                    if args.grid_method == "local_sd":
                        local_std_window = int(round(detector_scale))
                    elif args.grid_method == "log":
                        log_sigma = detector_scale

                    confidence, score, _bf_norm, _bf_hp = build_component_focus_confidence(
                        bf_crop,
                        method=args.grid_method,
                        confidence_mode="relative",
                        background_sigma=args.background_sigma,
                        highpass_sigma=args.highpass_sigma,
                        local_std_window=local_std_window,
                        log_sigma=log_sigma,
                    )
                    if base_confidence is None:
                        base_confidence = confidence
                        base_score = score

                    sigma_z = blur_sigma_xy / args.z_anisotropy
                    blurred_support = ndi.gaussian_filter(
                        confidence,
                        sigma=(sigma_z, blur_sigma_xy, blur_sigma_xy),
                    ).astype(np.float32)
                    product = (confidence * blurred_support).astype(np.float32)

                    central_images.append(product[central_z])
                    projection_images.append(product.sum(axis=0))
                    if args.grid_method == "local_sd":
                        labels.append(f"sd_w{int(round(label_scale))} b{blur_sigma_xy:g}")
                    elif args.grid_method == "log":
                        labels.append(f"log{label_scale:g} b{blur_sigma_xy:g}")
                    else:
                        labels.append(f"{args.grid_method} {label_scale:g} b{blur_sigma_xy:g}")

            if base_confidence is None or base_score is None:
                raise RuntimeError("No detector-blur variants were computed")

            ncols = len(args.blur_sigmas_xy)
            central_grid_png = png_dir / f"{sample_id}_z{central_z:02d}_{args.grid_method}_blur_grid.png"
            projection_grid_png = png_dir / f"{sample_id}_zsum_{args.grid_method}_blur_grid.png"
            save_labeled_grid_png(central_grid_png, central_images, labels, ncols=ncols)
            save_labeled_grid_png(projection_grid_png, projection_images, labels, ncols=ncols)

            if args.save_stacks:
                raw_tif = stack_dir / f"{sample_id}_raw_bf_{args.bf_channel}_crop{args.crop_size}_ZYX.tif"
                confidence_tif = stack_dir / f"{sample_id}_{args.grid_method}_confidence_crop{args.crop_size}_ZYX.tif"
                score_tif = stack_dir / f"{sample_id}_{args.grid_method}_score_crop{args.crop_size}_ZYX.tif"
                tf.imwrite(raw_tif, bf_crop)
                tf.imwrite(confidence_tif, base_confidence.astype("float32"))
                tf.imwrite(score_tif, base_score.astype("float32"))
                raw_tif_text = str(raw_tif)
                confidence_tif_text = str(confidence_tif)
                score_tif_text = str(score_tif)

            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "cell_line": row["cell_line"],
                    "source_dir": row["source_dir"],
                    "bf_channel": args.bf_channel,
                    "crop_size": args.crop_size,
                    "central_z": central_z,
                    "central_grid_png": str(central_grid_png),
                    "projection_grid_png": str(projection_grid_png),
                    "variant_names": ";".join(labels),
                    "raw_bf_crop_tif": raw_tif_text,
                    "focus_confidence_tif": confidence_tif_text,
                    "focus_score_tif": score_tif_text,
                    "extra_outputs_json": json.dumps(extra_outputs, sort_keys=True),
                },
            )
            continue

        variants = LOG_ONLY_VARIANTS if args.preset == "log-scales" else DEFAULT_VARIANTS
        if args.single_setting:
            variants = [
                {
                    "name": "single",
                    "method": args.method,
                    "background_sigma": args.background_sigma,
                    "highpass_sigma": args.highpass_sigma,
                    "local_std_window": args.local_std_window,
                    "log_sigma": args.log_sigma,
                }
            ]

        central_images = [bf_crop[central_z]]
        projection_images = [bf_crop.sum(axis=0)]
        labels = [f"{args.bf_channel} raw"]
        base_confidence = None
        base_score = None

        for variant in variants:
            confidence, score, _bf_norm, _bf_hp = build_component_focus_confidence(
                bf_crop,
                method=str(variant["method"]),
                confidence_mode="relative",
                background_sigma=float(variant["background_sigma"]),
                highpass_sigma=float(variant["highpass_sigma"]),
                local_std_window=int(variant["local_std_window"]),
                log_sigma=float(variant["log_sigma"]),
            )
            if base_confidence is None:
                base_confidence = confidence
                base_score = score
            central_images.append(confidence[central_z])
            projection_images.append(confidence.sum(axis=0))
            labels.append(str(variant["name"]))

        if base_confidence is None or base_score is None:
            raise RuntimeError("No focus variants were computed")

        central_grid_png = png_dir / f"{sample_id}_z{central_z:02d}_focus_grid_3x3.png"
        projection_grid_png = png_dir / f"{sample_id}_zsum_focus_grid_3x3.png"

        save_labeled_grid_png(central_grid_png, central_images, labels, ncols=3)
        save_labeled_grid_png(projection_grid_png, projection_images, labels, ncols=3)
        if args.save_stacks:
            raw_tif = stack_dir / f"{sample_id}_raw_bf_{args.bf_channel}_crop{args.crop_size}_ZYX.tif"
            confidence_tif = stack_dir / f"{sample_id}_focus_confidence_crop{args.crop_size}_ZYX.tif"
            score_tif = stack_dir / f"{sample_id}_focus_score_crop{args.crop_size}_ZYX.tif"
            tf.imwrite(raw_tif, bf_crop)
            tf.imwrite(confidence_tif, base_confidence.astype("float32"))
            tf.imwrite(score_tif, base_score.astype("float32"))
            raw_tif_text = str(raw_tif)
            confidence_tif_text = str(confidence_tif)
            score_tif_text = str(score_tif)

        append_summary(
            summary_path,
            {
                "sample_id": sample_id,
                "cell_line": row["cell_line"],
                "source_dir": row["source_dir"],
                "bf_channel": args.bf_channel,
                "crop_size": args.crop_size,
                "central_z": central_z,
                "central_grid_png": str(central_grid_png),
                "projection_grid_png": str(projection_grid_png),
                "variant_names": ";".join(labels),
                "raw_bf_crop_tif": raw_tif_text,
                "focus_confidence_tif": confidence_tif_text,
                "focus_score_tif": score_tif_text,
                "extra_outputs_json": json.dumps(extra_outputs, sort_keys=True),
            },
        )

    print(f"wrote summary to {summary_path}")


if __name__ == "__main__":
    main()
