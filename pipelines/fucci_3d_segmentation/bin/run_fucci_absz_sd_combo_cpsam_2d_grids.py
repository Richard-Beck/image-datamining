#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import time
from pathlib import Path

import numpy as np

from _bootstrap import add_src_to_path

add_src_to_path()

import tifffile as tf
from scipy import ndimage as ndi

from image_datamining.fucci_3d_segmentation.preprocess import percentile_normalize
from image_datamining.fucci_focus_channel.focus import (
    build_abs_zscore_smooth,
    focus_component_stack,
    relative_focus_confidence,
    robust01,
)


DEFAULT_RAW_CROP_DIR = Path(
    "pipelines/fucci_3d_segmentation/runs/jackson_fucci_crop256_bf_fucci_fused_3d/raw_crops"
)
DEFAULT_RUN_DIR = Path("pipelines/fucci_3d_segmentation/runs/fucci_absz_sd_combo_cpsam_2d_grids")
RAW_CROP_SUFFIX = "_raw3ch_crop256_ZCYX.tif"


def parse_int_list(text: str) -> list[int]:
    values = [int(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("Expected at least one comma-separated integer")
    return values


def sample_id_from_raw_crop(path: Path) -> str:
    if path.name.endswith(RAW_CROP_SUFFIX):
        return path.name[: -len(RAW_CROP_SUFFIX)]
    return re.sub(r"\.tiff?$", "", path.name, flags=re.IGNORECASE)


def discover_raw_crops(raw_crop_dir: Path, pattern: str) -> list[Path]:
    return sorted(raw_crop_dir.glob(pattern))


def norm(x: np.ndarray) -> np.ndarray:
    return percentile_normalize(x).astype(np.float32, copy=False)


def safe_product(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return norm(np.clip(a, 0, None) * np.clip(b, 0, None))


def build_base_channels(raw_zcyx: np.ndarray) -> dict[str, np.ndarray]:
    ch00 = norm(raw_zcyx[:, 0])
    bf = norm(raw_zcyx[:, 1])
    ch02 = norm(raw_zcyx[:, 2])
    fucci = np.maximum(ch00, ch02)
    fucci = ndi.gaussian_filter(fucci, sigma=(0.5, 0.8, 0.8)).astype(np.float32)

    bf_raw = raw_zcyx[:, 1].astype(np.float32, copy=False)
    absz = build_abs_zscore_smooth(bf_raw, blur_sigma_xy=20, z_anisotropy=1.44)
    abs_zscore = norm(absz["abs_zscore"])

    sd_score, _bf_norm, _bf_hp = focus_component_stack(
        bf_raw,
        method="local_sd",
        background_sigma=24,
        highpass_sigma=8,
        local_std_window=10,
        log_sigma=2.0,
    )
    sd_confidence = norm(relative_focus_confidence(sd_score))

    return {
        "fucci": fucci,
        "bf": bf,
        "absz": abs_zscore,
        "sd": sd_confidence,
    }


def build_artificial_channels(absz: np.ndarray, sd: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "absz": absz,
        "sd_w10": sd,
        "absz_sd_mean": norm((absz + sd) * 0.5),
        "absz_sd_geom": norm(np.sqrt(np.clip(absz, 0, None) * np.clip(sd, 0, None))),
        "absz_sd_product": safe_product(absz, sd),
        "absz_sd_min": norm(np.minimum(absz, sd)),
        "absz_sd_max": norm(np.maximum(absz, sd)),
        "absz2_sd_product": safe_product(absz * absz, sd),
        "absz_sd2_product": safe_product(absz, sd * sd),
    }


def build_variants(raw_zcyx: np.ndarray) -> dict[str, np.ndarray]:
    channels = build_base_channels(raw_zcyx)
    artificial = build_artificial_channels(channels["absz"], channels["sd"])

    variants: dict[str, np.ndarray] = {}
    for name, artificial_channel in artificial.items():
        variants[f"fucci_{name}"] = np.stack([channels["fucci"], artificial_channel], axis=1)
        variants[f"fucci_bf_{name}"] = np.stack([channels["fucci"], channels["bf"], artificial_channel], axis=1)
    return {name: value.astype(np.float32, copy=False) for name, value in variants.items()}


def colorize_labels(mask: np.ndarray) -> np.ndarray:
    mask = mask.astype(np.uint32, copy=False)
    out = np.zeros(mask.shape + (3,), dtype=np.uint8)
    ids = np.unique(mask)
    ids = ids[ids > 0]
    for label in ids:
        r = (label * 37 + 41) % 255
        g = (label * 67 + 89) % 255
        b = (label * 97 + 131) % 255
        out[mask == label] = (r, g, b)
    return out


def gray_rgb(x: np.ndarray) -> np.ndarray:
    u8 = (robust01(x) * 255).astype(np.uint8)
    return np.repeat(u8[..., None], 3, axis=2)


def save_grid_png(path: Path, brightfield: np.ndarray, masks: np.ndarray, z_labels: list[int]) -> None:
    try:
        from PIL import Image, ImageDraw
    except ImportError as exc:
        raise RuntimeError("Pillow is required for labeled preview PNG output") from exc

    tiles_top = [gray_rgb(plane) for plane in brightfield]
    tiles_bottom = [colorize_labels(mask) for mask in masks]
    tile_h, tile_w = tiles_top[0].shape[:2]
    label_h = 22
    gap = 8
    ncols = len(tiles_top)
    canvas_h = 2 * tile_h + label_h + gap
    canvas_w = ncols * tile_w + (ncols - 1) * gap
    canvas = np.full((canvas_h, canvas_w, 3), 255, dtype=np.uint8)

    for i, (top, bottom) in enumerate(zip(tiles_top, tiles_bottom)):
        x0 = i * (tile_w + gap)
        canvas[label_h : label_h + tile_h, x0 : x0 + tile_w] = top
        y1 = label_h + tile_h + gap
        canvas[y1 : y1 + tile_h, x0 : x0 + tile_w] = bottom

    img = Image.fromarray(canvas)
    draw = ImageDraw.Draw(img)
    for i, z_label in enumerate(z_labels):
        x0 = i * (tile_w + gap)
        draw.text((x0 + 4, 4), f"z={z_label}", fill=(0, 0, 0))
    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)


def eval_plane(model, plane_cyx: np.ndarray, min_size: int) -> np.ndarray:
    result = model.eval(plane_cyx, channel_axis=0, min_size=min_size)
    mask = result[0] if isinstance(result, tuple) else result
    return mask.astype(np.uint32, copy=False)


def append_summary(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "variant",
        "n_channels",
        "status",
        "raw_crop",
        "png_path",
        "objects_total",
        "objects_by_plane",
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
        description="Run 2D CPSAM PNG grids for FUCCI plus abs-zscore/SD artificial-channel combinations."
    )
    parser.add_argument("--raw-crop-dir", type=Path, default=DEFAULT_RAW_CROP_DIR)
    parser.add_argument("--raw-crop", action="append", type=Path, help="Specific raw ZCYX crop. May be repeated.")
    parser.add_argument("--pattern", default=f"*{RAW_CROP_SUFFIX}", help="Glob used with --raw-crop-dir.")
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument(
        "--z-labels",
        type=parse_int_list,
        default=parse_int_list("10,20,30,40,50,60,70"),
        help="1-based z labels to preview. z=70 maps to zero-based plane 69.",
    )
    parser.add_argument("--min-size", type=int, default=100)
    parser.add_argument("--limit-samples", type=int, default=0)
    parser.add_argument(
        "--variants",
        nargs="+",
        help="Optional subset of variant names. Example: fucci_absz_sd_geom fucci_bf_absz_sd_geom",
    )
    parser.add_argument("--resume", action="store_true", help="Skip variant/sample PNGs that already exist.")
    parser.add_argument("--cpu", action="store_true", help="Run CPSAM on CPU instead of the default GPU path.")
    args = parser.parse_args()

    raw_crops = args.raw_crop if args.raw_crop else discover_raw_crops(args.raw_crop_dir, args.pattern)
    if args.limit_samples:
        raw_crops = raw_crops[: args.limit_samples]
    if not raw_crops:
        raise ValueError("No raw crop TIFFs found.")

    z_indices = [z - 1 for z in args.z_labels]
    run_dir = args.run_dir
    png_dir = run_dir / "png"
    log_dir = run_dir / "logs"
    run_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    config = {
        "raw_crop_dir": str(args.raw_crop_dir),
        "raw_crops": [str(path) for path in raw_crops],
        "run_dir": str(run_dir),
        "z_labels_1_based": args.z_labels,
        "z_indices_0_based": z_indices,
        "min_size": args.min_size,
        "gpu": not args.cpu,
        "outputs": "PNG grids only plus TSV/config metadata; mask and input TIFFs are not saved.",
        "base_channels": {
            "fucci": "max(percentile-normalized ch00, percentile-normalized ch02), gaussian sigma=(0.5,0.8,0.8)",
            "bf": "percentile-normalized raw ch01 brightfield, only in fucci_bf_* variants",
            "absz": "abs(global z-score of raw ch01 brightfield), percentile-normalized",
            "sd_w10": "relative-focus confidence from local SD on high-pass ch01 brightfield, window=10",
        },
        "artificial_channels": {
            "absz": "absz baseline",
            "sd_w10": "SD detector baseline",
            "absz_sd_mean": "arithmetic mean of absz and sd_w10",
            "absz_sd_geom": "geometric mean sqrt(absz * sd_w10)",
            "absz_sd_product": "normalized product absz * sd_w10",
            "absz_sd_min": "pixelwise minimum consensus of absz and sd_w10",
            "absz_sd_max": "pixelwise maximum union of absz and sd_w10",
            "absz2_sd_product": "normalized product absz^2 * sd_w10",
            "absz_sd2_product": "normalized product absz * sd_w10^2",
        },
        "cpsam_input_axes": "CYX per selected z plane; segmentation is 2D only.",
        "variant_filter": args.variants or [],
    }
    (run_dir / "config.fucci_absz_sd_combo_cpsam_2d_grids.json").write_text(json.dumps(config, indent=2) + "\n")

    from cellpose import models

    model = models.CellposeModel(gpu=not args.cpu, pretrained_model="cpsam")
    summary_path = run_dir / "sample_variant_summary.tsv"

    for sample_idx, raw_crop_path in enumerate(raw_crops, start=1):
        sample_id = sample_id_from_raw_crop(raw_crop_path)
        try:
            raw = tf.imread(raw_crop_path)
            if raw.ndim != 4 or raw.shape[1] != 3:
                raise ValueError(f"Expected raw crop ZCYX with 3 channels, got {raw.shape}")
            if min(z_indices) < 0 or max(z_indices) >= raw.shape[0]:
                raise ValueError(f"Requested z labels {args.z_labels} outside available 1..{raw.shape[0]}")
            variants = build_variants(raw)
            variant_names = sorted(args.variants or variants)
            missing = [name for name in variant_names if name not in variants]
            if missing:
                raise ValueError(f"Unknown variant(s): {', '.join(missing)}")
        except Exception as exc:
            error_path = log_dir / f"{sample_id}_prepare_error.txt"
            error_path.write_text(f"{type(exc).__name__}: {exc}\n")
            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "variant": "",
                    "n_channels": "",
                    "status": "error",
                    "raw_crop": raw_crop_path,
                    "png_path": "",
                    "objects_total": "",
                    "objects_by_plane": "",
                    "elapsed_seconds": 0,
                    "error": str(exc),
                },
            )
            print(f"[{sample_idx}/{len(raw_crops)}] prepare error {sample_id}: {exc}", flush=True)
            continue

        for variant_idx, variant in enumerate(variant_names, start=1):
            png_path = png_dir / variant / f"{sample_id}_{variant}_z10-70_2x7_bf_mask.png"
            if args.resume and png_path.exists():
                append_summary(
                    summary_path,
                    {
                        "sample_id": sample_id,
                        "variant": variant,
                        "n_channels": variants[variant].shape[1],
                        "status": "skipped_existing",
                        "raw_crop": raw_crop_path,
                        "png_path": png_path,
                        "objects_total": "",
                        "objects_by_plane": "",
                        "elapsed_seconds": 0,
                        "error": "",
                    },
                )
                print(
                    f"[{sample_idx}/{len(raw_crops)} {variant_idx}/{len(variant_names)}] skipped {sample_id} {variant}",
                    flush=True,
                )
                continue

            start = time.time()
            try:
                img = variants[variant]
                print(
                    f"[{sample_idx}/{len(raw_crops)} {variant_idx}/{len(variant_names)}] "
                    f"segmenting {sample_id} {variant}",
                    flush=True,
                )
                masks = np.stack([eval_plane(model, img[z], args.min_size) for z in z_indices], axis=0)
                save_grid_png(png_path, raw[z_indices, 1], masks, args.z_labels)
                objects_by_plane = [int(mask.max()) for mask in masks]
                append_summary(
                    summary_path,
                    {
                        "sample_id": sample_id,
                        "variant": variant,
                        "n_channels": img.shape[1],
                        "status": "ok",
                        "raw_crop": raw_crop_path,
                        "png_path": png_path,
                        "objects_total": int(sum(objects_by_plane)),
                        "objects_by_plane": ",".join(str(v) for v in objects_by_plane),
                        "elapsed_seconds": round(time.time() - start, 3),
                        "error": "",
                    },
                )
                print(f"[{sample_idx}/{len(raw_crops)} {variant_idx}/{len(variant_names)}] wrote {png_path}", flush=True)
            except Exception as exc:
                error_path = log_dir / f"{sample_id}_{variant}_error.txt"
                error_path.write_text(f"{type(exc).__name__}: {exc}\n")
                append_summary(
                    summary_path,
                    {
                        "sample_id": sample_id,
                        "variant": variant,
                        "n_channels": variants[variant].shape[1],
                        "status": "error",
                        "raw_crop": raw_crop_path,
                        "png_path": png_path,
                        "objects_total": "",
                        "objects_by_plane": "",
                        "elapsed_seconds": round(time.time() - start, 3),
                        "error": str(exc),
                    },
                )
                print(f"[{sample_idx}/{len(raw_crops)} {variant_idx}/{len(variant_names)}] error: {exc}", flush=True)

    print(f"wrote summary: {summary_path}")


if __name__ == "__main__":
    main()
