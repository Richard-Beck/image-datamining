#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import time
import traceback
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
DEFAULT_RUN_DIR = Path("pipelines/fucci_3d_segmentation/runs/fucci_absz_sd2_product_cpsam_2d_fullstacks")
RAW_CROP_SUFFIX = "_raw3ch_crop256_ZCYX.tif"
INPUT_KIND = "fucci_absz_sd2_product"


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


def build_fucci_absz_sd2_product_input(raw_zcyx: np.ndarray) -> np.ndarray:
    ch00 = norm(raw_zcyx[:, 0])
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
    sd = norm(relative_focus_confidence(sd_score))
    artificial = norm(np.clip(abs_zscore, 0, None) * np.clip(sd, 0, None) * np.clip(sd, 0, None))
    return np.stack([fucci, artificial], axis=1).astype(np.float32)


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


def relabel_with_offset(mask: np.ndarray, offset: int) -> tuple[np.ndarray, dict[int, int], int]:
    labels = [int(label) for label in np.unique(mask) if label > 0]
    if not labels:
        return mask.astype(np.uint32, copy=True), {}, offset
    out = np.zeros(mask.shape, dtype=np.uint32)
    mapping: dict[int, int] = {}
    next_label = offset
    for local_label in labels:
        next_label += 1
        out[mask == local_label] = next_label
        mapping[local_label] = next_label
    return out, mapping, next_label


def object_rows(
    sample_id: str,
    z_index: int,
    z_label: int,
    local_mask: np.ndarray,
    global_mask: np.ndarray,
    mapping: dict[int, int],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for local_label, global_label in mapping.items():
        yy, xx = np.nonzero(local_mask == local_label)
        if yy.size == 0:
            continue
        rows.append(
            {
                "sample_id": sample_id,
                "z_index": z_index,
                "z_label_1_based": z_label,
                "local_label": local_label,
                "global_label": global_label,
                "area_px": int(yy.size),
                "bbox_ymin": int(yy.min()),
                "bbox_xmin": int(xx.min()),
                "bbox_ymax_exclusive": int(yy.max() + 1),
                "bbox_xmax_exclusive": int(xx.max() + 1),
                "centroid_y": float(yy.mean()),
                "centroid_x": float(xx.mean()),
            }
        )
    return rows


def write_object_table(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "sample_id",
        "z_index",
        "z_label_1_based",
        "local_label",
        "global_label",
        "area_px",
        "bbox_ymin",
        "bbox_xmin",
        "bbox_ymax_exclusive",
        "bbox_xmax_exclusive",
        "centroid_y",
        "centroid_x",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_mask_stack(path: Path, masks: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # ImageJ TIFF does not support uint32 labels. Standard tifffile metadata
    # preserves the axes while keeping the full label range for later stitching.
    tf.imwrite(path, masks.astype(np.uint32, copy=False), metadata={"axes": "ZYX"})


def append_summary(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "status",
        "raw_crop",
        "cellpose_input",
        "mask_path",
        "object_table",
        "png_path",
        "n_objects",
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
        description=(
            "Run full-stack 2D CPSAM for FUCCI fused plus product(abs z-score, SD squared), "
            "saving ZYX mask stacks for later 3D stitching."
        )
    )
    parser.add_argument("--raw-crop-dir", type=Path, default=DEFAULT_RAW_CROP_DIR)
    parser.add_argument("--raw-crop", action="append", type=Path, help="Specific raw ZCYX crop. May be repeated.")
    parser.add_argument("--pattern", default=f"*{RAW_CROP_SUFFIX}", help="Glob used with --raw-crop-dir.")
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument(
        "--z-labels",
        type=parse_int_list,
        default=parse_int_list("10,20,30,40,50,60,70"),
        help="1-based z labels for preview PNG. z=70 maps to zero-based plane 69.",
    )
    parser.add_argument("--min-size", type=int, default=100)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--resume", action="store_true", help="Skip samples whose mask, object table, and PNG already exist.")
    parser.add_argument("--save-inputs", action="store_true", help="Also save the prepared ZCYX CPSAM input TIFFs.")
    parser.add_argument("--cpu", action="store_true", help="Run CPSAM on CPU instead of the default GPU path.")
    args = parser.parse_args()

    raw_crops = args.raw_crop if args.raw_crop else discover_raw_crops(args.raw_crop_dir, args.pattern)
    if args.limit:
        raw_crops = raw_crops[: args.limit]
    if not raw_crops:
        raise ValueError("No raw crop TIFFs found.")

    z_indices = [z - 1 for z in args.z_labels]
    run_dir = args.run_dir
    input_dir = run_dir / "cellpose_inputs"
    mask_dir = run_dir / "masks"
    object_dir = run_dir / "objects"
    png_dir = run_dir / "png"
    log_dir = run_dir / "logs"
    run_dir.mkdir(parents=True, exist_ok=True)
    mask_dir.mkdir(parents=True, exist_ok=True)
    object_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    if args.save_inputs:
        input_dir.mkdir(parents=True, exist_ok=True)

    config = {
        "raw_crop_dir": str(args.raw_crop_dir),
        "raw_crops": [str(path) for path in raw_crops],
        "run_dir": str(run_dir),
        "input_kind": INPUT_KIND,
        "segmentation_mode": "2d_full_stack",
        "labeling": "2D masks are relabeled to globally unique labels across z within each sample.",
        "preprocessing": {
            "fucci": "max(percentile-normalized ch00, percentile-normalized ch02), gaussian sigma=(0.5,0.8,0.8)",
            "absz": "abs(global z-score of raw ch01 brightfield), percentile-normalized",
            "sd_w10": "relative-focus confidence from local SD on high-pass ch01 brightfield, window=10",
            "artificial": "percentile-normalized product absz * sd_w10^2",
            "cpsam_input_axes": "CYX per z plane; segmentation is independent 2D over all z slices",
        },
        "z_labels_1_based": args.z_labels,
        "z_indices_0_based": z_indices,
        "min_size": args.min_size,
        "gpu": not args.cpu,
        "save_inputs": args.save_inputs,
        "outputs": {
            "masks": "Full ZYX uint32 TIFF mask stack, globally unique 2D labels per sample.",
            "objects": "Per-slice object table with global labels and geometry for later 3D stitching.",
            "png": "2x7 preview sampled from the full 2D mask stack.",
        },
    }
    (run_dir / "config.fucci_absz_sd2_product_cpsam_2d_fullstacks.json").write_text(json.dumps(config, indent=2) + "\n")

    from cellpose import models

    model = models.CellposeModel(gpu=not args.cpu, pretrained_model="cpsam")
    summary_path = run_dir / "sample_summary.tsv"

    for i, raw_crop_path in enumerate(raw_crops, start=1):
        sample_id = sample_id_from_raw_crop(raw_crop_path)
        mask_path = mask_dir / f"{sample_id}_{INPUT_KIND}_2d_fullstack_masks.tif"
        object_path = object_dir / f"{sample_id}_{INPUT_KIND}_2d_objects.tsv"
        png_path = png_dir / f"{sample_id}_{INPUT_KIND}_2d_fullstack_z10-70_2x7_bf_mask.png"
        input_path = input_dir / f"{sample_id}_{INPUT_KIND}_ZCYX.tif"
        input_for_summary = str(input_path) if args.save_inputs else ""

        if args.resume and mask_path.exists() and object_path.exists() and png_path.exists():
            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "status": "skipped_existing",
                    "raw_crop": raw_crop_path,
                    "cellpose_input": input_for_summary,
                    "mask_path": mask_path,
                    "object_table": object_path,
                    "png_path": png_path,
                    "n_objects": "",
                    "elapsed_seconds": 0,
                    "error": "",
                },
            )
            print(f"[{i}/{len(raw_crops)}] skipped existing {sample_id}", flush=True)
            continue

        start = time.time()
        try:
            raw = tf.imread(raw_crop_path)
            if raw.ndim != 4 or raw.shape[1] != 3:
                raise ValueError(f"Expected raw crop ZCYX with 3 channels, got {raw.shape}")
            if min(z_indices) < 0 or max(z_indices) >= raw.shape[0]:
                raise ValueError(f"Requested z labels {args.z_labels} outside available 1..{raw.shape[0]}")

            img = build_fucci_absz_sd2_product_input(raw)
            if args.save_inputs:
                tf.imwrite(input_path, img, imagej=True, metadata={"axes": "ZCYX"})

            print(f"[{i}/{len(raw_crops)}] segmenting all z in 2D {sample_id}", flush=True)
            masks = np.zeros((img.shape[0], img.shape[2], img.shape[3]), dtype=np.uint32)
            rows: list[dict[str, object]] = []
            label_offset = 0
            for z in range(img.shape[0]):
                local_mask = eval_plane(model, img[z], args.min_size)
                global_mask, mapping, label_offset = relabel_with_offset(local_mask, label_offset)
                masks[z] = global_mask
                rows.extend(object_rows(sample_id, z, z + 1, local_mask, global_mask, mapping))
                print(f"[{i}/{len(raw_crops)}] {sample_id} z={z + 1}/{img.shape[0]}", flush=True)

            write_mask_stack(mask_path, masks)
            write_object_table(object_path, rows)
            save_grid_png(png_path, raw[z_indices, 1], masks[z_indices], args.z_labels)

            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "status": "ok",
                    "raw_crop": raw_crop_path,
                    "cellpose_input": input_for_summary,
                    "mask_path": mask_path,
                    "object_table": object_path,
                    "png_path": png_path,
                    "n_objects": len(rows),
                    "elapsed_seconds": round(time.time() - start, 3),
                    "error": "",
                },
            )
            print(f"[{i}/{len(raw_crops)}] wrote {mask_path}, {object_path}, and {png_path}", flush=True)
        except Exception:
            error_path = log_dir / f"{sample_id}_error.txt"
            error_path.write_text(traceback.format_exc())
            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "status": "error",
                    "raw_crop": raw_crop_path,
                    "cellpose_input": input_for_summary,
                    "mask_path": mask_path,
                    "object_table": object_path,
                    "png_path": png_path,
                    "n_objects": "",
                    "elapsed_seconds": round(time.time() - start, 3),
                    "error": str(error_path),
                },
            )
            print(f"[{i}/{len(raw_crops)}] error {sample_id}; see {error_path}", flush=True)

    print(f"wrote summary: {summary_path}")


if __name__ == "__main__":
    main()
