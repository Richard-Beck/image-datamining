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
from image_datamining.fucci_focus_channel.focus import build_abs_zscore_smooth, robust01


DEFAULT_RAW_CROP_DIR = Path(
    "pipelines/fucci_3d_segmentation/runs/jackson_fucci_crop256_bf_fucci_fused_3d/raw_crops"
)
DEFAULT_RUN_DIR = Path("pipelines/fucci_3d_segmentation/runs/fucci_absz_cpsam_2d_grids")
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


def build_fucci_absz_input(raw_zcyx: np.ndarray) -> np.ndarray:
    ch00 = norm(raw_zcyx[:, 0])
    ch02 = norm(raw_zcyx[:, 2])
    fucci = np.maximum(ch00, ch02)
    fucci = ndi.gaussian_filter(fucci, sigma=(0.5, 0.8, 0.8)).astype(np.float32)

    bf_raw = raw_zcyx[:, 1].astype(np.float32, copy=False)
    absz = build_abs_zscore_smooth(bf_raw, blur_sigma_xy=20, z_anisotropy=1.44)
    abs_zscore = norm(absz["abs_zscore"])
    return np.stack([fucci, abs_zscore], axis=1).astype(np.float32)


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
        "status",
        "raw_crop",
        "input_path",
        "mask_path",
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
        description="Run 2D CPSAM grids for FUCCI fused + brightfield abs z-score across raw crop TIFFs."
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
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--resume", action="store_true", help="Skip samples whose mask and PNG already exist.")
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
    png_dir = run_dir / "png"
    log_dir = run_dir / "logs"
    run_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    config = {
        "raw_crop_dir": str(args.raw_crop_dir),
        "raw_crops": [str(path) for path in raw_crops],
        "run_dir": str(run_dir),
        "input_kind": "fucci_absz",
        "preprocessing": {
            "fucci": "max(percentile-normalized ch00, percentile-normalized ch02), gaussian sigma=(0.5,0.8,0.8)",
            "absz": "abs(global z-score of raw ch01 brightfield), percentile-normalized",
            "cpsam_input_axes": "CYX per selected z plane",
        },
        "z_labels_1_based": args.z_labels,
        "z_indices_0_based": z_indices,
        "min_size": args.min_size,
        "gpu": not args.cpu,
    }
    (run_dir / "config.fucci_absz_cpsam_2d_grids.json").write_text(json.dumps(config, indent=2) + "\n")

    from cellpose import models

    model = models.CellposeModel(gpu=not args.cpu, pretrained_model="cpsam")
    summary_path = run_dir / "sample_summary.tsv"

    for i, raw_crop_path in enumerate(raw_crops, start=1):
        sample_id = sample_id_from_raw_crop(raw_crop_path)
        input_path = input_dir / f"{sample_id}_fucci_absz_z10-70_ZCYX.tif"
        mask_path = mask_dir / f"{sample_id}_fucci_absz_z10-70_2d_cpsam_masks.tif"
        png_path = png_dir / f"{sample_id}_fucci_absz_z10-70_2x7_bf_mask.png"

        if args.resume and mask_path.exists() and png_path.exists():
            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "status": "skipped_existing",
                    "raw_crop": raw_crop_path,
                    "input_path": input_path,
                    "mask_path": mask_path,
                    "png_path": png_path,
                    "objects_total": "",
                    "objects_by_plane": "",
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

            img = build_fucci_absz_input(raw)
            input_dir.mkdir(parents=True, exist_ok=True)
            mask_dir.mkdir(parents=True, exist_ok=True)
            tf.imwrite(input_path, img[z_indices], imagej=True, metadata={"axes": "ZCYX"})

            print(f"[{i}/{len(raw_crops)}] segmenting {sample_id}", flush=True)
            masks = np.stack([eval_plane(model, img[z], args.min_size) for z in z_indices], axis=0)
            tf.imwrite(mask_path, masks.astype(np.uint32))
            save_grid_png(png_path, raw[z_indices, 1], masks, args.z_labels)

            objects_by_plane = [int(mask.max()) for mask in masks]
            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "status": "ok",
                    "raw_crop": raw_crop_path,
                    "input_path": input_path,
                    "mask_path": mask_path,
                    "png_path": png_path,
                    "objects_total": int(sum(objects_by_plane)),
                    "objects_by_plane": ",".join(str(v) for v in objects_by_plane),
                    "elapsed_seconds": round(time.time() - start, 3),
                    "error": "",
                },
            )
            print(f"[{i}/{len(raw_crops)}] wrote {png_path}", flush=True)
        except Exception as exc:
            error_path = log_dir / f"{sample_id}_error.txt"
            error_path.write_text(f"{type(exc).__name__}: {exc}\n")
            append_summary(
                summary_path,
                {
                    "sample_id": sample_id,
                    "status": "error",
                    "raw_crop": raw_crop_path,
                    "input_path": input_path,
                    "mask_path": mask_path,
                    "png_path": png_path,
                    "objects_total": "",
                    "objects_by_plane": "",
                    "elapsed_seconds": round(time.time() - start, 3),
                    "error": str(exc),
                },
            )
            print(f"[{i}/{len(raw_crops)}] error {sample_id}: {exc}", flush=True)

    print(f"wrote summary: {summary_path}")


if __name__ == "__main__":
    main()
