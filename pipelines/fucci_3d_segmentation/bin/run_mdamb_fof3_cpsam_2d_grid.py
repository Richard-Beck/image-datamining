#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
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
    build_sd_log_combo,
    focus_component_stack,
    relative_focus_confidence,
    robust01,
)


SAMPLE_ID = "MDAMB_240904_MDAMB_FUCCI_MDAMB231_30K_FoF3_fluorescent_nucleus"
DEFAULT_RAW_CROP = Path(
    "pipelines/fucci_3d_segmentation/runs/jackson_fucci_crop256_bf_fucci_fused_3d/raw_crops/"
    f"{SAMPLE_ID}_raw3ch_crop256_ZCYX.tif"
)
DEFAULT_RUN_DIR = Path("pipelines/fucci_3d_segmentation/runs/mdamb_fof3_cpsam_2d_input_combos")


def parse_int_list(text: str) -> list[int]:
    values = [int(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("Expected at least one comma-separated integer")
    return values


def norm(x: np.ndarray) -> np.ndarray:
    return percentile_normalize(x).astype(np.float32, copy=False)


def colorize_labels(mask: np.ndarray) -> np.ndarray:
    mask = mask.astype(np.uint32, copy=False)
    out = np.zeros(mask.shape + (3,), dtype=np.uint8)
    ids = np.unique(mask)
    ids = ids[ids > 0]
    for label in ids:
        # Deterministic compact color hash; avoids storing a large LUT.
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


def build_variants(raw_zcyx: np.ndarray) -> dict[str, np.ndarray]:
    ch00 = norm(raw_zcyx[:, 0])
    bf = norm(raw_zcyx[:, 1])
    ch02 = norm(raw_zcyx[:, 2])
    fucci = np.maximum(ch00, ch02)
    fucci_smooth = ndi.gaussian_filter(fucci, sigma=(0.5, 0.8, 0.8)).astype(np.float32)

    bf_raw = raw_zcyx[:, 1].astype(np.float32, copy=False)
    absz = build_abs_zscore_smooth(bf_raw, blur_sigma_xy=20, z_anisotropy=1.44)
    abs_zscore = norm(absz["abs_zscore"])
    absz_product = norm(absz["abs_zscore"] * absz["smooth"])

    sd_score, _bf_norm, _bf_hp = focus_component_stack(
        bf_raw,
        method="local_sd",
        background_sigma=24,
        highpass_sigma=8,
        local_std_window=10,
        log_sigma=2.0,
    )
    log_score, _bf_norm, _bf_hp = focus_component_stack(
        bf_raw,
        method="log",
        background_sigma=24,
        highpass_sigma=8,
        local_std_window=10,
        log_sigma=2.0,
    )
    sd_conf = norm(relative_focus_confidence(sd_score))
    log_conf = norm(relative_focus_confidence(log_score))
    combo = build_sd_log_combo(
        bf_raw,
        sd_window=10,
        log_sigma=2.0,
        blur_sigma_xy=20,
        z_anisotropy=1.44,
        background_sigma=24,
        highpass_sigma=8,
    )

    scalar_channels = {
        "absz": abs_zscore,
        "absz_product": absz_product,
        "sd_w10": sd_conf,
        "log2": log_conf,
        "sd_log_combined": norm(combo["combined"]),
        "sd_log_product": norm(combo["product"]),
        "sd_log_otsu_product": norm(combo["otsu_product"]),
    }

    variants: dict[str, np.ndarray] = {
        "raw_ch00": ch00[:, None],
        "raw_ch01_bf": bf[:, None],
        "raw_ch02": ch02[:, None],
        "raw_ch00_ch02": np.stack([ch00, ch02], axis=1),
        "raw_all_ch00_ch01_ch02": np.stack([ch00, bf, ch02], axis=1),
        "fucci": fucci_smooth[:, None],
        "fucci_bf": np.stack([fucci_smooth, bf], axis=1),
    }
    for name, extra in scalar_channels.items():
        variants[f"fucci_{name}"] = np.stack([fucci_smooth, extra], axis=1)
        variants[f"fucci_bf_{name}"] = np.stack([fucci_smooth, bf, extra], axis=1)
    return variants


def eval_plane(model, plane_cyx: np.ndarray, min_size: int) -> np.ndarray:
    result = model.eval(plane_cyx, channel_axis=0, min_size=min_size)
    mask = result[0] if isinstance(result, tuple) else result
    return mask.astype(np.uint32, copy=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run 2D CPSAM input-combo sweep for the MDAMB FoF3 crop.")
    parser.add_argument("--raw-crop", type=Path, default=DEFAULT_RAW_CROP)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument(
        "--z-labels",
        type=parse_int_list,
        default=parse_int_list("10,20,30,40,50,60,70"),
        help="1-based z labels to preview. z=70 maps to zero-based plane 69.",
    )
    parser.add_argument("--min-size", type=int, default=100)
    parser.add_argument("--cpu", action="store_true", help="Run CPSAM on CPU instead of the default GPU path.")
    parser.add_argument(
        "--limit-variants",
        type=parse_int_list,
        help="Debug helper: run only variants by zero-based index in the sorted variant list.",
    )
    args = parser.parse_args()

    raw = tf.imread(args.raw_crop)
    if raw.ndim != 4 or raw.shape[1] != 3:
        raise ValueError(f"Expected raw crop ZCYX with 3 channels, got {raw.shape}")

    z_indices = [z - 1 for z in args.z_labels]
    if min(z_indices) < 0 or max(z_indices) >= raw.shape[0]:
        raise ValueError(f"Requested z labels {args.z_labels} outside available 1..{raw.shape[0]}")

    run_dir = args.run_dir
    input_dir = run_dir / "cellpose_inputs"
    mask_dir = run_dir / "masks"
    png_dir = run_dir / "png"
    run_dir.mkdir(parents=True, exist_ok=True)
    config = {
        "sample_id": SAMPLE_ID,
        "raw_crop": str(args.raw_crop),
        "run_dir": str(run_dir),
        "z_labels_1_based": args.z_labels,
        "z_indices_0_based": z_indices,
        "min_size": args.min_size,
        "gpu": not args.cpu,
        "note": "2D CPSAM sweep; z labels are 1-based for the requested 10,20,...,70 planes.",
    }
    (run_dir / "config.mdamb_fof3_cpsam_2d_grid.json").write_text(json.dumps(config, indent=2) + "\n")

    variants = build_variants(raw)
    variant_names = sorted(variants)
    if args.limit_variants is not None:
        variant_names = [variant_names[i] for i in args.limit_variants]

    from cellpose import models

    model = models.CellposeModel(gpu=not args.cpu, pretrained_model="cpsam")
    bf_preview = raw[z_indices, 1]

    summary_path = run_dir / "variant_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "variant",
                "n_channels",
                "input_path",
                "mask_path",
                "png_path",
                "objects_total",
                "objects_by_plane",
                "elapsed_seconds",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for i, name in enumerate(variant_names, start=1):
            start = time.time()
            img = variants[name].astype(np.float32, copy=False)
            input_path = input_dir / f"{SAMPLE_ID}_{name}_z10-70_CYX.tif"
            mask_path = mask_dir / f"{SAMPLE_ID}_{name}_z10-70_2d_cpsam_masks.tif"
            png_path = png_dir / f"{SAMPLE_ID}_{name}_z10-70_2x7_bf_mask.png"

            input_dir.mkdir(parents=True, exist_ok=True)
            mask_dir.mkdir(parents=True, exist_ok=True)
            tf.imwrite(input_path, img[z_indices], imagej=True, metadata={"axes": "ZCYX"})

            print(f"[{i}/{len(variant_names)}] {name}: {img.shape[1]} channel(s)", flush=True)
            masks = np.stack([eval_plane(model, img[z], args.min_size) for z in z_indices], axis=0)
            tf.imwrite(mask_path, masks.astype(np.uint32))
            save_grid_png(png_path, bf_preview, masks, args.z_labels)

            objects_by_plane = [int(mask.max()) for mask in masks]
            writer.writerow(
                {
                    "variant": name,
                    "n_channels": img.shape[1],
                    "input_path": input_path,
                    "mask_path": mask_path,
                    "png_path": png_path,
                    "objects_total": int(sum(objects_by_plane)),
                    "objects_by_plane": ",".join(str(v) for v in objects_by_plane),
                    "elapsed_seconds": round(time.time() - start, 3),
                }
            )
            handle.flush()
            print(f"[{i}/{len(variant_names)}] wrote {png_path}", flush=True)

    print(f"wrote summary: {summary_path}")


if __name__ == "__main__":
    main()
