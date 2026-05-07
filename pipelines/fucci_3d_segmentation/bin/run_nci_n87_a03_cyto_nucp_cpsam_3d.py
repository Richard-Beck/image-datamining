#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import time
import traceback
from pathlib import Path

import numpy as np

from _bootstrap import add_src_to_path

add_src_to_path()

import tifffile as tf

from image_datamining.fucci_3d_segmentation.preprocess import centered_xy_crop, percentile_normalize
from image_datamining.fucci_focus_channel.focus import robust01


DEFAULT_SOURCE_DIR = Path(
    "/share/lab_crd/lab_crd/MeasuringFitnessPerClone/data/GastricCancerCLs/"
    "3Dbrightfield/NCI-N87/A03_allenModel/FoF4_220228_fluorescent.cytoplasm"
)
DEFAULT_RUN_DIR = Path("pipelines/fucci_3d_segmentation/runs/nci_n87_a03_fof4_cyto_t_nuc_p_cpsam_3d")
SAMPLE_ID = "NCI-N87_FoF4_220228_cytoplasm_t_nucleus_p"


def parse_int_list(text: str) -> list[int]:
    values = [int(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("Expected at least one comma-separated integer")
    return values


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


def overlay_labels_on_gray(gray_plane: np.ndarray, mask: np.ndarray, alpha: float = 0.45) -> np.ndarray:
    base = gray_rgb(gray_plane).astype(np.float32)
    colors = colorize_labels(mask).astype(np.float32)
    keep = mask > 0
    out = base.copy()
    out[keep] = (1.0 - alpha) * base[keep] + alpha * colors[keep]
    return np.clip(out, 0, 255).astype(np.uint8)


def save_grid_png(
    path: Path,
    brightfield: np.ndarray,
    masks: np.ndarray,
    z_labels: list[int],
    overlay: bool,
) -> None:
    try:
        from PIL import Image, ImageDraw
    except ImportError as exc:
        raise RuntimeError("Pillow is required for labeled preview PNG output") from exc

    tiles_top = [gray_rgb(plane) for plane in brightfield]
    if overlay:
        tiles_bottom = [overlay_labels_on_gray(plane, mask) for plane, mask in zip(brightfield, masks)]
    else:
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


def summarize_masks(masks: np.ndarray) -> dict[str, object]:
    ids, counts = np.unique(masks[masks > 0], return_counts=True)
    summary: dict[str, object] = {"objects": int(len(ids))}
    if len(ids):
        summary["max_label"] = int(ids.max())
        summary["voxel_min"] = int(counts.min())
        summary["voxel_median"] = float(np.median(counts))
        summary["voxel_max"] = int(counts.max())
    else:
        summary["max_label"] = 0
        summary["voxel_min"] = 0
        summary["voxel_median"] = 0
        summary["voxel_max"] = 0
    return summary


def run_cpsam_3d(img_zcyx: np.ndarray, anisotropy: float, min_size: int, gpu: bool) -> np.ndarray:
    from cellpose import models

    model = models.CellposeModel(gpu=gpu, pretrained_model="cpsam")
    result = model.eval(
        img_zcyx,
        do_3D=True,
        z_axis=0,
        channel_axis=1,
        anisotropy=anisotropy,
        min_size=min_size,
    )
    masks = result[0] if isinstance(result, tuple) else result
    return masks.astype(np.uint32, copy=False)


def write_summary(path: Path, row: dict[str, object]) -> None:
    fieldnames = [
        "sample_id",
        "status",
        "source_dir",
        "raw_crop",
        "cellpose_input",
        "mask_path",
        "grid_png",
        "overlay_png",
        "objects",
        "max_label",
        "voxel_min",
        "voxel_median",
        "voxel_max",
        "elapsed_seconds",
        "error",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Crop NCI-N87 A03 cytoplasm.t plus nucleus.p, run 3D CellposeSAM, "
            "and write masks plus 2x7 brightfield/mask review PNGs."
        )
    )
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--crop-size", type=int, default=256)
    parser.add_argument("--z-labels", type=parse_int_list, default=parse_int_list("10,20,30,40,50,60,70"))
    parser.add_argument("--anisotropy", type=float, default=1.44)
    parser.add_argument("--min-size", type=int, default=100)
    parser.add_argument("--cpu", action="store_true", help="Run CPSAM on CPU instead of the default GPU path.")
    parser.add_argument("--prepare-only", action="store_true", help="Write crops and config without running CPSAM.")
    args = parser.parse_args()

    start = time.time()
    run_dir = args.run_dir
    crop_dir = run_dir / "raw_crops"
    input_dir = run_dir / "cellpose_inputs"
    mask_dir = run_dir / "masks"
    png_dir = run_dir / "png"
    log_dir = run_dir / "logs"
    for directory in (crop_dir, input_dir, mask_dir, png_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)

    z_indices = [z - 1 for z in args.z_labels]
    raw_crop_path = crop_dir / f"{SAMPLE_ID}_raw3ch_crop{args.crop_size}_ZCYX.tif"
    input_path = input_dir / f"{SAMPLE_ID}_cytoT_nucP_crop{args.crop_size}_ZCYX.tif"
    mask_path = mask_dir / f"{SAMPLE_ID}_cytoT_nucP_crop{args.crop_size}_3d_cpsam_masks.tif"
    grid_png_path = png_dir / f"{SAMPLE_ID}_cytoT_nucP_crop{args.crop_size}_3d_z10-70_2x7_bf_mask.png"
    overlay_png_path = png_dir / f"{SAMPLE_ID}_cytoT_nucP_crop{args.crop_size}_3d_z10-70_2x7_bf_mask_overlay.png"
    summary_path = run_dir / "sample_summary.tsv"

    config = {
        "source_dir": str(args.source_dir),
        "sample_id": SAMPLE_ID,
        "channels": {
            "cellpose_input_ch0": "cytoplasm.t.tif observed target",
            "cellpose_input_ch1": "nucleus.p.tif Allen/fnet prediction",
            "raw_crop_ch0": "cytoplasm.t.tif observed target",
            "raw_crop_ch1": "nucleus.p.tif Allen/fnet prediction",
            "raw_crop_ch2": "tifs/0_signal.tif brightfield signal",
        },
        "crop_size": args.crop_size,
        "z_labels_1_based": args.z_labels,
        "z_indices_0_based": z_indices,
        "anisotropy": args.anisotropy,
        "min_size": args.min_size,
        "gpu": not args.cpu,
        "prepare_only": args.prepare_only,
        "axes": {
            "raw_crop": "ZCYX",
            "cellpose_input": "ZCYX",
            "mask": "ZYX",
        },
    }
    (run_dir / "config.nci_n87_a03_cyto_t_nuc_p_cpsam_3d.json").write_text(json.dumps(config, indent=2) + "\n")

    try:
        cyt = tf.imread(args.source_dir / "cytoplasm.t.tif")
        nuc = tf.imread(args.source_dir / "nucleus.p.tif")
        bf = tf.imread(args.source_dir / "tifs" / "0_signal.tif")
        if cyt.shape != nuc.shape or cyt.shape != bf.shape:
            raise ValueError(f"Channel shapes differ: cyt={cyt.shape}, nuc={nuc.shape}, bf={bf.shape}")
        if cyt.ndim != 3:
            raise ValueError(f"Expected 3D ZYX TIFF inputs, got cytoplasm shape {cyt.shape}")
        if min(z_indices) < 0 or max(z_indices) >= cyt.shape[0]:
            raise ValueError(f"Requested z labels {args.z_labels} outside available 1..{cyt.shape[0]}")

        cyt_crop = centered_xy_crop(cyt, args.crop_size)
        nuc_crop = centered_xy_crop(nuc, args.crop_size)
        bf_crop = centered_xy_crop(bf, args.crop_size)

        raw_crop = np.stack([cyt_crop, nuc_crop, bf_crop], axis=1).astype(np.float32)
        tf.imwrite(raw_crop_path, raw_crop, imagej=True, metadata={"axes": "ZCYX"})

        img = np.stack(
            [
                percentile_normalize(cyt_crop),
                percentile_normalize(nuc_crop),
            ],
            axis=1,
        ).astype(np.float32)
        tf.imwrite(input_path, img, imagej=True, metadata={"axes": "ZCYX"})

        if args.prepare_only:
            write_summary(
                summary_path,
                {
                    "sample_id": SAMPLE_ID,
                    "status": "prepared",
                    "source_dir": args.source_dir,
                    "raw_crop": raw_crop_path,
                    "cellpose_input": input_path,
                    "mask_path": "",
                    "grid_png": "",
                    "overlay_png": "",
                    "objects": "",
                    "max_label": "",
                    "voxel_min": "",
                    "voxel_median": "",
                    "voxel_max": "",
                    "elapsed_seconds": round(time.time() - start, 3),
                    "error": "",
                },
            )
            print(f"prepared raw crop: {raw_crop_path}")
            print(f"prepared CPSAM input: {input_path}")
            print(f"wrote summary: {summary_path}")
            return

        print(f"segmenting {SAMPLE_ID} with 3D CPSAM", flush=True)
        masks = run_cpsam_3d(img, anisotropy=args.anisotropy, min_size=args.min_size, gpu=not args.cpu)
        tf.imwrite(mask_path, masks.astype(np.uint32), metadata={"axes": "ZYX"})
        save_grid_png(grid_png_path, bf_crop[z_indices], masks[z_indices], args.z_labels, overlay=False)
        save_grid_png(overlay_png_path, bf_crop[z_indices], masks[z_indices], args.z_labels, overlay=True)

        mask_summary = summarize_masks(masks)
        write_summary(
            summary_path,
            {
                "sample_id": SAMPLE_ID,
                "status": "ok",
                "source_dir": args.source_dir,
                "raw_crop": raw_crop_path,
                "cellpose_input": input_path,
                "mask_path": mask_path,
                "grid_png": grid_png_path,
                "overlay_png": overlay_png_path,
                "elapsed_seconds": round(time.time() - start, 3),
                "error": "",
                **mask_summary,
            },
        )
        print(f"wrote raw crop: {raw_crop_path}")
        print(f"wrote CPSAM input: {input_path}")
        print(f"wrote masks: {mask_path}")
        print(f"wrote grid PNG: {grid_png_path}")
        print(f"wrote overlay PNG: {overlay_png_path}")
        print(f"wrote summary: {summary_path}")
        print(json.dumps(mask_summary, indent=2))
    except Exception:
        error_path = log_dir / f"{SAMPLE_ID}_error.txt"
        error_path.write_text(traceback.format_exc())
        write_summary(
            summary_path,
            {
                "sample_id": SAMPLE_ID,
                "status": "error",
                "source_dir": args.source_dir,
                "raw_crop": raw_crop_path,
                "cellpose_input": input_path,
                "mask_path": mask_path,
                "grid_png": grid_png_path,
                "overlay_png": overlay_png_path,
                "objects": "",
                "max_label": "",
                "voxel_min": "",
                "voxel_median": "",
                "voxel_max": "",
                "elapsed_seconds": round(time.time() - start, 3),
                "error": error_path,
            },
        )
        print(f"error; see {error_path}", flush=True)
        raise


if __name__ == "__main__":
    main()

