#!/usr/bin/env python3
"""Run a small CellposeSAM sample on K00 registered image frames.

This is intended as a quick segmentation-quality smoke test for the registered
Images_40Frames TIFF stacks used by the K00 tracking-overlay workflow.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import imageio.v3 as iio
import numpy as np
import tifffile as tf
from cellpose import models


DEFAULT_IMAGES_DIR = Path(
    "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/"
    "SUM-159/K00_GemcitabineExposure_033023/"
    "New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/"
    "Final_Tracking_analysis/Images_40Frames"
)
DEFAULT_OUT_DIR = Path(
    "analyses/K00_GemcitabineExposure_033023/cpsam_frame_sample"
)
SITE_RE = re.compile(r"^(?P<well>[A-H]\d{1,2})_(?P<position>\d+)\.tiff?$", re.IGNORECASE)


@dataclass(frozen=True)
class FrameSample:
    image_path: Path
    well: str
    position: str
    frame_index0: int
    frame_label: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sample frames from K00 Images_40Frames, run CellposeSAM, and write "
            "mask TIFFs plus PNG outline overlays."
        )
    )
    parser.add_argument("--images-dir", type=Path, default=DEFAULT_IMAGES_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--n-frames", type=int, default=24)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument(
        "--frame-base",
        type=int,
        default=0,
        help="Frame numbering base to report in filenames/CSV. K00 tracking uses 0.",
    )
    parser.add_argument(
        "--max-sites",
        type=int,
        default=24,
        help="Maximum number of TIFF stacks to sample before taking frames.",
    )
    parser.add_argument(
        "--channel-mode",
        choices=("rgb", "gray", "red", "green", "blue"),
        default="rgb",
        help="Input representation passed to CPSAM.",
    )
    parser.add_argument(
        "--diameter",
        type=float,
        default=None,
        help="Optional Cellpose diameter. Omit for CPSAM automatic sizing.",
    )
    parser.add_argument("--flow-threshold", type=float, default=0.4)
    parser.add_argument("--cellprob-threshold", type=float, default=0.0)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument(
        "--cpu",
        action="store_true",
        help="Run on CPU. GPU is the default and is strongly preferred for CPSAM.",
    )
    parser.add_argument(
        "--no-montage",
        action="store_true",
        help="Skip writing the combined overlay contact sheet.",
    )
    return parser.parse_args()


def site_sort_key(path: Path) -> tuple[str, int, str]:
    match = SITE_RE.match(path.name)
    if match:
        return (match.group("well"), int(match.group("position")), path.name)
    return ("", 0, path.name)


def list_sites(images_dir: Path) -> list[Path]:
    sites = sorted(images_dir.glob("*.tiff"), key=site_sort_key)
    if not sites:
        sites = sorted(images_dir.glob("*.tif"), key=site_sort_key)
    if not sites:
        raise FileNotFoundError(f"No TIFF files found in {images_dir}")
    return sites


def site_metadata(path: Path) -> tuple[str, str]:
    match = SITE_RE.match(path.name)
    if match:
        return match.group("well").upper(), match.group("position")
    return path.stem, "NA"


def stack_frame_count(path: Path) -> int:
    with tf.TiffFile(path) as tif:
        if tif.series:
            shape = tif.series[0].shape
            axes = tif.series[0].axes
            if "I" in axes:
                return int(shape[axes.index("I")])
            if len(tif.pages) > 1:
                return len(tif.pages)
        return len(tif.pages)


def choose_samples(sites: list[Path], n_frames: int, max_sites: int, seed: int, frame_base: int) -> list[FrameSample]:
    rng = np.random.default_rng(seed)
    n_sites = min(max_sites, n_frames, len(sites))
    site_indices = np.linspace(0, len(sites) - 1, n_sites, dtype=int)
    chosen_sites = [sites[i] for i in site_indices]

    samples: list[FrameSample] = []
    per_site = [n_frames // n_sites] * n_sites
    for i in range(n_frames % n_sites):
        per_site[i] += 1

    for path, count in zip(chosen_sites, per_site):
        n_stack_frames = stack_frame_count(path)
        if count >= n_stack_frames:
            frame_indices = np.arange(n_stack_frames)
        else:
            frame_indices = np.linspace(0, n_stack_frames - 1, count, dtype=int)
            # Move off exact endpoints for a less edge-heavy sample when possible.
            if n_stack_frames > 4 and count == 1:
                frame_indices = np.array([int(rng.integers(1, n_stack_frames - 1))])
        well, position = site_metadata(path)
        for idx0 in frame_indices:
            samples.append(
                FrameSample(
                    image_path=path,
                    well=well,
                    position=position,
                    frame_index0=int(idx0),
                    frame_label=int(idx0) + frame_base,
                )
            )
    return samples[:n_frames]


def read_frame(sample: FrameSample) -> np.ndarray:
    with tf.TiffFile(sample.image_path) as tif:
        frame = tif.pages[sample.frame_index0].asarray()
    if frame.ndim == 2:
        frame = np.repeat(frame[:, :, None], 3, axis=2)
    if frame.ndim != 3:
        raise ValueError(f"Expected 2D or RGB frame from {sample.image_path}, got {frame.shape}")
    if frame.shape[2] > 3:
        frame = frame[:, :, :3]
    return frame


def prepare_model_input(frame: np.ndarray, channel_mode: str) -> np.ndarray:
    if channel_mode == "rgb":
        return frame
    channel_index = {"red": 0, "green": 1, "blue": 2}.get(channel_mode)
    if channel_index is None:
        gray = np.mean(frame[:, :, :3].astype(np.float32), axis=2)
    else:
        gray = frame[:, :, channel_index].astype(np.float32)
    return gray


def normalize_u8(frame: np.ndarray) -> np.ndarray:
    arr = frame.astype(np.float32)
    lo, hi = np.nanpercentile(arr, [1, 99.5])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
    arr = np.clip((arr - lo) / max(hi - lo, 1e-6), 0, 1)
    return (arr * 255).astype(np.uint8)


def label_boundaries(labels: np.ndarray) -> np.ndarray:
    labels = np.asarray(labels)
    boundaries = np.zeros(labels.shape, dtype=bool)
    boundaries[1:, :] |= labels[1:, :] != labels[:-1, :]
    boundaries[:-1, :] |= labels[:-1, :] != labels[1:, :]
    boundaries[:, 1:] |= labels[:, 1:] != labels[:, :-1]
    boundaries[:, :-1] |= labels[:, :-1] != labels[:, 1:]
    return boundaries & (labels > 0)


def outline_overlay(frame: np.ndarray, masks: np.ndarray) -> np.ndarray:
    if frame.ndim == 2:
        rgb = np.repeat(normalize_u8(frame)[:, :, None], 3, axis=2)
    else:
        rgb = normalize_u8(frame[:, :, :3])
    outline = label_boundaries(masks)
    rgb[outline, :] = np.array([255, 35, 0], dtype=np.uint8)
    return rgb


def save_montage(overlays: list[np.ndarray], out_path: Path, n_cols: int = 6) -> None:
    if not overlays:
        return
    h, w = overlays[0].shape[:2]
    thumb_w = 360
    thumb_h = max(1, int(round(h * thumb_w / w)))
    thumbs = [resize_nearest(img, thumb_h, thumb_w) for img in overlays]
    n_rows = math.ceil(len(thumbs) / n_cols)
    canvas = np.zeros((n_rows * thumb_h, n_cols * thumb_w, 3), dtype=np.uint8)
    for i, thumb in enumerate(thumbs):
        row = i // n_cols
        col = i % n_cols
        canvas[row * thumb_h:(row + 1) * thumb_h, col * thumb_w:(col + 1) * thumb_w, :] = thumb
    iio.imwrite(out_path, canvas)


def resize_nearest(img: np.ndarray, out_h: int, out_w: int) -> np.ndarray:
    y_idx = np.linspace(0, img.shape[0] - 1, out_h).astype(int)
    x_idx = np.linspace(0, img.shape[1] - 1, out_w).astype(int)
    return img[y_idx][:, x_idx]


def write_manifest(samples: list[FrameSample], rows: list[dict[str, object]], out_path: Path) -> None:
    fieldnames = [
        "sample_id",
        "well",
        "position",
        "frame",
        "frame_index0",
        "image_path",
        "mask_tif",
        "overlay_png",
        "n_objects",
        "mask_max",
    ]
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for sample, row in zip(samples, rows):
            writer.writerow(
                {
                    "sample_id": row["sample_id"],
                    "well": sample.well,
                    "position": sample.position,
                    "frame": sample.frame_label,
                    "frame_index0": sample.frame_index0,
                    "image_path": str(sample.image_path),
                    "mask_tif": row["mask_tif"],
                    "overlay_png": row["overlay_png"],
                    "n_objects": row["n_objects"],
                    "mask_max": row["mask_max"],
                }
            )


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    masks_dir = args.out_dir / "masks"
    overlays_dir = args.out_dir / "overlays"
    masks_dir.mkdir(exist_ok=True)
    overlays_dir.mkdir(exist_ok=True)

    sites = list_sites(args.images_dir)
    samples = choose_samples(sites, args.n_frames, args.max_sites, args.seed, args.frame_base)
    frames = [read_frame(sample) for sample in samples]
    model_inputs = [prepare_model_input(frame, args.channel_mode) for frame in frames]

    print(f"Images dir: {args.images_dir}")
    print(f"Sites available: {len(sites)}")
    print(f"Frames sampled: {len(samples)}")
    print(f"Output dir: {args.out_dir}")
    print(f"Running CPSAM with gpu={not args.cpu}")

    model = models.CellposeModel(gpu=not args.cpu, pretrained_model="cpsam")
    result = model.eval(
        model_inputs,
        diameter=args.diameter,
        flow_threshold=args.flow_threshold,
        cellprob_threshold=args.cellprob_threshold,
        batch_size=args.batch_size,
        channel_axis=-1 if args.channel_mode == "rgb" else None,
    )
    masks = result[0] if isinstance(result, tuple) else result
    if isinstance(masks, np.ndarray) and masks.ndim == 3:
        mask_list = [masks[i] for i in range(masks.shape[0])]
    else:
        mask_list = list(masks)

    overlays: list[np.ndarray] = []
    manifest_rows: list[dict[str, object]] = []
    for i, (sample, frame, mask) in enumerate(zip(samples, frames, mask_list), start=1):
        sample_id = f"{i:02d}_{sample.well}_{sample.position}_frame{sample.frame_label:03d}"
        mask_path = masks_dir / f"{sample_id}_cpsam_mask.tiff"
        overlay_path = overlays_dir / f"{sample_id}_cpsam_outline.png"
        tf.imwrite(mask_path, np.asarray(mask).astype(np.uint32), photometric="minisblack")
        overlay = outline_overlay(frame, np.asarray(mask))
        iio.imwrite(overlay_path, overlay)
        overlays.append(overlay)
        positive = np.asarray(mask)[np.asarray(mask) > 0]
        manifest_rows.append(
            {
                "sample_id": sample_id,
                "mask_tif": str(mask_path),
                "overlay_png": str(overlay_path),
                "n_objects": int(len(np.unique(positive))) if positive.size else 0,
                "mask_max": int(np.max(mask)) if np.size(mask) else 0,
            }
        )
        print(f"Wrote {overlay_path}")

    write_manifest(samples, manifest_rows, args.out_dir / "sample_manifest.csv")
    if not args.no_montage:
        save_montage(overlays, args.out_dir / "cpsam_outline_montage.png")
    print(f"Wrote manifest: {args.out_dir / 'sample_manifest.csv'}")


if __name__ == "__main__":
    main()
