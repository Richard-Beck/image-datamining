#!/usr/bin/env python3
"""Run CellposeSAM on one K00 Images_40Frames TIFF stack."""

from __future__ import annotations

import argparse
import csv
import re
import time
from pathlib import Path

import numpy as np
import tifffile as tf
from cellpose import models


DEFAULT_OUT_DIR = Path("analyses/K00_GemcitabineExposure_033023/cpsam_full_stacks")
SITE_RE = re.compile(r"^(?P<well>[A-H]\d{1,2})_(?P<position>\d+)\.tiff?$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run CPSAM on every frame in one K00 registered Images_40Frames "
            "TIFF stack and save one compressed label-mask TIFF."
        )
    )
    parser.add_argument("--input-tiff", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
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
    parser.add_argument("--min-size", type=int, default=15)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--cpu",
        action="store_true",
        help="Run on CPU. GPU is the default and should be used for production.",
    )
    return parser.parse_args()


def site_metadata(path: Path) -> tuple[str, str, str]:
    match = SITE_RE.match(path.name)
    if match:
        well = match.group("well").upper()
        position = match.group("position")
        return well, position, f"{well}_{position}"
    return path.stem, "NA", path.stem


def read_stack(path: Path) -> np.ndarray:
    with tf.TiffFile(path) as tif:
        arr = tif.asarray()
        axes = tif.series[0].axes if tif.series else ""

    if arr.ndim == 3:
        if axes == "YXS" or arr.shape[-1] in (3, 4):
            arr = arr[None, :, :, :]
        else:
            arr = arr[:, :, :, None]
    if arr.ndim != 4:
        raise ValueError(f"Expected a frame stack with shape TYXC/IYXS, got {arr.shape} from {path}")
    if arr.shape[-1] > 3:
        arr = arr[:, :, :, :3]
    if arr.shape[-1] == 1:
        arr = np.repeat(arr, 3, axis=-1)
    return arr


def prepare_model_inputs(stack: np.ndarray, channel_mode: str) -> list[np.ndarray]:
    if channel_mode == "rgb":
        return [stack[i] for i in range(stack.shape[0])]
    channel_index = {"red": 0, "green": 1, "blue": 2}.get(channel_mode)
    if channel_index is None:
        gray = np.mean(stack[:, :, :, :3].astype(np.float32), axis=3)
    else:
        gray = stack[:, :, :, channel_index].astype(np.float32)
    return [gray[i] for i in range(gray.shape[0])]


def coerce_mask_stack(masks: object) -> np.ndarray:
    if isinstance(masks, np.ndarray):
        if masks.ndim == 2:
            return masks[None, :, :]
        if masks.ndim == 3:
            return masks
    mask_list = [np.asarray(mask) for mask in masks]  # type: ignore[arg-type]
    if not mask_list:
        raise ValueError("CPSAM returned no masks")
    return np.stack(mask_list, axis=0)


def frame_counts(mask_stack: np.ndarray) -> list[int]:
    counts: list[int] = []
    for frame_mask in mask_stack:
        labels = np.unique(frame_mask[frame_mask > 0])
        counts.append(int(labels.size))
    return counts


def write_mask_tiff(path: Path, masks: np.ndarray) -> str:
    max_label = int(np.max(masks)) if masks.size else 0
    dtype = np.uint16 if max_label <= np.iinfo(np.uint16).max else np.uint32
    mask_out = masks.astype(dtype, copy=False)
    tf.imwrite(
        path,
        mask_out,
        photometric="minisblack",
        compression="zlib",
        metadata={"axes": "TYX"},
    )
    return np.dtype(dtype).name


def write_manifest_row(path: Path, row: dict[str, object]) -> None:
    fieldnames = [
        "site_id",
        "well",
        "position",
        "input_tiff",
        "mask_tiff",
        "n_frames",
        "height",
        "width",
        "mask_dtype",
        "max_label",
        "total_objects",
        "frame_object_counts",
        "channel_mode",
        "diameter",
        "flow_threshold",
        "cellprob_threshold",
        "min_size",
        "batch_size",
        "gpu",
        "elapsed_seconds",
        "status",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    args = parse_args()
    input_tiff = args.input_tiff.resolve()
    well, position, site_id = site_metadata(input_tiff)

    masks_dir = args.out_dir / "masks"
    manifest_dir = args.out_dir / "manifest_rows"
    masks_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir.mkdir(parents=True, exist_ok=True)

    mask_path = masks_dir / f"{site_id}_cpsam_masks.tiff"
    manifest_path = manifest_dir / f"{site_id}_cpsam_manifest.csv"
    if mask_path.exists() and manifest_path.exists() and not args.overwrite:
        print(f"Existing outputs found for {site_id}; use --overwrite to rerun.")
        return

    started = time.time()
    stack = read_stack(input_tiff)
    model_inputs = prepare_model_inputs(stack, args.channel_mode)

    print(f"Input: {input_tiff}")
    print(f"Stack shape TYXC: {stack.shape}")
    print(f"Output mask: {mask_path}")
    print(f"Running CPSAM with gpu={not args.cpu}")

    model = models.CellposeModel(gpu=not args.cpu, pretrained_model="cpsam")
    result = model.eval(
        model_inputs,
        diameter=args.diameter,
        flow_threshold=args.flow_threshold,
        cellprob_threshold=args.cellprob_threshold,
        batch_size=args.batch_size,
        min_size=args.min_size,
        channel_axis=-1 if args.channel_mode == "rgb" else None,
    )
    masks = coerce_mask_stack(result[0] if isinstance(result, tuple) else result)
    counts = frame_counts(masks)
    mask_dtype = write_mask_tiff(mask_path, masks)
    elapsed = time.time() - started

    row = {
        "site_id": site_id,
        "well": well,
        "position": position,
        "input_tiff": str(input_tiff),
        "mask_tiff": str(mask_path),
        "n_frames": int(masks.shape[0]),
        "height": int(masks.shape[1]),
        "width": int(masks.shape[2]),
        "mask_dtype": mask_dtype,
        "max_label": int(np.max(masks)) if masks.size else 0,
        "total_objects": int(sum(counts)),
        "frame_object_counts": ";".join(str(x) for x in counts),
        "channel_mode": args.channel_mode,
        "diameter": "" if args.diameter is None else args.diameter,
        "flow_threshold": args.flow_threshold,
        "cellprob_threshold": args.cellprob_threshold,
        "min_size": args.min_size,
        "batch_size": args.batch_size,
        "gpu": not args.cpu,
        "elapsed_seconds": round(elapsed, 3),
        "status": "ok",
    }
    write_manifest_row(manifest_path, row)
    print(f"Wrote {mask_path}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
