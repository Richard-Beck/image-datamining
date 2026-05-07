#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from _bootstrap import add_src_to_path

add_src_to_path()

import tifffile as tf

from image_datamining.fucci_3d_segmentation.segment import run_cellpose_cpsam


def main() -> None:
    parser = argparse.ArgumentParser(description="Run CellposeSAM on a prepared FUCCI input.")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mode", choices=["3d", "2d-per-z"], required=True)
    parser.add_argument("--anisotropy", type=float, default=1.44)
    parser.add_argument("--min-size", type=int, default=2000)
    parser.add_argument("--gpu", action="store_true")
    args = parser.parse_args()

    img = tf.imread(args.input)
    if img.ndim != 4:
        raise ValueError(f"Expected prepared ZCYX input, got {img.shape}")

    masks = run_cellpose_cpsam(
        img,
        mode=args.mode,
        anisotropy=args.anisotropy,
        min_size=args.min_size,
        gpu=args.gpu,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    tf.imwrite(args.output, masks.astype("uint32"))
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
