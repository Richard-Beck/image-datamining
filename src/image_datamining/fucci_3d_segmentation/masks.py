from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import tifffile as tf


def summarize_mask(mask_path: Path, sample_id: str = "") -> list[dict[str, object]]:
    masks = tf.imread(mask_path)
    labels, counts = np.unique(masks[masks > 0], return_counts=True)
    rows: list[dict[str, object]] = []
    for label, count in zip(labels, counts):
        coords = np.argwhere(masks == label)
        mins = coords.min(axis=0)
        maxs = coords.max(axis=0)
        row: dict[str, object] = {
            "sample_id": sample_id,
            "mask_path": str(mask_path),
            "label": int(label),
            "voxel_count": int(count),
        }
        if masks.ndim == 3:
            row.update(
                {
                    "z_min": int(mins[0]),
                    "z_max": int(maxs[0]),
                    "y_min": int(mins[1]),
                    "y_max": int(maxs[1]),
                    "x_min": int(mins[2]),
                    "x_max": int(maxs[2]),
                }
            )
        elif masks.ndim == 2:
            row.update(
                {
                    "z_min": "",
                    "z_max": "",
                    "y_min": int(mins[0]),
                    "y_max": int(maxs[0]),
                    "x_min": int(mins[1]),
                    "x_max": int(maxs[1]),
                }
            )
        rows.append(row)
    return rows


def write_object_summary(rows: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id",
        "mask_path",
        "label",
        "voxel_count",
        "z_min",
        "z_max",
        "y_min",
        "y_max",
        "x_min",
        "x_max",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

