from __future__ import annotations

from pathlib import Path

import numpy as np
import tifffile as tf


def read_channel_file_list(paths_text: str) -> list[Path]:
    return [Path(item) for item in paths_text.split(";") if item]


def read_zcyx_from_manifest_row(row: dict[str, str]) -> np.ndarray:
    channels = []
    for channel in ("ch00_files", "ch01_files", "ch02_files"):
        z_planes = [tf.imread(path) for path in read_channel_file_list(row[channel])]
        channels.append(np.stack(z_planes, axis=0))
    return np.stack(channels, axis=1)


def write_zcyx_tiff(path: Path, img: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tf.imwrite(path, img, imagej=True, metadata={"axes": "ZCYX"})


def write_mask_tiff(path: Path, masks: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tf.imwrite(path, masks.astype(np.uint32))

