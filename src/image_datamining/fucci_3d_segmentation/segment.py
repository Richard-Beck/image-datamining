from __future__ import annotations

from typing import Literal

import numpy as np


def run_cellpose_cpsam(
    img: np.ndarray,
    mode: Literal["3d", "2d-per-z"],
    anisotropy: float,
    min_size: int,
    gpu: bool,
) -> np.ndarray:
    from cellpose import models

    model = models.CellposeModel(gpu=gpu, pretrained_model="cpsam")
    if mode == "3d":
        result = model.eval(
            img,
            do_3D=True,
            z_axis=0,
            channel_axis=1,
            anisotropy=anisotropy,
            min_size=min_size,
        )
        masks = result[0] if isinstance(result, tuple) else result
        return masks.astype(np.uint32, copy=False)

    masks = []
    label_offset = 0
    for z in range(img.shape[0]):
        plane = img[z]
        result = model.eval(plane, channel_axis=0, min_size=min_size)
        mask = result[0] if isinstance(result, tuple) else result
        mask = mask.astype(np.uint32, copy=False)
        positive = mask > 0
        if positive.any():
            mask[positive] += label_offset
            label_offset = int(mask.max())
        masks.append(mask)
    return np.stack(masks, axis=0)

