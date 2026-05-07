from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def percentile_normalize(x: np.ndarray, lo: float = 1, hi: float = 99.8, eps: float = 1e-6) -> np.ndarray:
    x = x.astype(np.float32, copy=False)
    a, b = np.percentile(x, (lo, hi))
    return np.clip((x - a) / (b - a + eps), 0, 1)


def centered_xy_crop(img: np.ndarray, crop_size: int) -> np.ndarray:
    if crop_size <= 0:
        return img
    if img.ndim < 2:
        raise ValueError(f"Expected at least 2 dimensions for XY crop, got {img.shape}")
    y, x = img.shape[-2:]
    if crop_size > y or crop_size > x:
        raise ValueError(f"Crop size {crop_size} is larger than image XY shape {(y, x)}")
    cy, cx = y // 2, x // 2
    half = crop_size // 2
    return img[..., cy - half : cy - half + crop_size, cx - half : cx - half + crop_size]


def build_bf_fucci_fused_input(
    zcyx: np.ndarray,
    bf_channel: int = 0,
    red_channel: int = 1,
    green_channel: int = 2,
    invert_bf: bool = False,
    smooth_sigma: tuple[float, float, float] | None = (0.5, 0.8, 0.8),
) -> np.ndarray:
    if zcyx.ndim != 4:
        raise ValueError(f"Expected ZCYX input, got {zcyx.shape}")

    bf = percentile_normalize(zcyx[:, bf_channel, :, :])
    if invert_bf:
        bf = 1 - bf

    red = percentile_normalize(zcyx[:, red_channel, :, :])
    green = percentile_normalize(zcyx[:, green_channel, :, :])
    fucci_fused = np.maximum(red, green)
    if smooth_sigma is not None:
        fucci_fused = ndi.gaussian_filter(fucci_fused, sigma=smooth_sigma)

    return np.stack([bf, fucci_fused], axis=1).astype(np.float32)


def build_fucci_bf_absz_product_input(
    zcyx: np.ndarray,
    bf_channel: int = 1,
    red_channel: int = 0,
    green_channel: int = 2,
    blur_sigma_xy: float = 20,
    z_anisotropy: float = 1.44,
) -> np.ndarray:
    if zcyx.ndim != 4:
        raise ValueError(f"Expected ZCYX input, got {zcyx.shape}")

    from image_datamining.fucci_focus_channel.focus import build_abs_zscore_smooth

    red = percentile_normalize(zcyx[:, red_channel, :, :])
    green = percentile_normalize(zcyx[:, green_channel, :, :])
    fucci_merge = np.maximum(red, green)

    bf_raw = zcyx[:, bf_channel, :, :]
    bf = percentile_normalize(bf_raw)

    focus = build_abs_zscore_smooth(
        bf_raw,
        blur_sigma_xy=blur_sigma_xy,
        z_anisotropy=z_anisotropy,
    )
    focus_product = percentile_normalize(focus["abs_zscore"] * focus["smooth"])

    return np.stack([fucci_merge, bf, focus_product], axis=1).astype(np.float32)
