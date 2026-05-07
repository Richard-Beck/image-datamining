from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def robust01(x: np.ndarray, lo: float = 1, hi: float = 99.5, eps: float = 1e-6) -> np.ndarray:
    x = x.astype(np.float32, copy=False)
    a, b = np.percentile(x, (lo, hi))
    return np.clip((x - a) / (b - a + eps), 0, 1)


def zscore_volume(x: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    return (x - np.mean(x)) / (np.std(x) + eps)


def normalize_brightfield(
    bf_zyx: np.ndarray,
    background_sigma: float = 24,
    highpass_sigma: float = 8,
    eps: float = 1e-6,
) -> tuple[np.ndarray, np.ndarray]:
    bf = bf_zyx.astype(np.float32, copy=False)
    background = ndi.gaussian_filter(bf, sigma=(0, background_sigma, background_sigma))
    centered = bf - background
    local_scale = ndi.gaussian_filter(np.abs(centered), sigma=(0, background_sigma, background_sigma))
    bf_norm = centered / (local_scale + eps)
    bf_hp = bf_norm - ndi.gaussian_filter(bf_norm, sigma=(0, highpass_sigma, highpass_sigma))
    return bf_norm.astype(np.float32), bf_hp.astype(np.float32)


def local_std_xy(x: np.ndarray, window: int = 9, eps: float = 1e-6) -> np.ndarray:
    size = (1, window, window)
    mean = ndi.uniform_filter(x, size=size)
    mean_sq = ndi.uniform_filter(x * x, size=size)
    var = np.maximum(mean_sq - mean * mean, 0)
    return np.sqrt(var + eps).astype(np.float32)


def tenengrad_xy(x: np.ndarray) -> np.ndarray:
    sx = ndi.sobel(x, axis=2)
    sy = ndi.sobel(x, axis=1)
    return (sx * sx + sy * sy).astype(np.float32)


def abs_log_xy(x: np.ndarray, sigma: float = 1.2) -> np.ndarray:
    return np.abs(ndi.gaussian_laplace(x, sigma=(0, sigma, sigma))).astype(np.float32)


def focus_score_stack(
    bf_zyx: np.ndarray,
    background_sigma: float = 24,
    highpass_sigma: float = 8,
    local_std_window: int = 9,
    log_sigma: float = 1.2,
    smooth_sigma: tuple[float, float, float] = (0.5, 1.0, 1.0),
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    bf_norm, bf_hp = normalize_brightfield(
        bf_zyx,
        background_sigma=background_sigma,
        highpass_sigma=highpass_sigma,
    )
    local_sd = local_std_xy(bf_hp, window=local_std_window)
    ten = tenengrad_xy(bf_hp)
    log_abs = abs_log_xy(bf_hp, sigma=log_sigma)
    score = zscore_volume(local_sd) + zscore_volume(ten) + 0.5 * zscore_volume(log_abs)
    score = ndi.gaussian_filter(score, sigma=smooth_sigma).astype(np.float32)
    return score, bf_norm, bf_hp


def focus_component_stack(
    bf_zyx: np.ndarray,
    method: str,
    background_sigma: float = 24,
    highpass_sigma: float = 8,
    local_std_window: int = 9,
    log_sigma: float = 1.2,
    smooth_sigma: tuple[float, float, float] = (0.5, 1.0, 1.0),
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    bf_norm, bf_hp = normalize_brightfield(
        bf_zyx,
        background_sigma=background_sigma,
        highpass_sigma=highpass_sigma,
    )
    if method == "local_sd":
        score = zscore_volume(local_std_xy(bf_hp, window=local_std_window))
    elif method == "tenengrad":
        score = zscore_volume(tenengrad_xy(bf_hp))
    elif method == "log":
        score = zscore_volume(abs_log_xy(bf_hp, sigma=log_sigma))
    elif method == "composite":
        local_sd = local_std_xy(bf_hp, window=local_std_window)
        ten = tenengrad_xy(bf_hp)
        log_abs = abs_log_xy(bf_hp, sigma=log_sigma)
        score = zscore_volume(local_sd) + zscore_volume(ten) + 0.5 * zscore_volume(log_abs)
    else:
        raise ValueError(f"Unknown focus method: {method}")
    score = ndi.gaussian_filter(score, sigma=smooth_sigma).astype(np.float32)
    return score, bf_norm, bf_hp


def relative_focus_confidence(score_zyx: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    shifted = score_zyx - np.min(score_zyx, axis=0, keepdims=True)
    denom = np.max(shifted, axis=0, keepdims=True)
    return (shifted / (denom + eps)).astype(np.float32)


def softmax_focus_confidence(score_zyx: np.ndarray, tau: float = 1.0) -> np.ndarray:
    scaled = score_zyx / tau
    scaled = scaled - np.max(scaled, axis=0, keepdims=True)
    exp = np.exp(scaled)
    return (exp / np.sum(exp, axis=0, keepdims=True)).astype(np.float32)


def build_focus_confidence(
    bf_zyx: np.ndarray,
    confidence_mode: str = "relative",
    tau: float = 1.0,
    **kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    score, bf_norm, bf_hp = focus_score_stack(bf_zyx, **kwargs)
    if confidence_mode == "relative":
        confidence = relative_focus_confidence(score)
    elif confidence_mode == "softmax":
        confidence = softmax_focus_confidence(score, tau=tau)
    else:
        raise ValueError(f"Unknown confidence_mode: {confidence_mode}")
    return confidence, score, bf_norm, bf_hp


def build_component_focus_confidence(
    bf_zyx: np.ndarray,
    method: str,
    confidence_mode: str = "relative",
    tau: float = 1.0,
    **kwargs,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    score, bf_norm, bf_hp = focus_component_stack(bf_zyx, method=method, **kwargs)
    if confidence_mode == "relative":
        confidence = relative_focus_confidence(score)
    elif confidence_mode == "softmax":
        confidence = softmax_focus_confidence(score, tau=tau)
    else:
        raise ValueError(f"Unknown confidence_mode: {confidence_mode}")
    return confidence, score, bf_norm, bf_hp


def build_sd_blur_product(
    bf_zyx: np.ndarray,
    local_std_window: int = 12,
    blur_sigma_xy: float = 20,
    z_anisotropy: float = 1.44,
    background_sigma: float = 24,
    highpass_sigma: float = 8,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    score, _bf_norm, _bf_hp = focus_component_stack(
        bf_zyx,
        method="local_sd",
        background_sigma=background_sigma,
        highpass_sigma=highpass_sigma,
        local_std_window=local_std_window,
    )
    sd_confidence = relative_focus_confidence(score)
    sigma_z = blur_sigma_xy / z_anisotropy
    blur = ndi.gaussian_filter(sd_confidence, sigma=(sigma_z, blur_sigma_xy, blur_sigma_xy)).astype(np.float32)
    product = (sd_confidence * blur).astype(np.float32)
    return blur, sd_confidence, product, score


def otsu_threshold(x: np.ndarray, nbins: int = 256) -> float:
    values = x[np.isfinite(x)].astype(np.float32, copy=False)
    if values.size == 0:
        return 0.0
    vmin = float(values.min())
    vmax = float(values.max())
    if vmax <= vmin:
        return vmin
    counts, edges = np.histogram(values, bins=nbins, range=(vmin, vmax))
    centers = (edges[:-1] + edges[1:]) / 2
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    mean1 = np.cumsum(counts * centers) / np.maximum(weight1, 1)
    mean2 = (np.cumsum((counts * centers)[::-1]) / np.maximum(weight2[::-1], 1))[::-1]
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2
    if variance12.size == 0:
        return vmin
    return float(centers[:-1][np.argmax(variance12)])


def build_sd_log_combo(
    bf_zyx: np.ndarray,
    sd_window: int = 10,
    log_sigma: float = 2.0,
    blur_sigma_xy: float = 20,
    z_anisotropy: float = 1.44,
    background_sigma: float = 24,
    highpass_sigma: float = 8,
) -> dict[str, np.ndarray | float]:
    sd_score, _bf_norm, _bf_hp = focus_component_stack(
        bf_zyx,
        method="local_sd",
        background_sigma=background_sigma,
        highpass_sigma=highpass_sigma,
        local_std_window=sd_window,
        log_sigma=log_sigma,
    )
    log_score, _bf_norm, _bf_hp = focus_component_stack(
        bf_zyx,
        method="log",
        background_sigma=background_sigma,
        highpass_sigma=highpass_sigma,
        local_std_window=sd_window,
        log_sigma=log_sigma,
    )
    sd_confidence = relative_focus_confidence(sd_score)
    log_confidence = relative_focus_confidence(log_score)
    combined = np.sqrt(np.clip(sd_confidence, 0, None) * np.clip(log_confidence, 0, None)).astype(np.float32)
    sigma_z = blur_sigma_xy / z_anisotropy
    blurred = ndi.gaussian_filter(combined, sigma=(sigma_z, blur_sigma_xy, blur_sigma_xy)).astype(np.float32)
    product = (combined * blurred).astype(np.float32)
    threshold = otsu_threshold(blurred)
    otsu_mask = (blurred >= threshold).astype(np.uint8)
    otsu_product = (product * otsu_mask).astype(np.float32)
    return {
        "sd_confidence": sd_confidence,
        "log_confidence": log_confidence,
        "combined": combined,
        "blurred": blurred,
        "product": product,
        "otsu_mask": otsu_mask,
        "otsu_product": otsu_product,
        "sd_score": sd_score,
        "log_score": log_score,
        "otsu_threshold": threshold,
    }


def build_abs_zscore_smooth(
    bf_zyx: np.ndarray,
    blur_sigma_xy: float = 20,
    z_anisotropy: float = 1.44,
    eps: float = 1e-6,
) -> dict[str, np.ndarray]:
    bf = bf_zyx.astype(np.float32, copy=False)
    zscore = (bf - np.mean(bf)) / (np.std(bf) + eps)
    abs_zscore = np.abs(zscore).astype(np.float32)
    sigma_z = blur_sigma_xy / z_anisotropy
    smooth = ndi.gaussian_filter(abs_zscore, sigma=(sigma_z, blur_sigma_xy, blur_sigma_xy)).astype(np.float32)
    return {
        "zscore": zscore.astype(np.float32),
        "abs_zscore": abs_zscore,
        "smooth": smooth,
    }
