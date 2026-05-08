from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi


@dataclass(frozen=True)
class NearestMaskGraph:
    labels: np.ndarray
    nearest_other_dist: np.ndarray
    nearest_other_label: np.ndarray
    source_y: np.ndarray
    source_x: np.ndarray
    neighbor_y: np.ndarray
    neighbor_x: np.ndarray


def nearest_mask_pixels(label_image: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return per-pixel distance and label ID for the nearest nonzero mask."""
    labels = _validate_label_image(label_image)
    occupied = labels > 0
    dist, inds = ndi.distance_transform_edt(~occupied, return_indices=True)
    nearest_label = labels[tuple(inds)]
    dist[occupied] = 0.0
    nearest_label[occupied] = labels[occupied]
    return dist, nearest_label


def nearest_mask_voronoi(label_image: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return distance, nearest label, and nearest mask-pixel coordinates per pixel."""
    labels = _validate_label_image(label_image)
    occupied = labels > 0
    dist, inds = ndi.distance_transform_edt(~occupied, return_indices=True)
    nearest_label = labels[tuple(inds)]
    dist[occupied] = 0.0
    nearest_label[occupied] = labels[occupied]
    return dist, nearest_label, inds


def nearest_mask_neighbors(
    label_image: np.ndarray,
    *,
    connectivity: int = 8,
) -> NearestMaskGraph:
    """Return exact nearest other mask distances from EDT/Voronoi adjacencies."""
    labels = _validate_label_image(label_image)
    max_label = int(labels.max(initial=0))
    nearest_other_dist = np.full(max_label + 1, np.inf, dtype=np.float64)
    nearest_other_label = np.full(max_label + 1, -1, dtype=np.int64)
    source_y = np.full(max_label + 1, -1, dtype=np.int64)
    source_x = np.full(max_label + 1, -1, dtype=np.int64)
    neighbor_y = np.full(max_label + 1, -1, dtype=np.int64)
    neighbor_x = np.full(max_label + 1, -1, dtype=np.int64)
    object_labels = np.unique(labels[labels > 0]).astype(np.int64, copy=False)
    if object_labels.size <= 1:
        return NearestMaskGraph(
            object_labels,
            nearest_other_dist,
            nearest_other_label,
            source_y,
            source_x,
            neighbor_y,
            neighbor_x,
        )

    _, nearest_label, inds = nearest_mask_voronoi(labels)
    offsets = [(1, 0), (0, 1)]
    if connectivity == 8:
        offsets.extend([(1, 1), (1, -1)])
    elif connectivity != 4:
        raise ValueError("connectivity must be 4 or 8")

    for dy, dx in offsets:
        y0, y1 = _offset_slices(dy, labels.shape[0])
        x0, x1 = _offset_slices(dx, labels.shape[1])
        lab_a = nearest_label[y0, x0]
        lab_b = nearest_label[y1, x1]
        different = (lab_a > 0) & (lab_b > 0) & (lab_a != lab_b)
        if not np.any(different):
            continue

        src_a_y = inds[0, y0, x0][different].astype(np.int64, copy=False)
        src_a_x = inds[1, y0, x0][different].astype(np.int64, copy=False)
        src_b_y = inds[0, y1, x1][different].astype(np.int64, copy=False)
        src_b_x = inds[1, y1, x1][different].astype(np.int64, copy=False)
        labels_a = lab_a[different].astype(np.int64, copy=False)
        labels_b = lab_b[different].astype(np.int64, copy=False)
        dists = np.hypot(src_a_y - src_b_y, src_a_x - src_b_x)

        _update_best(
            labels_a,
            labels_b,
            dists,
            src_a_y,
            src_a_x,
            src_b_y,
            src_b_x,
            nearest_other_dist,
            nearest_other_label,
            source_y,
            source_x,
            neighbor_y,
            neighbor_x,
        )
        _update_best(
            labels_b,
            labels_a,
            dists,
            src_b_y,
            src_b_x,
            src_a_y,
            src_a_x,
            nearest_other_dist,
            nearest_other_label,
            source_y,
            source_x,
            neighbor_y,
            neighbor_x,
        )

    return NearestMaskGraph(
        object_labels,
        nearest_other_dist,
        nearest_other_label,
        source_y,
        source_x,
        neighbor_y,
        neighbor_x,
    )


def _validate_label_image(label_image: np.ndarray) -> np.ndarray:
    labels = np.asarray(label_image)
    if labels.ndim != 2:
        raise ValueError(f"Expected a 2D label image, got shape {labels.shape}")
    if np.any(labels < 0):
        raise ValueError("Label image must use nonnegative labels")
    return labels


def _offset_slices(delta: int, size: int) -> tuple[slice, slice]:
    if delta > 0:
        return slice(None, -delta), slice(delta, None)
    if delta < 0:
        return slice(-delta, None), slice(None, delta)
    return slice(None), slice(None)


def _update_best(
    source_labels: np.ndarray,
    other_labels: np.ndarray,
    dists: np.ndarray,
    y0: np.ndarray,
    x0: np.ndarray,
    y1: np.ndarray,
    x1: np.ndarray,
    nearest_other_dist: np.ndarray,
    nearest_other_label: np.ndarray,
    source_y: np.ndarray,
    source_x: np.ndarray,
    neighbor_y: np.ndarray,
    neighbor_x: np.ndarray,
) -> None:
    order = np.lexsort((other_labels, dists, source_labels))
    sorted_source = source_labels[order]
    first_for_source = np.r_[True, sorted_source[1:] != sorted_source[:-1]]
    best_idx = order[first_for_source]
    best_source = source_labels[best_idx]
    improved = dists[best_idx] < nearest_other_dist[best_source]
    best_idx = best_idx[improved]
    best_source = best_source[improved]
    nearest_other_dist[best_source] = dists[best_idx]
    nearest_other_label[best_source] = other_labels[best_idx]
    source_y[best_source] = y0[best_idx]
    source_x[best_source] = x0[best_idx]
    neighbor_y[best_source] = y1[best_idx]
    neighbor_x[best_source] = x1[best_idx]
