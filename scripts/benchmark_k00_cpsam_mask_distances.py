#!/usr/bin/env python3
"""Benchmark nearest-mask distance utilities on sampled K00 CPSAM masks."""

from __future__ import annotations

import argparse
import csv
import random
import sys
import time
from pathlib import Path

import numpy as np
import tifffile as tf
from PIL import Image, ImageDraw

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from image_datamining.mask_distances import nearest_mask_neighbors, nearest_mask_pixels


DEFAULT_MASK_DIR = Path("analyses/K00_GemcitabineExposure_033023/cpsam_frame_sample/masks")
DEFAULT_OUT = Path("analyses/K00_GemcitabineExposure_033023/benchmarks/cpsam_mask_distance_benchmark.csv")
DEFAULT_SAMPLE_MANIFEST = Path("analyses/K00_GemcitabineExposure_033023/cpsam_frame_sample/sample_manifest.csv")
DEFAULT_LINK_OVERLAY_DIR = Path("analyses/K00_GemcitabineExposure_033023/cpsam_nucleus_tracking_links/overlays")
DEFAULT_VALIDATION_PNG = Path(
    "analyses/K00_GemcitabineExposure_033023/benchmarks/cpsam_mask_distance_max_pair_validation.png"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark EDT nearest-pixel and Voronoi-adjacency nearest-object calculations on 2D CPSAM masks."
    )
    parser.add_argument("--mask-dir", type=Path, default=DEFAULT_MASK_DIR)
    parser.add_argument("--glob", default="*_cpsam_mask.tif*")
    parser.add_argument("-n", "--n", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20260508)
    parser.add_argument("--connectivity", type=int, choices=(4, 8), default=8)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--sample-manifest", type=Path, default=DEFAULT_SAMPLE_MANIFEST)
    parser.add_argument("--link-overlay-dir", type=Path, default=DEFAULT_LINK_OVERLAY_DIR)
    parser.add_argument("--validation-png", type=Path, default=DEFAULT_VALIDATION_PNG)
    return parser.parse_args()


def read_2d_mask(path: Path) -> np.ndarray:
    mask = tf.imread(path)
    if mask.ndim != 2:
        raise ValueError(f"Expected a 2D mask TIFF, got shape {mask.shape} from {path}")
    return np.asarray(mask)


def benchmark_mask(path: Path, connectivity: int) -> dict[str, object]:
    loaded_at = time.perf_counter()
    mask = read_2d_mask(path)
    load_seconds = time.perf_counter() - loaded_at

    labels = np.unique(mask[mask > 0])
    edt_started = time.perf_counter()
    dist, nearest_label = nearest_mask_pixels(mask)
    edt_seconds = time.perf_counter() - edt_started

    graph_started = time.perf_counter()
    graph = nearest_mask_neighbors(mask, connectivity=connectivity)
    graph_seconds = time.perf_counter() - graph_started

    finite_neighbor = np.isfinite(graph.nearest_other_dist[graph.labels])
    finite_dists = graph.nearest_other_dist[graph.labels][finite_neighbor]
    max_label = -1
    if finite_dists.size:
        finite_labels = graph.labels[finite_neighbor]
        max_label = int(finite_labels[np.argmax(graph.nearest_other_dist[finite_labels])])
    return {
        "mask_tiff": str(path),
        "height": int(mask.shape[0]),
        "width": int(mask.shape[1]),
        "pixels": int(mask.size),
        "max_label": int(mask.max(initial=0)),
        "n_labels_present": int(labels.size),
        "load_seconds": load_seconds,
        "edt_seconds": edt_seconds,
        "graph_seconds": graph_seconds,
        "total_compute_seconds": edt_seconds + graph_seconds,
        "dist_dtype": str(dist.dtype),
        "nearest_label_dtype": str(nearest_label.dtype),
        "n_labels_with_neighbor": int(finite_neighbor.sum()),
        "n_labels_without_neighbor": int((~finite_neighbor).sum()),
        "min_nearest_other_px": float(np.min(finite_dists)) if finite_dists.size else "",
        "median_nearest_other_px": float(np.median(finite_dists)) if finite_dists.size else "",
        "p90_nearest_other_px": float(np.quantile(finite_dists, 0.90)) if finite_dists.size else "",
        "p99_nearest_other_px": float(np.quantile(finite_dists, 0.99)) if finite_dists.size else "",
        "max_nearest_other_px": float(np.max(finite_dists)) if finite_dists.size else "",
        "max_pair_label_a": max_label if max_label >= 0 else "",
        "max_pair_label_b": int(graph.nearest_other_label[max_label]) if max_label >= 0 else "",
        "max_pair_y0": int(graph.source_y[max_label]) if max_label >= 0 else "",
        "max_pair_x0": int(graph.source_x[max_label]) if max_label >= 0 else "",
        "max_pair_y1": int(graph.neighbor_y[max_label]) if max_label >= 0 else "",
        "max_pair_x1": int(graph.neighbor_x[max_label]) if max_label >= 0 else "",
    }


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "mask_tiff",
        "height",
        "width",
        "pixels",
        "max_label",
        "n_labels_present",
        "load_seconds",
        "edt_seconds",
        "graph_seconds",
        "total_compute_seconds",
        "dist_dtype",
        "nearest_label_dtype",
        "n_labels_with_neighbor",
        "n_labels_without_neighbor",
        "min_nearest_other_px",
        "median_nearest_other_px",
        "p90_nearest_other_px",
        "p99_nearest_other_px",
        "max_nearest_other_px",
        "max_pair_label_a",
        "max_pair_label_b",
        "max_pair_y0",
        "max_pair_x0",
        "max_pair_y1",
        "max_pair_x1",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict[str, object]], output: Path) -> None:
    totals = np.asarray([float(row["total_compute_seconds"]) for row in rows], dtype=float)
    edts = np.asarray([float(row["edt_seconds"]) for row in rows], dtype=float)
    graphs = np.asarray([float(row["graph_seconds"]) for row in rows], dtype=float)
    max_dists = np.asarray([float(row["max_nearest_other_px"]) for row in rows if row["max_nearest_other_px"]], dtype=float)
    labels = np.asarray([int(row["n_labels_present"]) for row in rows], dtype=int)
    unresolved = np.asarray([int(row["n_labels_without_neighbor"]) for row in rows], dtype=int)
    pixels = np.asarray([int(row["pixels"]) for row in rows], dtype=int)
    print(f"Benchmarked {len(rows)} masks")
    print(f"Pixels per mask: min={pixels.min()} median={int(np.median(pixels))} max={pixels.max()}")
    print(f"Labels per mask: min={labels.min()} median={int(np.median(labels))} max={labels.max()}")
    print(f"EDT seconds: median={np.median(edts):.4f} max={edts.max():.4f}")
    print(f"Voronoi graph seconds: median={np.median(graphs):.4f} max={graphs.max():.4f}")
    print(f"Total compute seconds: median={np.median(totals):.4f} max={totals.max():.4f}")
    print(f"Labels without neighbor: max={unresolved.max()}")
    if max_dists.size:
        print(f"Per-mask max nearest-other px: median={np.median(max_dists):.3f} max={max_dists.max():.3f}")
    print(f"Wrote CSV: {output}")


def read_sample_manifest(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {Path(row["mask_tif"]).name: row for row in rows}


def validation_base_image(row: dict[str, object], manifest: dict[str, dict[str, str]], link_overlay_dir: Path) -> np.ndarray:
    mask_name = Path(str(row["mask_tiff"])).name
    sample = manifest.get(mask_name)
    if sample:
        site_id = f"{sample['well']}_{sample['position']}"
        link_overlay = link_overlay_dir / f"{site_id}_cpsam_nucleus_track_overlay_framesall_excluding_frame0.tiff"
        frame = int(sample["frame"])
        if link_overlay.exists() and frame > 0:
            with tf.TiffFile(link_overlay) as tif:
                return np.asarray(tif.pages[frame - 1].asarray())
        overlay_png = Path(sample["overlay_png"])
        if overlay_png.exists():
            return np.asarray(Image.open(overlay_png).convert("RGB"))
    return np.asarray(read_2d_mask(Path(str(row["mask_tiff"]))))


def render_validation_png(rows: list[dict[str, object]], manifest_path: Path, link_overlay_dir: Path, output: Path) -> None:
    candidates = [row for row in rows if row["max_nearest_other_px"]]
    if not candidates:
        return
    row = max(candidates, key=lambda item: float(item["max_nearest_other_px"]))
    manifest = read_sample_manifest(manifest_path)
    base = validation_base_image(row, manifest, link_overlay_dir)
    if base.ndim == 2:
        lo, hi = np.percentile(base, [1, 99])
        scaled = np.clip((base.astype(np.float32) - lo) / max(hi - lo, 1), 0, 1)
        rgb = np.repeat((scaled * 255).astype(np.uint8)[:, :, None], 3, axis=2)
    elif base.ndim == 3 and base.shape[-1] >= 3:
        rgb = base[:, :, :3].astype(np.uint8, copy=False)
    else:
        raise ValueError(f"Cannot render validation base image with shape {base.shape}")

    image = Image.fromarray(rgb).convert("RGB")
    draw = ImageDraw.Draw(image)
    x0 = int(row["max_pair_x0"])
    y0 = int(row["max_pair_y0"])
    x1 = int(row["max_pair_x1"])
    y1 = int(row["max_pair_y1"])
    draw.line((x0, y0, x1, y1), fill=(255, 255, 0), width=4)
    r = 7
    draw.ellipse((x0 - r, y0 - r, x0 + r, y0 + r), outline=(255, 0, 0), width=3)
    draw.ellipse((x1 - r, y1 - r, x1 + r, y1 + r), outline=(0, 255, 255), width=3)
    text = (
        f"{Path(str(row['mask_tiff'])).name} "
        f"{row['max_pair_label_a']}->{row['max_pair_label_b']} "
        f"d={float(row['max_nearest_other_px']):.2f}px"
    )
    draw.rectangle((8, 8, 8 + 8 * len(text), 32), fill=(0, 0, 0))
    draw.text((12, 12), text, fill=(255, 255, 255))
    output.parent.mkdir(parents=True, exist_ok=True)
    image.save(output)
    print(f"Wrote validation PNG: {output}")


def main() -> None:
    args = parse_args()
    paths = sorted(args.mask_dir.glob(args.glob))
    if not paths:
        raise FileNotFoundError(f"No masks matched {args.mask_dir / args.glob}")
    rng = random.Random(args.seed)
    selected = rng.sample(paths, k=min(args.n, len(paths)))

    rows = []
    for path in selected:
        row = benchmark_mask(path, args.connectivity)
        rows.append(row)
        print(
            f"{Path(str(row['mask_tiff'])).name}: "
            f"labels={row['n_labels_present']} "
            f"edt={float(row['edt_seconds']):.4f}s "
            f"graph={float(row['graph_seconds']):.4f}s",
            flush=True,
        )
    write_rows(args.out, rows)
    render_validation_png(rows, args.sample_manifest, args.link_overlay_dir, args.validation_png)
    print_summary(rows, args.out)


if __name__ == "__main__":
    main()
