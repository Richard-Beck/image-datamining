#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile as tf


DEFAULT_RUN_DIR = Path("pipelines/fucci_3d_segmentation/runs/fucci_absz_sd2_product_cpsam_2d_fullstacks")
DEFAULT_RAW_CROP_DIR = Path(
    "pipelines/fucci_3d_segmentation/runs/jackson_fucci_crop256_bf_fucci_fused_3d/raw_crops"
)
INPUT_KIND = "fucci_absz_sd2_product"
MASK_SUFFIX = f"_{INPUT_KIND}_2d_fullstack_masks.tif"
RAW_CROP_SUFFIX = "_raw3ch_crop256_ZCYX.tif"


class UnionFind:
    def __init__(self, labels: list[int]) -> None:
        self.parent = {label: label for label in labels}

    def find(self, label: int) -> int:
        parent = self.parent[label]
        if parent != label:
            self.parent[label] = self.find(parent)
        return self.parent[label]

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if rb < ra:
            ra, rb = rb, ra
        self.parent[rb] = ra


def parse_int_list(text: str) -> list[int]:
    values = [int(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("Expected at least one comma-separated integer")
    return values


def sample_id_from_mask(path: Path) -> str:
    if path.name.endswith(MASK_SUFFIX):
        return path.name[: -len(MASK_SUFFIX)]
    return re.sub(r"\.tiff?$", "", path.name, flags=re.IGNORECASE)


def discover_masks(run_dir: Path, sample_id: str, limit: int) -> list[Path]:
    mask_dir = run_dir / "masks"
    if sample_id:
        path = mask_dir / f"{sample_id}{MASK_SUFFIX}"
        if not path.exists():
            raise FileNotFoundError(path)
        return [path]
    masks = sorted(mask_dir.glob(f"*{MASK_SUFFIX}"))
    if limit:
        masks = masks[:limit]
    if not masks:
        raise FileNotFoundError(f"No masks found under {mask_dir}")
    return masks


def robust_u8(x: np.ndarray) -> np.ndarray:
    x = x.astype(np.float32, copy=False)
    lo, hi = np.percentile(x, [1, 99])
    if hi <= lo:
        return np.zeros(x.shape, dtype=np.uint8)
    y = np.clip((x - lo) / (hi - lo), 0, 1)
    return (y * 255).astype(np.uint8)


def label_color(label: int) -> tuple[int, int, int]:
    return (
        int((label * 37 + 41) % 255),
        int((label * 67 + 89) % 255),
        int((label * 97 + 131) % 255),
    )


def colorize_labels(mask: np.ndarray) -> np.ndarray:
    out = np.zeros(mask.shape + (3,), dtype=np.uint8)
    for label in np.unique(mask):
        if label == 0:
            continue
        out[mask == label] = label_color(int(label))
    return out


def component_stats(masks: np.ndarray) -> dict[int, dict[str, float]]:
    stats: dict[int, dict[str, float]] = {}
    for label in np.unique(masks):
        if label == 0:
            continue
        zz, yy, xx = np.nonzero(masks == label)
        stats[int(label)] = {
            "area_px": float(zz.size),
            "z_min": float(zz.min()),
            "z_max": float(zz.max()),
            "centroid_y": float(yy.mean()),
            "centroid_x": float(xx.mean()),
        }
    return stats


def find_links(
    masks: np.ndarray,
    max_dz: int,
    min_overlap_fraction: float,
) -> list[dict[str, object]]:
    max_label = int(masks.max())
    if max_label == 0:
        return []
    areas = np.bincount(masks.reshape(-1), minlength=max_label + 1)
    links: list[dict[str, object]] = []

    for dz in range(1, max_dz + 1):
        for z in range(0, masks.shape[0] - dz):
            a = masks[z]
            b = masks[z + dz]
            overlap = (a > 0) & (b > 0)
            if not np.any(overlap):
                continue
            pair_codes, counts = np.unique(a[overlap].astype(np.uint64) * (max_label + 1) + b[overlap], return_counts=True)
            for code, intersection in zip(pair_codes, counts):
                label_a = int(code // (max_label + 1))
                label_b = int(code % (max_label + 1))
                smaller_area = min(int(areas[label_a]), int(areas[label_b]))
                union = int(areas[label_a]) + int(areas[label_b]) - int(intersection)
                overlap_fraction = float(intersection) / float(smaller_area) if smaller_area else 0.0
                if overlap_fraction < min_overlap_fraction:
                    continue
                links.append(
                    {
                        "z_from": z,
                        "z_to": z + dz,
                        "dz": dz,
                        "label_from": label_a,
                        "label_to": label_b,
                        "intersection_px": int(intersection),
                        "area_from_px": int(areas[label_a]),
                        "area_to_px": int(areas[label_b]),
                        "overlap_fraction": overlap_fraction,
                        "iou": float(intersection) / float(union) if union else 0.0,
                    }
                )
    return links


def stitch_masks(masks: np.ndarray, links: list[dict[str, object]]) -> tuple[np.ndarray, dict[int, int]]:
    labels = [int(label) for label in np.unique(masks) if label > 0]
    uf = UnionFind(labels)
    for link in links:
        uf.union(int(link["label_from"]), int(link["label_to"]))

    root_to_component: dict[int, int] = {}
    label_to_component: dict[int, int] = {}
    next_component = 1
    for label in labels:
        root = uf.find(label)
        if root not in root_to_component:
            root_to_component[root] = next_component
            next_component += 1
        label_to_component[label] = root_to_component[root]

    stitched = np.zeros_like(masks, dtype=np.uint32)
    for label, component in label_to_component.items():
        stitched[masks == label] = component
    return stitched, label_to_component


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_component_table(path: Path, sample_id: str, stitched: np.ndarray, label_to_component: dict[int, int]) -> None:
    stats = component_stats(stitched)
    slices_by_component: dict[int, set[int]] = defaultdict(set)
    labels_by_component: dict[int, list[int]] = defaultdict(list)
    zz, yy, xx = np.nonzero(stitched)
    for component in np.unique(stitched[zz, yy, xx]):
        if component > 0:
            slices_by_component[int(component)] = set(np.unique(zz[stitched[zz, yy, xx] == component]).astype(int))
    for label, component in label_to_component.items():
        labels_by_component[component].append(label)

    rows: list[dict[str, object]] = []
    for component in sorted(stats):
        labels = sorted(labels_by_component[component])
        z_slices = sorted(slices_by_component[component])
        rows.append(
            {
                "sample_id": sample_id,
                "stitched_object_id": component,
                "n_2d_masks": len(labels),
                "n_z_slices": len(z_slices),
                "z_min": int(stats[component]["z_min"]),
                "z_max": int(stats[component]["z_max"]),
                "z_labels_1_based": ",".join(str(z + 1) for z in z_slices),
                "area_px_total": int(stats[component]["area_px"]),
                "centroid_y": stats[component]["centroid_y"],
                "centroid_x": stats[component]["centroid_x"],
                "source_global_labels": ",".join(str(label) for label in labels),
            }
        )
    write_tsv(
        path,
        rows,
        [
            "sample_id",
            "stitched_object_id",
            "n_2d_masks",
            "n_z_slices",
            "z_min",
            "z_max",
            "z_labels_1_based",
            "area_px_total",
            "centroid_y",
            "centroid_x",
            "source_global_labels",
        ],
    )


def save_7x2_qc(path: Path, brightfield: np.ndarray, stitched: np.ndarray, z_labels: list[int]) -> None:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError as exc:
        raise RuntimeError("Pillow is required for QC PNG output") from exc

    z_indices = [z - 1 for z in z_labels]
    top_tiles = [np.repeat(robust_u8(brightfield[z])[..., None], 3, axis=2) for z in z_indices]
    bottom_tiles = [colorize_labels(stitched[z]) for z in z_indices]
    tile_h, tile_w = top_tiles[0].shape[:2]
    label_h = 22
    gap = 8
    ncols = len(z_indices)
    canvas_h = label_h + 2 * tile_h + gap
    canvas_w = ncols * tile_w + (ncols - 1) * gap
    canvas = np.full((canvas_h, canvas_w, 3), 255, dtype=np.uint8)

    for i, (top, bottom) in enumerate(zip(top_tiles, bottom_tiles)):
        x0 = i * (tile_w + gap)
        canvas[label_h : label_h + tile_h, x0 : x0 + tile_w] = top
        y1 = label_h + tile_h + gap
        canvas[y1 : y1 + tile_h, x0 : x0 + tile_w] = bottom

    img = Image.fromarray(canvas)
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    for i, z_label in enumerate(z_labels):
        x0 = i * (tile_w + gap)
        draw.text((x0 + 4, 4), f"z={z_label}", fill=(0, 0, 0), font=font)
        y_offset = label_h + tile_h + gap
        mask = stitched[z_indices[i]]
        for label in np.unique(mask):
            if label == 0:
                continue
            yy, xx = np.nonzero(mask == label)
            if yy.size == 0:
                continue
            text = str(int(label))
            tx = int(x0 + xx.mean())
            ty = int(y_offset + yy.mean())
            draw.text((tx - 1, ty), text, fill=(0, 0, 0), font=font)
            draw.text((tx + 1, ty), text, fill=(0, 0, 0), font=font)
            draw.text((tx, ty - 1), text, fill=(0, 0, 0), font=font)
            draw.text((tx, ty + 1), text, fill=(0, 0, 0), font=font)
            draw.text((tx, ty), text, fill=(255, 255, 255), font=font)

    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path)


def process_one(mask_path: Path, run_dir: Path, raw_crop_dir: Path, args: argparse.Namespace) -> dict[str, object]:
    sample_id = sample_id_from_mask(mask_path)
    raw_crop_path = raw_crop_dir / f"{sample_id}{RAW_CROP_SUFFIX}"
    if not raw_crop_path.exists():
        raise FileNotFoundError(raw_crop_path)

    masks = tf.imread(mask_path).astype(np.uint32, copy=False)
    raw = tf.imread(raw_crop_path)
    if raw.ndim != 4 or raw.shape[1] < 2:
        raise ValueError(f"Expected raw crop ZCYX with at least two channels, got {raw.shape}")
    if masks.shape != raw[:, 1].shape:
        raise ValueError(f"Mask stack shape {masks.shape} does not match brightfield shape {raw[:, 1].shape}")

    links = find_links(masks, max_dz=args.max_dz, min_overlap_fraction=args.min_overlap_fraction)
    stitched, label_to_component = stitch_masks(masks, links)

    out_dir = run_dir / "stitched_simple"
    mask_out = out_dir / "masks" / f"{sample_id}_{INPUT_KIND}_stitched_simple_masks.tif"
    links_out = out_dir / "links" / f"{sample_id}_{INPUT_KIND}_stitched_simple_links.tsv"
    objects_out = out_dir / "objects" / f"{sample_id}_{INPUT_KIND}_stitched_simple_objects.tsv"
    png_out = out_dir / "png" / f"{sample_id}_{INPUT_KIND}_stitched_simple_z10-70_7x2_bf_mask_ids.png"

    mask_out.parent.mkdir(parents=True, exist_ok=True)
    tf.imwrite(mask_out, stitched, metadata={"axes": "ZYX"})
    write_tsv(
        links_out,
        links,
        [
            "z_from",
            "z_to",
            "dz",
            "label_from",
            "label_to",
            "intersection_px",
            "area_from_px",
            "area_to_px",
            "overlap_fraction",
            "iou",
        ],
    )
    write_component_table(objects_out, sample_id, stitched, label_to_component)
    save_7x2_qc(png_out, raw[:, 1], stitched, args.z_labels)

    return {
        "sample_id": sample_id,
        "mask_path": mask_path,
        "stitched_mask_path": mask_out,
        "links_path": links_out,
        "objects_path": objects_out,
        "png_path": png_out,
        "n_input_2d_masks": len(label_to_component),
        "n_links": len(links),
        "n_stitched_objects": len(set(label_to_component.values())),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Simple first-pass stitching of globally labeled 2D full-stack masks into merged z objects."
    )
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--raw-crop-dir", type=Path, default=DEFAULT_RAW_CROP_DIR)
    parser.add_argument("--sample-id", default="", help="Process one sample by sample id.")
    parser.add_argument("--limit", type=int, default=0, help="Process the first N masks when --sample-id is not set.")
    parser.add_argument("--max-dz", type=int, default=1)
    parser.add_argument("--min-overlap-fraction", type=float, default=0.35)
    parser.add_argument("--z-labels", type=parse_int_list, default=parse_int_list("10,20,30,40,50,60,70"))
    args = parser.parse_args()

    mask_paths = discover_masks(args.run_dir, args.sample_id, args.limit)
    out_dir = args.run_dir / "stitched_simple"
    out_dir.mkdir(parents=True, exist_ok=True)
    config = {
        "run_dir": str(args.run_dir),
        "raw_crop_dir": str(args.raw_crop_dir),
        "max_dz": args.max_dz,
        "min_overlap_fraction": args.min_overlap_fraction,
        "link_score": "intersection_px / min(area_from_px, area_to_px)",
        "resolution": "All links above threshold are merged by connected components; many-to-one and one-to-many links are allowed.",
        "z_labels_1_based": args.z_labels,
    }
    (out_dir / "config.stitched_simple.json").write_text(json.dumps(config, indent=2) + "\n")

    summary_rows = [process_one(path, args.run_dir, args.raw_crop_dir, args) for path in mask_paths]
    write_tsv(
        out_dir / "sample_summary.tsv",
        summary_rows,
        [
            "sample_id",
            "mask_path",
            "stitched_mask_path",
            "links_path",
            "objects_path",
            "png_path",
            "n_input_2d_masks",
            "n_links",
            "n_stitched_objects",
        ],
    )
    for row in summary_rows:
        print(
            f"{row['sample_id']}: {row['n_input_2d_masks']} 2D masks, "
            f"{row['n_links']} links, {row['n_stitched_objects']} stitched objects"
        )


if __name__ == "__main__":
    main()
