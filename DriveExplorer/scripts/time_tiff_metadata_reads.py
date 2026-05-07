#!/usr/bin/env python3
"""Time progressively more expensive TIFF metadata access patterns."""

from __future__ import annotations

import argparse
import re
import statistics
import time
from pathlib import Path

import tifffile


OME_SIZE_PATTERN = re.compile(r'\b(Size[ZTCXY])="([^"]+)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark TIFF metadata access patterns.")
    parser.add_argument("tiff", help="TIFF/OME-TIFF path to inspect.")
    parser.add_argument("--repeat", type=int, default=3, help="Repeats per access pattern.")
    parser.add_argument(
        "--include-pixel-read",
        action="store_true",
        help="Also time reading one small pixel slice. This reads image data.",
    )
    return parser.parse_args()


def time_call(label: str, repeat: int, func) -> None:
    durations: list[float] = []
    last_result = None
    for _ in range(repeat):
        started = time.perf_counter()
        last_result = func()
        durations.append(time.perf_counter() - started)

    print(
        "\t".join(
            [
                label,
                f"min={min(durations):.4f}s",
                f"median={statistics.median(durations):.4f}s",
                f"max={max(durations):.4f}s",
                f"result={last_result}",
            ]
        )
    )


def open_close(path: Path) -> str:
    with tifffile.TiffFile(path):
        return "opened"


def first_page_basic(path: Path) -> str:
    with tifffile.TiffFile(path) as tif:
        page0 = tif.pages[0]
        return f"shape={page0.shape}, dtype={page0.dtype}, samples={getattr(page0, 'samplesperpixel', None)}"


def first_page_description(path: Path) -> str:
    with tifffile.TiffFile(path) as tif:
        desc = tif.pages[0].description or ""
        return f"description_chars={len(desc)}, ome={'OME' in desc[:256]}"


def first_page_ome_sizes(path: Path) -> str:
    with tifffile.TiffFile(path) as tif:
        desc = tif.pages[0].description or ""
        sizes = dict(OME_SIZE_PATTERN.findall(desc))
        return ", ".join(f"{key}={value}" for key, value in sorted(sizes.items())) or "no OME sizes"


def series_metadata(path: Path) -> str:
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        return f"shape={series.shape}, axes={series.axes}, dtype={series.dtype}"


def page_count(path: Path) -> str:
    with tifffile.TiffFile(path) as tif:
        return f"pages={len(tif.pages)}"


def first_pixel_page(path: Path) -> str:
    with tifffile.TiffFile(path) as tif:
        arr = tif.pages[0].asarray()
        return f"shape={arr.shape}, dtype={arr.dtype}"


def main() -> int:
    args = parse_args()
    path = Path(args.tiff)
    if not path.is_file():
        raise SystemExit(f"not a file: {path}")

    tests = [
        ("open_close", open_close),
        ("first_page_basic", first_page_basic),
        ("first_page_description", first_page_description),
        ("first_page_ome_sizes", first_page_ome_sizes),
        ("series_metadata", series_metadata),
        ("page_count", page_count),
    ]
    if args.include_pixel_read:
        tests.append(("first_pixel_page", first_pixel_page))

    print(f"path={path}")
    print(f"repeat={args.repeat}")
    for label, func in tests:
        time_call(label, args.repeat, lambda func=func: func(path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
