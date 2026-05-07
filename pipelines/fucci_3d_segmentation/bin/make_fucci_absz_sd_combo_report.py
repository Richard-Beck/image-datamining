#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import html
import re
from collections import defaultdict
from pathlib import Path


DEFAULT_RUN_DIR = Path("pipelines/fucci_3d_segmentation/runs/fucci_absz_sd_combo_cpsam_2d_grids")

METHOD_DESCRIPTIONS = {
    "fucci_absz": "FUCCI fused + abs z-score",
    "fucci_sd_w10": "FUCCI fused + local SD detector",
    "fucci_absz_sd_mean": "FUCCI fused + mean(abs z-score, SD)",
    "fucci_absz_sd_geom": "FUCCI fused + geometric mean(abs z-score, SD)",
    "fucci_absz_sd_product": "FUCCI fused + product(abs z-score, SD)",
    "fucci_absz_sd_min": "FUCCI fused + min(abs z-score, SD)",
    "fucci_absz_sd_max": "FUCCI fused + max(abs z-score, SD)",
    "fucci_absz2_sd_product": "FUCCI fused + product(abs z-score squared, SD)",
    "fucci_absz_sd2_product": "FUCCI fused + product(abs z-score, SD squared)",
    "fucci_bf_absz": "FUCCI fused + brightfield + abs z-score",
    "fucci_bf_sd_w10": "FUCCI fused + brightfield + local SD detector",
    "fucci_bf_absz_sd_mean": "FUCCI fused + brightfield + mean(abs z-score, SD)",
    "fucci_bf_absz_sd_geom": "FUCCI fused + brightfield + geometric mean(abs z-score, SD)",
    "fucci_bf_absz_sd_product": "FUCCI fused + brightfield + product(abs z-score, SD)",
    "fucci_bf_absz_sd_min": "FUCCI fused + brightfield + min(abs z-score, SD)",
    "fucci_bf_absz_sd_max": "FUCCI fused + brightfield + max(abs z-score, SD)",
    "fucci_bf_absz2_sd_product": "FUCCI fused + brightfield + product(abs z-score squared, SD)",
    "fucci_bf_absz_sd2_product": "FUCCI fused + brightfield + product(abs z-score, SD squared)",
}


METHOD_ORDER = [
    "fucci_absz",
    "fucci_sd_w10",
    "fucci_absz_sd_mean",
    "fucci_absz_sd_geom",
    "fucci_absz_sd_product",
    "fucci_absz_sd_min",
    "fucci_absz_sd_max",
    "fucci_absz2_sd_product",
    "fucci_absz_sd2_product",
    "fucci_bf_absz",
    "fucci_bf_sd_w10",
    "fucci_bf_absz_sd_mean",
    "fucci_bf_absz_sd_geom",
    "fucci_bf_absz_sd_product",
    "fucci_bf_absz_sd_min",
    "fucci_bf_absz_sd_max",
    "fucci_bf_absz2_sd_product",
    "fucci_bf_absz_sd2_product",
]


def sample_from_png(path: Path, variant: str) -> str:
    suffix = f"_{variant}_z10-70_2x7_bf_mask.png"
    if path.name.endswith(suffix):
        return path.name[: -len(suffix)]
    return re.sub(r"_z\d+-\d+_2x\d+_bf_mask\.png$", "", path.stem)


def rows_from_summary(run_dir: Path) -> list[dict[str, str]]:
    summary = run_dir / "sample_variant_summary.tsv"
    if not summary.exists():
        return []
    rows: list[dict[str, str]] = []
    with summary.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            png = Path(row.get("png_path", ""))
            if row.get("status") in {"ok", "skipped_existing"} and png.exists():
                rows.append(
                    {
                        "sample_id": row["sample_id"],
                        "variant": row["variant"],
                        "png_path": str(png),
                    }
                )
    return rows


def rows_from_pngs(run_dir: Path) -> list[dict[str, str]]:
    rows = []
    for png in sorted((run_dir / "png").glob("*/*.png")):
        variant = png.parent.name
        rows.append(
            {
                "sample_id": sample_from_png(png, variant),
                "variant": variant,
                "png_path": str(png),
            }
        )
    return rows


def rel(path: Path, start: Path) -> str:
    return path.resolve().relative_to(start.resolve()).as_posix()


def method_sort_key(variant: str) -> tuple[int, str]:
    try:
        return METHOD_ORDER.index(variant), variant
    except ValueError:
        return len(METHOD_ORDER), variant


def build_html(run_dir: Path, output: Path, rows: list[dict[str, str]]) -> str:
    grouped: dict[str, dict[str, str]] = defaultdict(dict)
    for row in rows:
        grouped[row["sample_id"]][row["variant"]] = row["png_path"]

    sample_ids = sorted(grouped)
    total_images = sum(len(v) for v in grouped.values())
    sections = []
    for sample_id in sample_ids:
        methods = sorted(grouped[sample_id], key=method_sort_key)
        cards = []
        for variant in methods:
            png_path = Path(grouped[sample_id][variant])
            src = html.escape(rel(png_path, output.parent))
            title = html.escape(METHOD_DESCRIPTIONS.get(variant, variant.replace("_", " ")))
            variant_text = html.escape(variant)
            cards.append(
                f"""
        <article class="method">
          <div class="method-title">{title}</div>
          <div class="method-code"><code>{variant_text}</code></div>
          <img src="{src}" alt="{html.escape(sample_id)} {variant_text}">
        </article>"""
            )
        sections.append(
            f"""
      <section class="sample">
        <h2>{html.escape(sample_id)}</h2>
        <div class="methods">
          {''.join(cards)}
        </div>
      </section>"""
        )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>FUCCI Abs Z-Score / SD Combo CPSAM Grids</title>
  <style>
    :root {{
      --ink: #1f2528;
      --muted: #5f6b70;
      --line: #d8dee2;
      --bg: #f6f7f8;
      --panel: #ffffff;
      --accent: #145f66;
    }}
    body {{
      margin: 0;
      background: var(--bg);
      color: var(--ink);
      font: 15px/1.45 system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    header {{
      position: sticky;
      top: 0;
      z-index: 5;
      background: rgba(246, 247, 248, 0.96);
      border-bottom: 1px solid var(--line);
      padding: 14px 22px;
      backdrop-filter: blur(8px);
    }}
    h1 {{
      margin: 0 0 4px;
      font-size: 21px;
      letter-spacing: 0;
    }}
    .meta {{
      color: var(--muted);
      max-width: 1180px;
    }}
    main {{
      max-width: 1440px;
      margin: 0 auto;
      padding: 18px 22px 44px;
    }}
    .sample {{
      margin: 0 0 26px;
      border: 1px solid var(--line);
      border-radius: 8px;
      background: var(--panel);
      overflow: hidden;
    }}
    h2 {{
      position: sticky;
      top: 75px;
      z-index: 2;
      margin: 0;
      padding: 10px 14px;
      font-size: 17px;
      letter-spacing: 0;
      background: #ffffff;
      border-bottom: 1px solid var(--line);
    }}
    .methods {{
      padding: 12px;
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(520px, 1fr));
      gap: 14px;
    }}
    .method {{
      border: 1px solid var(--line);
      border-radius: 7px;
      overflow: hidden;
      background: #fff;
    }}
    .method-title {{
      padding: 9px 10px 2px;
      font-weight: 700;
    }}
    .method-code {{
      padding: 0 10px 9px;
      color: var(--muted);
      font-size: 12px;
    }}
    code {{
      font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
      background: #eef2f3;
      border: 1px solid #e2e7e9;
      border-radius: 4px;
      padding: 1px 4px;
    }}
    img {{
      display: block;
      width: 100%;
      height: auto;
      border-top: 1px solid var(--line);
      background: #fff;
    }}
    @media (max-width: 700px) {{
      .methods {{
        grid-template-columns: 1fr;
      }}
      h2 {{
        top: 98px;
      }}
    }}
  </style>
</head>
<body>
  <header>
    <h1>FUCCI Abs Z-Score / SD Combo CPSAM Grids</h1>
    <div class="meta">
      Grouped by image first, then method. Each PNG shows z=10,20,...,70 with raw ch01 brightfield on top and 2D CPSAM masks on bottom.
      Included {len(sample_ids)} images and {total_images} method grids from <code>{html.escape(str(run_dir))}</code>.
    </div>
  </header>
  <main>
    {''.join(sections)}
  </main>
</body>
</html>
"""


def main() -> None:
    parser = argparse.ArgumentParser(description="Build image-first HTML report for FUCCI abs-zscore/SD combo PNG grids.")
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--output", type=Path, help="HTML output path. Defaults to RUN_DIR/report_by_image.html.")
    parser.add_argument(
        "--from-pngs",
        action="store_true",
        help="Ignore sample_variant_summary.tsv and discover PNGs directly.",
    )
    args = parser.parse_args()

    output = args.output or args.run_dir / "report_by_image.html"
    rows = [] if args.from_pngs else rows_from_summary(args.run_dir)
    if not rows:
        rows = rows_from_pngs(args.run_dir)
    if not rows:
        raise ValueError(f"No PNG rows found under {args.run_dir}")

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_html(args.run_dir, output, rows))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
