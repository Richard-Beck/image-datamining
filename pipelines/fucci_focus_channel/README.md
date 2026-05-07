# FUCCI Focus Channel Exploration

This pipeline explores derived focus-support channels from FUCCI brightfield
z-stacks. It is separate from segmentation and writes visual review outputs
for small cropped regions.

## Dependencies

Requires `numpy`, `scipy`, and `tifffile`. PNG writing uses the first available
option among `imageio`, `Pillow`, or `matplotlib`.

## One-Command Example

```bash
python pipelines/fucci_focus_channel/bin/batch_focus_preview.py \
  --root /share/andor_lab/Jackson/FUCCI \
  --run-dir pipelines/fucci_focus_channel/runs/preview_crop256_ch01 \
  --crop-size 256 \
  --limit 6
```

By default this uses `ch01` as brightfield, computes focus confidence across z,
and writes two comparison figures per z-stack:

- central-slice 3x3 grid: raw `ch01` plus 8 relative focus-confidence methods
- z-sum projection 3x3 grid: raw `ch01` sum plus the same 8 methods

The default 8 tiles break out local standard deviation, Tenengrad, LoG, and
composite scores at a few scales. Only the two grid PNGs are written by
default. Add `--save-stacks` to also write the raw brightfield crop, the
first/default focus-confidence stack, and the first/default focus-score stack
as TIFFs for audit.

To compare only LoG scales:

```bash
python pipelines/fucci_focus_channel/bin/batch_focus_preview.py \
  --root /share/andor_lab/Jackson/FUCCI \
  --run-dir pipelines/fucci_focus_channel/runs/preview_crop256_ch01_log_scales \
  --crop-size 256 \
  --limit 6 \
  --preset log-scales
```

To try the local-SD plus large 3D Gaussian blur detector:

```bash
python pipelines/fucci_focus_channel/bin/batch_focus_preview.py \
  --root /share/andor_lab/Jackson/FUCCI \
  --run-dir pipelines/fucci_focus_channel/runs/preview_crop256_ch01_sd_blur \
  --crop-size 256 \
  --limit 6 \
  --preset sd-blur \
  --sd-window 12 \
  --blur-sigma-xy 20 \
  --z-anisotropy 1.44
```

This writes two 2x2 grids per stack:

- central z-slice
- z-sum projection

The tiles are raw `ch01`, blurred SD confidence, unblurred SD confidence, and
their product.

To evaluate detector scale crossed with blur scale, using LoG by default:

```bash
python pipelines/fucci_focus_channel/bin/batch_focus_preview.py \
  --root /share/andor_lab/Jackson/FUCCI \
  --run-dir pipelines/fucci_focus_channel/runs/preview_crop256_ch01_log_blur_grid \
  --crop-size 256 \
  --limit 6 \
  --preset detector-blur-grid \
  --grid-method log \
  --detector-scales 1,2,4 \
  --blur-sigmas-xy 5,10,20,30,40 \
  --z-anisotropy 1.44
```

This writes central-slice and z-sum projection grids. Rows are detector scale;
columns are blur scale. Each tile is:

```text
relative_focus_confidence * GaussianBlur3D(relative_focus_confidence)
```

To combine local SD and LoG with moderate 3D blur and Otsu gating:

```bash
python pipelines/fucci_focus_channel/bin/batch_focus_preview.py \
  --root /share/andor_lab/Jackson/FUCCI \
  --run-dir pipelines/fucci_focus_channel/runs/preview_crop256_ch01_sd_log_combo \
  --crop-size 256 \
  --limit 6 \
  --preset sd-log-combo \
  --combo-sd-window 10 \
  --combo-log-sigma 2 \
  --blur-sigma-xy 20 \
  --z-anisotropy 1.44
```

This preset writes preview grids and TIFF stacks for raw brightfield, SD
confidence, LoG confidence, combined confidence, blurred support, product, Otsu
mask on the blur, and Otsu-gated product.

To try global absolute z-score followed by large 3D smoothing:

```bash
python pipelines/fucci_focus_channel/bin/batch_focus_preview.py \
  --root /share/andor_lab/Jackson/FUCCI \
  --run-dir pipelines/fucci_focus_channel/runs/preview_crop256_ch01_abs_zscore_smooth \
  --crop-size 256 \
  --limit 6 \
  --preset abs-zscore-smooth \
  --blur-sigma-xy 20 \
  --z-anisotropy 1.44
```

This computes a global z-score over the cropped brightfield stack, takes
absolute values, smooths in 3D with anisotropy-aware sigma, and writes raw,
absolute z-score, smoothed absolute z-score, and Otsu mask previews plus TIFF
stacks.

Use `--single-setting` to disable the 8-setting comparison grid and use only
the explicitly supplied parameters:

```bash
python pipelines/fucci_focus_channel/bin/batch_focus_preview.py \
  --root /share/andor_lab/Jackson/FUCCI \
  --run-dir pipelines/fucci_focus_channel/runs/single_tenengrad_hp16 \
  --crop-size 256 \
  --limit 6 \
  --single-setting \
  --method tenengrad \
  --highpass-sigma 16
```
