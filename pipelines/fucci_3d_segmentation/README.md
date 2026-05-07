# FUCCI 3D Segmentation Pipeline

This pipeline is for Jackson FUCCI z-stack segmentation experiments. It is
manifest-driven and keeps segmentation, mask collection, summarization, and
visualization as separate stages.

## Scope

- Inputs are z-stacks with split channels `ch00`, `ch01`, and `ch02`.
- 3D segmentation uses the full z-stack as a single volume.
- 2D segmentation means segmenting each z-plane independently, not segmenting
  max projections.
- Pixel-to-object mask TIFFs are the primary durable segmentation artifact.
- Visualization is a later standalone review step.

## Layout

```text
pipelines/fucci_3d_segmentation/
  configs/
  manifests/
  runs/
  bin/
```

Reusable implementation lives in `src/image_datamining/fucci_3d_segmentation/`.
The files in `bin/` are thin command-line wrappers for this pipeline.

## Dependencies

Manifest and command generation use the Python standard library. Input
preparation and mask summarization require `numpy` and `tifffile`; FUCCI fused
preprocessing also requires `scipy`. Running generated CellposeSAM commands
requires `cellpose` in the execution environment.

## One-Command Streaming Run

This is the preferred route for CellposeSAM experiments because it processes
one field at a time:

1. discover z-stack fields
2. load one field
3. make the centered crop
4. run CellposeSAM
5. write the pixel-to-object mask TIFF
6. continue to the next field

```bash
python pipelines/fucci_3d_segmentation/bin/batch_segment.py \
  --root /share/andor_lab/Jackson/FUCCI \
  --run-dir pipelines/fucci_3d_segmentation/runs/jackson_fucci_crop256_bf_fucci_fused_3d \
  --mode 3d \
  --crop-size 256 \
  --resume
```

By default this writes all three audit artifacts per field:

- `raw_crops/`: 3-channel cropped `ZCYX` TIFFs
- `cellpose_inputs/`: 2-channel merged Cellpose input `ZCYX` TIFFs
- `masks/`: pixel-to-object segmentation mask TIFFs

Use `--discard-prepared-inputs` only for a temporary low-disk run. Use
`--plan-only` to write the manifest, config, and planned output paths without
importing or running CellposeSAM.

## FUCCI/BF/Focus Product Input

To build 3-channel CellposeSAM inputs with FUCCI merge, raw brightfield, and
the product of global absolute z-score with its large 3D smooth:

```bash
python pipelines/fucci_3d_segmentation/bin/batch_segment.py \
  --root /share/andor_lab/Jackson/FUCCI \
  --run-dir pipelines/fucci_3d_segmentation/runs/jackson_fucci_crop256_fucci_bf_absz_product_3d \
  --mode 3d \
  --crop-size 256 \
  --input-kind fucci_bf_absz_product \
  --focus-blur-sigma-xy 20 \
  --focus-z-anisotropy 1.44 \
  --anisotropy 1.44 \
  --resume
```

This still writes the original 3-channel crop under `raw_crops/`, the derived
3-channel model input under `cellpose_inputs/`, and the pixel-to-object mask
under `masks/`.

## Staged Flow

```bash
python pipelines/fucci_3d_segmentation/bin/make_manifest.py \
  --root /share/andor_lab/Jackson/FUCCI \
  --recursive \
  --output pipelines/fucci_3d_segmentation/manifests/jackson_fucci_zstacks.tsv

python pipelines/fucci_3d_segmentation/bin/prepare_inputs.py \
  --manifest pipelines/fucci_3d_segmentation/manifests/jackson_fucci_zstacks.tsv \
  --run-dir pipelines/fucci_3d_segmentation/runs/example_crop256 \
  --crop-size 256

python pipelines/fucci_3d_segmentation/bin/write_cellpose_commands.py \
  --prepared-manifest pipelines/fucci_3d_segmentation/runs/example_crop256/prepared_inputs.tsv \
  --run-dir pipelines/fucci_3d_segmentation/runs/example_crop256 \
  --mode 3d
```

The command manifest does not run CellposeSAM. It records commands for manual
execution on an appropriate machine.
