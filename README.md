# image-datamining

This repo supports methods development and exploratory biology from existing lab microscopy data. The immediate goal is to find existing datasets in `/share/andor_lab` and `/share/lab_crd/lab_crd` that are suitable for two analysis directions:

1. Estimating volumes of different cell lines, which requires z-stack data.
2. Evaluating whether there is enough longitudinal data for cell tracking, which requires repeated observations at sufficiently frequent time intervals.

The drive-inventory tools live in `DriveExplorer/`. They map directory structure, count direct child entries, sample filenames, identify timeout-prone image-heavy folders, read cheap TIFF stack metadata when available, and rank candidate image directories for follow-up. At this stage, follow-up should prioritize filename, path, and lightweight metadata evidence for z, time, channel, well, field/site, and instrument-specific acquisition structure. Expensive steps, such as segmentation, volume measurement, or tracking, should run only after this triage identifies promising datasets.

## Current Inventory

- `data_inventory/andor_lab/`: inventory outputs from `/share/andor_lab`.
- `data_inventory/lab_crd/`: inventory outputs from `/share/lab_crd/lab_crd`.
- `data/JacksonIncucyte_230623_Chemotaxis/`: curated notes and extracted plate map for the 2023-06-23 chemotaxis experiment.

The walk uses per-directory timeouts. Directories that time out are marked `timeout_blackhole`, sampled for filenames, and not explored further.

## Triage Priorities

Filename and path parsing should focus on evidence that a directory supports one of the two analysis goals:

- Z-stack / volume candidates: filenames or paths with z indices, planes, slices, channels, and repeated fields. Examples include Zeiss-style `zNN` naming and OperaPhenix-style plane/position tokens such as `pNN`.
- Longitudinal / tracking candidates: filenames or paths with explicit time indices, elapsed time, acquisition dates, or repeated timestamped observations. Examples include Incucyte-style `YYYYyMMmDDd_HHhMMm`, elapsed-time names such as `04d12h00m`, and Zeiss-style `tNNN` naming.
- Cell line and experiment context: parent directories that encode cell line names, treatment conditions, assay names, or experiment dates.
- Raw-vs-derived status: distinguish raw image folders from masks, segmentation outputs, transformed images, bounding boxes, and other processed derivatives.

The first parsing pass uses existing filepaths, sampled filenames, and cheap TIFF metadata. It does not read image pixels.

`DriveExplorer/scripts/profile_acquisition_candidates.py` implements this pass. It scores directories using transparent filename/path/metadata evidence columns rather than attempting to fully parse every instrument format. Current evidence includes image-like extensions, explicit `zNN` tokens, OperaPhenix-style plane tokens, Zeiss-style `tNNN` tokens, elapsed-time tokens such as `04d12h00m`, Incucyte-style acquisition timestamps, well/site structure, channel labels, cell-line-like path tokens, raw-vs-derived path hints, cheap TIFF stack metadata, and tracking-focused path matches against `chemotaxis`, `migration`, `motility`, and `tracking`.

The candidate generator excludes known screening/control folders that are not useful for the current volume/tracking mining goals. Current exclusions include path tokens matching `myco`.

The z-stack candidate output also records sampled channel multiplicity with `distinct_channels`, `has_multiple_channels`, and `channel_z_or_plane_coverage`.

`DriveExplorer/scripts/summarize_z_stack_experiments.py` groups z-stack candidates into likely experiment roots so sibling fields, sites, channels, and condition folders can be reviewed as coherent volume-analysis datasets.

`DriveExplorer/scripts/summarize_tracking_experiments.py` groups the focused tracking candidates into likely experiment roots so sibling channel, modality, and condition folders can be evaluated together.

`DriveExplorer/scripts/summarize_tracking_csvs.py` summarizes existing tracking CSV exports by field, including frame coverage, track counts, track-length thresholds, objects per frame, lineage/link counts, merger labels, and missingness in key columns.

## Regenerate

See `DriveExplorer/README.md` for full regeneration commands for both source drives.

```bash
python3 DriveExplorer/scripts/summarize_tracking_csvs.py \
  /share/andor_lab/Jackson/FUCCI_TimepointAnalysis/QI_Core_Analysis_Results/Tracking_Data/Tiffs \
  --output data_inventory/andor_lab/fucci_tracking_csv_qc.tsv
```
