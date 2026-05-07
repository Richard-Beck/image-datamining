# DriveExplorer

DriveExplorer is the lightweight shared-drive inventory and triage toolkit for this repository.
It is designed for large microscopy folders where full image reads are too expensive for a
first pass.

## What It Does

- Walks a directory tree with per-directory timeouts.
- Records direct-child counts, extension counts, and sampled filenames.
- Ranks image-heavy directories.
- Scores candidate z-stack / volume datasets from filename, path, and cheap TIFF metadata evidence.
- Scores candidate longitudinal / tracking datasets from filename and path evidence.
- Collapses per-directory candidates into likely experiment-level groups.

The profiler does not read image pixels. For TIFF stack detection it first reads only the first
page `ImageDescription` and parses OME `SizeZ`, `SizeT`, `SizeC`, `SizeX`, and `SizeY`. If OME
sizes are absent, it falls back to `len(tif.pages)` and marks those hits as ambiguous multipage
TIFFs.

## Scripts

- `scripts/walk_directory_tree.py`: bounded directory walk with timeout black-hole handling.
- `scripts/rank_image_dirs.py`: rank directories by image-like file counts.
- `scripts/profile_acquisition_candidates.py`: generate z-stack and tracking candidate TSVs.
- `scripts/summarize_z_stack_experiments.py`: group z-stack candidates into experiment roots.
- `scripts/summarize_tracking_experiments.py`: group tracking candidates into experiment roots.
- `scripts/summarize_tracking_csvs.py`: summarize existing tracking CSV exports.
- `scripts/time_tiff_metadata_reads.py`: benchmark TIFF metadata access patterns.

## Current Output Layout

Inventory outputs are grouped by source drive:

```text
data_inventory/
  andor_lab/
  lab_crd/
```

Each source folder contains the tree walk, ranked image directories, candidate TSVs, and grouped
experiment summaries generated for that drive.

## Regenerate `/share/andor_lab`

```bash
python3 DriveExplorer/scripts/walk_directory_tree.py /share/andor_lab \
  --jsonl-output data_inventory/andor_lab/andor_lab_tree_walk.jsonl \
  --tsv-output data_inventory/andor_lab/andor_lab_tree_walk.tsv \
  --scan-timeout 2 \
  --sample-timeout 2 \
  --sample-limit 50 \
  --progress-every 100

python3 DriveExplorer/scripts/rank_image_dirs.py \
  data_inventory/andor_lab/andor_lab_tree_walk.jsonl \
  --output data_inventory/andor_lab/image_dirs_by_count.tsv

python3 DriveExplorer/scripts/profile_acquisition_candidates.py \
  data_inventory/andor_lab/andor_lab_tree_walk.jsonl \
  --z-output data_inventory/andor_lab/z_stack_candidates.tsv \
  --tracking-output data_inventory/andor_lab/tracking_candidates.tsv \
  --focused-tracking-output data_inventory/andor_lab/focused_tracking_candidates.tsv

python3 DriveExplorer/scripts/summarize_z_stack_experiments.py \
  data_inventory/andor_lab/z_stack_candidates.tsv \
  --output data_inventory/andor_lab/z_stack_experiment_groups.tsv

python3 DriveExplorer/scripts/summarize_tracking_experiments.py \
  data_inventory/andor_lab/focused_tracking_candidates.tsv \
  --output data_inventory/andor_lab/tracking_experiment_groups.tsv
```

## Regenerate `/share/lab_crd/lab_crd`

```bash
python3 DriveExplorer/scripts/walk_directory_tree.py /share/lab_crd/lab_crd \
  --jsonl-output data_inventory/lab_crd/lab_crd_tree_walk.jsonl \
  --tsv-output data_inventory/lab_crd/lab_crd_tree_walk.tsv \
  --scan-timeout 2 \
  --sample-timeout 2 \
  --sample-limit 50 \
  --progress-every 100

python3 DriveExplorer/scripts/rank_image_dirs.py \
  data_inventory/lab_crd/lab_crd_tree_walk.jsonl \
  --output data_inventory/lab_crd/lab_crd_image_dirs_by_count.tsv

python3 DriveExplorer/scripts/profile_acquisition_candidates.py \
  data_inventory/lab_crd/lab_crd_tree_walk.jsonl \
  --z-output data_inventory/lab_crd/lab_crd_z_stack_candidates.tsv \
  --tracking-output data_inventory/lab_crd/lab_crd_tracking_candidates.tsv \
  --focused-tracking-output data_inventory/lab_crd/lab_crd_focused_tracking_candidates.tsv

python3 DriveExplorer/scripts/summarize_z_stack_experiments.py \
  data_inventory/lab_crd/lab_crd_z_stack_candidates.tsv \
  --output data_inventory/lab_crd/lab_crd_z_stack_experiment_groups.tsv

python3 DriveExplorer/scripts/summarize_tracking_experiments.py \
  data_inventory/lab_crd/lab_crd_focused_tracking_candidates.tsv \
  --output data_inventory/lab_crd/lab_crd_tracking_experiment_groups.tsv
```

