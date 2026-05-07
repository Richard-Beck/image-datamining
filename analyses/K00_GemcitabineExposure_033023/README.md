# K00 Gemcitabine Exposure Tracking Analysis

This analysis evaluates pre-existing SUM-159 gemcitabine exposure tracking
results. The folder has been aggressively cleaned to keep only the current
preprocessing, isolated-track, OU-tracking, data-overview, and overlay-helper
workflow.

## Input Sources

Tracking CSVs:

```text
/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Tracking_CSVs
```

Registered image stacks for visualization:

```text
/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Images_40Frames
```

Use the registered `Images_40Frames` stacks for coordinate overlays. The
tracking coordinates align to that registered image space, not directly to the
canonical raw TIFFs.

## Kept Files

- `Gemcitabine_PlateMap_20240111.xlsx`: source treatment map.
- `platemap_conditions.csv`: parsed plate conditions.
- `data/tracking_data.rds`: full loaded tracking data with metadata and
  isolation summaries.
- `data/tracking_data_isolated_2x_min3.rds`: isolated tracks using a 2x
  isolation threshold and minimum 3-frame contiguous segments.
- `data/tracking_summaries.rds`: frame-level tracking summaries.
- `data/isolation_summaries.rds`: isolation-threshold summaries.
- `data/migration_summaries_isolated_2x_min3.rds`: migration summaries from the
  filtered isolated tracks.
- `data/ou_tracking_manifest.csv`: condition-level OU fit manifest.
- `data/ou_tracking_fits.csv`: condition-level OU fit outputs.
- `data_overview.Rmd` and `data_overview.html`: current data overview report.
- `ou_tracking_fit_report.Rmd` and `ou_tracking_fit_report.html`: rough report
  plotting condition-level OU fit results.
- `R/preprocessing.R`: tracking CSV loading, platemap parsing, and frame
  summaries.
- `R/migration_stats.R`: isolated-track filtering and migration summaries.
- `R/ou_velocity_model.R`: OU velocity model likelihood and fitting helpers.
- `R/overlays.R`: reusable PNG/GIF tracking overlay helpers.
- `R/k00_batch_utils.R`: shared CLI and file helper utilities.
- `workflow/load_preprocess_tracks.R`: builds `tracking_data.rds`,
  `tracking_summaries.rds`, and `isolation_summaries.rds`.
- `workflow/filter_and_summarize_migration.R`: builds
  `tracking_data_isolated_2x_min3.rds` and
  `migration_summaries_isolated_2x_min3.rds`.
- `workflow/make_ou_tracking_manifest.R`: builds the OU condition manifest.
- `workflow/fit_ou_tracking_manifest_row.R`: fits one OU manifest row.
- `workflow/submit_ou_tracking_manifest_array.slurm`: Slurm launcher for the OU
  manifest.
- `workflow/benchmark_ou_velocity_rcpp.R`: validates/times R vs Rcpp OU
  likelihood implementations.
- `workflow/render_data_overview.R`: renders `data_overview.Rmd`.
- `workflow/render_ou_tracking_fit_report.R`: renders
  `ou_tracking_fit_report.Rmd`.
- `workflow/render_ou_fit_sanity_overlays.R`: renders selected GIF overlays for
  rough visual QC of OU fit contrasts.

## Regenerate Current Data

Run through the repository R container from the repository root:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/load_preprocess_tracks.R
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/filter_and_summarize_migration.R
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/make_ou_tracking_manifest.R
```

Fit the condition-level OU manifest locally for one row:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/fit_ou_tracking_manifest_row.R \
  --file_index=1
```

Or submit the array:

```bash
sbatch analyses/K00_GemcitabineExposure_033023/workflow/submit_ou_tracking_manifest_array.slurm
```

Render the data overview:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/render_data_overview.R
```

Render the OU tracking fit report:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/render_ou_tracking_fit_report.R
```

## OU Tracking Fits

`data/ou_tracking_fits.csv` is generated from the isolated-track data, not from
the raw tracking CSVs directly. The rough lineage is:

```text
Tracking_CSVs + Gemcitabine_PlateMap_20240111.xlsx
  -> workflow/load_preprocess_tracks.R
  -> data/tracking_data.rds
  -> workflow/filter_and_summarize_migration.R
  -> data/tracking_data_isolated_2x_min3.rds
  -> workflow/make_ou_tracking_manifest.R
  -> data/ou_tracking_manifest.csv
  -> workflow/fit_ou_tracking_manifest_row.R or submit_ou_tracking_manifest_array.slurm
  -> data/ou_tracking_fits.csv
```

The manifest has one row per `ploidy`/`Gemcitabine` condition. In the current
data this is 20 rows: 2 ploidy states by 10 gemcitabine doses. The Slurm
launcher fans this out across starts with `#SBATCH --array=1-500`, corresponding
to 20 conditions by 25 starts. If `N_STARTS` or the manifest row count changes,
submit with an explicit Slurm array range, for example
`sbatch --array=1-1000 ...` for 20 conditions by 50 starts.

Current default OU fit settings are written into the manifest:

- `frame_interval = 2`
- `min_segment_frames = 3`
- `max_tracks = 0`, meaning no track subsampling
- `max_segments = 0`, meaning no segment subsampling
- `n_starts = 25`, meaning one empirical start plus 24 seeded random starts
- `seed = 17`, used for random starts and any optional subsampling

Each fit job filters `tracking_data_isolated_2x_min3.rds` to one manifest
condition, optionally applies frame/window/subsampling filters, re-splits tracks
into contiguous post-filter segments, drops segments shorter than
`min_segment_frames`, converts those segments to OU velocity-model inputs, fits
one requested start, and appends one output row to `data/ou_tracking_fits.csv`.
With the default Slurm launcher, `data/ou_tracking_fits.csv` should have 500
rows plus a header after all tasks finish. The report selects the
highest-likelihood start for each condition before plotting condition-level
parameters.

Important rough-edge note: `fit_ou_tracking_manifest_row.R` appends to
`data/ou_tracking_fits.csv`; it does not clear old results. For a clean rerun,
remove or rename the old output CSV first, or pass a different `--out_csv` when
creating the manifest. Row order in the finished CSV is job-completion order,
not manifest order, so sort by `job_id` when comparing conditions.

The output columns include the manifest metadata plus OU fit results:
`tau`, `velocity_scale`, `effective_diffusivity`, `obs_noise`,
`log_likelihood`, `fit_status`, `n_segments`, and `total_track_time`. Multistart
provenance columns include `n_starts`, `seed`, `start_id`, `start_tau`,
`start_velocity_scale`, and `start_obs_noise`. Slurm task logs are written under
`logs/ou_tracking_manifest/condition_${FILE_INDEX}_start_${START_ID}.log`.

## Tracking Overlay Helpers

Reusable registered-frame overlay tools live in:

```text
analyses/K00_GemcitabineExposure_033023/R/overlays.R
```

Use these helpers for visual QC of tracking data in the registered
`Final_Tracking_analysis/Images_40Frames` coordinate space. The input tracking
object should already be filtered to one well/position site and should use the
same column names as `data/tracking_data_isolated_2x_min3.rds`. By default,
coordinates are interpreted as `Center_of_the_object_0 = x` and
`Center_of_the_object_1 = y`.

Render one PNG overlay:

```r
source("analyses/K00_GemcitabineExposure_033023/R/overlays.R")

tracks <- readRDS("analyses/K00_GemcitabineExposure_033023/data/tracking_data_isolated_2x_min3.rds")
site_tracks <- subset(tracks, well == "E7" & position == 2 & frame %in% 11:13)

tracking_overlay_png(
  tracks = site_tracks,
  images_dir = "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Images_40Frames",
  out_png = "analyses/K00_GemcitabineExposure_033023/overlay_checks/E7_2_frame012.png",
  display_frame = 12,
  arrow_mode = "all"
)
```

`tracking_overlay_png()` draws points only for `display_frame`. Arrows are
controlled by `arrow_mode`: use `"all"` for every consecutive step in the
passed tracks, `"current_to_next"` for arrows beginning on `display_frame`, or
`"none"` for point-only overlays. Cropping can be set with `crop = "tracks"` and
`crop_buffer_px`, or with explicit pixel bounds:

```r
crop_bounds = c(xmin = 100, xmax = 400, ymin = 250, ymax = 550)
```

Render a PNG sequence with stable crop bounds and convert it to a GIF:

```r
source("analyses/K00_GemcitabineExposure_033023/R/overlays.R")

tracks <- readRDS("analyses/K00_GemcitabineExposure_033023/data/tracking_data_isolated_2x_min3.rds")
one_track <- subset(tracks, migration_track_id == "A10_1_108_1")

frames <- tracking_overlay_png_sequence(
  tracks = one_track,
  images_dir = "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Images_40Frames",
  out_dir = "analyses/K00_GemcitabineExposure_033023/overlay_checks/A10_1_108_1_frames",
  crop = "tracks",
  crop_buffer_px = 30,
  arrow_mode = "current_to_next"
)

png_sequence_gif(
  frame_paths = frames$frame_paths,
  out_gif = "analyses/K00_GemcitabineExposure_033023/overlay_checks/A10_1_108_1.gif",
  fps = 2,
  cleanup_frames = FALSE
)
```

`tracking_overlay_png_sequence()` defaults to rendering all frames present in
the supplied tracks. It computes one crop window across all supplied
coordinates, so the GIF does not jump from frame to frame. Set
`cleanup_frames = TRUE` in `png_sequence_gif()` when the intermediate PNGs are
not needed.

## OU Fit Visual Sanity Check

The current rough visual QC script is:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/render_ou_fit_sanity_overlays.R
```

It reads the current `data/ou_tracking_fits.csv`, selects the best
highest-likelihood start per condition, ranks 2N-vs-4N differences by
`effective_diffusivity`, and renders GIF overlays for:

- the 0 nM gemcitabine condition
- the nonzero dose with the largest 2N-vs-4N fitted difference

For each selected dose, it renders low/median/high path-length examples for
both 2N and 4N isolated tracks. Outputs are written here:

```text
analyses/K00_GemcitabineExposure_033023/overlay_checks/ou_fit_sanity
```

The script writes one flat GIF folder per dose, plus:

- `ou_fit_2N_4N_contrasts.csv`
- `selected_tracks.csv`
- `rendered_overlays.csv`

Current note from visual inspection: the `dose_100` 4N overlays revealed clear
tracking errors consistent with tracking label switching. This means at least
some large fitted 4N motility differences may be driven or inflated by tracking
artifacts rather than true cell motion.

Next step: add more rigorous filtering and QC of the tracking data before
relying on the OU fits. In particular, filter or flag likely label-switching
events using track-level motion diagnostics, such as implausibly large single
frame jumps, large path length with modest net displacement, abrupt angle
changes, site-level outlier rates, and/or manual overlay review of tracks that
dominate high fitted motility conditions.
