# K00 Gemcitabine Exposure Tracking Analysis

This directory now keeps the current yellow reconstructed CellposeSAM/nucleus
tracking analysis for the SUM-159 gemcitabine exposure experiment. Earlier
nucleus-only fitting outputs, btrack prototypes, temporary smoke tests, and
superseded comparison products have been removed.

## Active Path

The primary migration analysis path is:

```text
cpsam_nucleus_tracking_links/yellow_track_segments_len3
  -> workflow/build_yellow_reconstructed_ou_tracks.R
  -> data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds
  -> workflow/make_ou_tracking_manifest.R
  -> data/ou_tracking_manifest_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv
  -> workflow/fit_ou_tracking_manifest_row.R or workflow/submit_ou_tracking_manifest_array.slurm
  -> data/ou_tracking_fits_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv
  -> workflow/precompute_ou_report_diagnostics.R
  -> data/ou_report_diagnostics_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds
  -> ou_tracking_fit_report.Rmd
```

The rendered primary report is `ou_tracking_fit_report.html`.

## Kept Inputs and Generated Artifacts

- `Gemcitabine_PlateMap_20240111.xlsx`: source treatment map.
- `platemap_conditions.csv`: parsed plate conditions retained for review.
- `cpsam_full_stacks/`: user-generated CPSAM mask stack outputs needed to
  regenerate CPSAM/nucleus links locally. Codex should not run these GPU jobs.
- `cpsam_nucleus_tracking_links/`: CPSAM/nucleus link outputs. The
  `yellow_track_segments_len3/` subfolder is the direct upstream input for the
  current primary analysis.
- `cpsam_frame_sample/`: compact CPSAM visual QC artifact set.
- `overlay_checks/ou_fit_sanity/`: retained OU fit visual sanity overlays.
- `data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds`:
  primary staged track object.
- `data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3_summary.csv`:
  compact summary of the primary staged track object.
- `data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3_track_summary.csv`:
  per-track summary of the primary staged track object.
- `data/ou_tracking_manifest_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv`:
  current OU fit manifest.
- `data/ou_tracking_fits_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv`:
  current multistart OU fit results.
- `data/ou_report_diagnostics_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds`:
  compact precomputed report diagnostics. This keeps report rendering fast.

## Regeneration Commands

Run R commands through the repository wrapper from the repository root.

Build the staged yellow reconstructed track object:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/build_yellow_reconstructed_ou_tracks.R
```

Build the OU manifest:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/make_ou_tracking_manifest.R
```

Fit one OU manifest row/start locally:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/fit_ou_tracking_manifest_row.R \
  --manifest=analyses/K00_GemcitabineExposure_033023/data/ou_tracking_manifest_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv \
  --file_index=1 \
  --start_id=1
```

Submit the OU multistart array:

```bash
sbatch analyses/K00_GemcitabineExposure_033023/workflow/submit_ou_tracking_manifest_array.slurm
```

Precompute report diagnostics. This uses 16 local R workers by default:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/precompute_ou_report_diagnostics.R
```

Render the primary report:

```bash
scripts/agentRrunner.sh analyses/K00_GemcitabineExposure_033023/workflow/render_ou_tracking_fit_report.R
```

## CellposeSAM and Link Generation

CellposeSAM execution is GPU-backed and should be run by the user, not by Codex.
The full-stack launcher is retained for user execution:

```bash
sbatch analyses/K00_GemcitabineExposure_033023/workflow/submit_cpsam_full_stacks_array.slurm
```

After CPSAM masks exist, regenerate CPSAM/nucleus links with:

```bash
sbatch analyses/K00_GemcitabineExposure_033023/workflow/submit_cpsam_nucleus_tracking_links_array.slurm
```

The link workflow excludes CSV frame `0` because the registered TIFF stacks
contain a duplicated first frame. Clean yellow track segments are defined as
consecutive 3-frame segments where each frame has exactly one nucleus track
point inside a non-touching CPSAM object.

## Nearest Object Distance

Current predictive checks in `dev/predictive_checks/track_censoring.R` show
that track censoring can have a large effect on MSD, step-autocorrelation, and
related migration summaries. Simulated observation boxes suggest that censoring
can depress both MSD and ACF relative to the raw OU process. Additional
proactive censoring based on nearest-object distance has not resolved this
bias, so the right modelling or filtering strategy remains unsettled.

The nearest-object-distance workpackage adds a stricter isolation covariate for
future track subsetting: not only whether a CPSAM object is currently touching
another object, but how far it is from the nearest other object and therefore
how likely it is to become overlap-censored at the next timepoint. This remains
diagnostic rather than a validated correction for censoring-biased migration
summaries.

The reusable implementation lives in `src/image_datamining/mask_distances.py`.
It uses one exact Euclidean distance transform per 2D mask frame to construct a
labeled Voronoi map, then scans neighboring Voronoi territories to recover each
object's nearest other object and the witness mask-pixel pair. Distances are
reported as mask-pixel-center to mask-pixel-center distances; the companion
`nearest_empty_gap_px` column is `max(distance - 1, 0)` for an empty-pixel gap
interpretation.

Generate per-object nearest-distance tables for all CPSAM full-stack masks with:

```bash
sbatch analyses/K00_GemcitabineExposure_033023/workflow/submit_cpsam_nearest_distances_array.slurm
```

Outputs are written under:

```text
cpsam_full_stacks/nearest_distances/
  site_manifest.tsv
  object_tables/*_nearest_object_distances.tsv
  summaries/*_nearest_object_distances_summary.json
```

Each object-table row is keyed by `site_id`, `frame`, and `cpsam_label` and
includes `nearest_cpsam_label`, `nearest_mask_center_distance_px`,
`nearest_empty_gap_px`, and the two endpoint coordinates defining the nearest
object vector. These tables can be joined into the CPSAM/nucleus link or
yellow-track staging workflow for censoring diagnostics and sensitivity checks,
but nearest-distance filtering should not yet be treated as a validated
correction before fitting OU/MSD summaries.

## Literature Positioning

The local migration literature review in
`literature_reviews/migration_literature_review.txt` supports the current
OU/PRW framing as a reasonable first quantitative language for cell tracking
data. Gail and Boone, Stokes et al., and Martens et al. establish the standard
use of persistent-random-walk or OU-style summaries for separating movement
magnitude, persistence time, and long-time diffusivity
(`literature/Gail_1970.txt`, `literature/Stokes_1991.txt`,
`literature/Martens_2006.txt`). Our current workflow follows that tradition,
but fits an OU velocity model directly to trajectories with observation noise
rather than only fitting a PRW equation to an MSD curve.

The review also argues against treating a good MSD fit as sufficient model
validation. Selmeczi et al. and Wu et al. show that velocity distributions,
velocity autocorrelation, conditional turning behavior, anisotropy, and
cell-to-cell heterogeneity can reveal failures of a simple OU/PRW model even
when MSD is well described (`literature/Selmeczi_2005.txt`,
`literature/Wu_2014.txt`). The promoted report therefore uses fitted OU
parameters as compact phenomenological summaries, then checks them against
observed-vs-simulated tracks, MSD, one-frame step lengths, directional
autocorrelation, conditional turning, track-length distributions, and
minimum-track-length sensitivity.

The experimental-design lesson from the review is equally important here:
classical motility studies often used low density, explicit interaction
exclusions, stable imaging intervals, adequate field of view, and deliberately
chosen trajectory inclusion rules. This K00 analysis is more opportunistic, so
the staged yellow-track inclusion criteria, overlap censoring diagnostics,
nearest-object-distance work, and IPCW stress checks should be treated as part
of the biological interpretation rather than only as preprocessing details.
Tao et al. provide a useful precedent for using persistence-style summaries in
complex imaging data while interpreting them alongside other biological
evidence (`literature/Tao_2019.txt`).

Overall, the current promoted workflow is best described as a modern OU/PRW
phenotyping workflow for accepted yellow trajectories, with stronger diagnostic
and censoring-awareness layers than a basic MSD-fit pipeline. Its main
remaining gaps relative to the review are hierarchical modeling of well, site,
and cell-level heterogeneity, plus deferred alive/dead and proliferation
adjustment.

## Interpretation Notes

Current data support greater 4N dispersal in the staged yellow-track input.
The OU fits suggest increased persistence, but tau is the most vulnerable
parameter because it depends on short tracks, filtering choices, and unmodeled
time-varying cell state. Step-length and MSD-like dispersal summaries should be
treated as more stable than tau alone.

Alive/dead classification and proliferation adjustment remain deferred until
staged joinable inputs are supplied.
