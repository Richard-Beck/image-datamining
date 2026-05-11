# Censoring Audit Simulation Workflow

This folder separates expensive censoring-model fitting and OU simulation from the
HTML report. Run all R scripts through the repository wrapper:

```bash
sbatch analyses/K00_GemcitabineExposure_033023/workflow/submit_yellow_censoring_sites_array.slurm
sbatch analyses/K00_GemcitabineExposure_033023/workflow/submit_yellow_censoring_model_fit.slurm
sbatch analyses/K00_GemcitabineExposure_033023/workflow/submit_ou_censoring_conditions_array.slurm
sbatch analyses/K00_GemcitabineExposure_033023/workflow/submit_ou_censoring_combine_report.slurm
```

Main artifacts:

- `staged/sites/*_yellow_censoring_at_risk.rds`: site-local raw nucleus-detection
  censoring tables. These are built in parallel from `trackpoint_links`, one
  site per Slurm task.
- `staged/yellow_censoring_at_risk.rds`: combined raw at-risk model input used by
  `01_fit_censoring_models.R` when present.
- `staged/yellow_censoring_stage_summary.csv`: raw-to-yellow stage counts for
  checking where detections are removed.
- `staged/yellow_censoring_reason_summary.csv`: next-frame censoring reason
  counts.
- `artifacts/censoring_models.rds`: context-only and full mixed logistic dropout models.
- `artifacts/ou_censoring_simulation_conditions/*.rds`: per-condition OU
  simulation artifacts. These are built in parallel by Slurm, one condition per
  array task, with each task using `parallel::mclapply()` across batches.
- `artifacts/ou_censoring_simulation.rds`: observed, OU-uncensored, context-censored, and full-censored tracks plus comparison metrics.
- `ou_censoring_simulation_report.html`: rendered comparison report.

The censoring model fit now prefers the combined raw at-risk table at
`staged/yellow_censoring_at_risk.rds`. If that file is absent, it falls back to
the older already-staged yellow-track input.

The condition simulation no longer uses target censored-track or censored-step
quotas. Each condition simulates one proposed OU start for each observed
point-level yellow start in that condition, then applies the min-3 consecutive
frame analysis-track rule after censoring. Use `BATCH_TRACKS` and
`PARALLEL_BATCHES` to tune within-task batching; the default Slurm runner uses
20 condition tasks and 8 CPUs per task.

The final combine/report step keeps the report artifact lightweight: it stores
comparison metrics, track counts, per-track lengths, and per-step speeds, but not
the full bound observed/simulated track table. Its Slurm runner requests 32 cores
and parallelizes summary calculations within R using `PARALLEL_CORES`. The
default report lag horizon is 5 frames for MSD and velocity autocorrelation.
