#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
  library(tibble)
})

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
source(file.path(dirname(dirname(script_file)), "R/censoring_audit_utils.R"))

usage <- paste0(
  "Usage: combine_ou_censoring_conditions.R [options]\n\n",
  "Options:\n",
  "  --analysis_dir=/path/to/analyses/K00_GemcitabineExposure_033023\n",
  "  --models_rds=/path/to/censoring_models.rds\n",
  "  --condition_dir=/path/to/condition_artifacts\n",
  "  --out_rds=/path/to/ou_censoring_simulation.rds\n",
  "  --max_lag=10\n",
  "  --parallel_cores=0  Use 0 for SLURM_CPUS_PER_TASK or parallel::detectCores().\n"
)

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
audit_dir <- audit_dir_from_script()
analysis_dir <- normalizePath(args$analysis_dir %||% analysis_dir_from_audit_dir(audit_dir), mustWork = TRUE)
models_rds <- normalizePath(args$models_rds %||% file.path(analysis_dir, "data/yellow_censoring_models.rds"), mustWork = TRUE)
condition_dir <- normalizePath(
  args$condition_dir %||% file.path(analysis_dir, "data/ou_censoring_simulation_conditions"),
  mustWork = TRUE
)
out_rds <- normalizePath(
  args$out_rds %||% file.path(analysis_dir, "data/ou_censoring_simulation.rds"),
  mustWork = FALSE
)
max_lag <- as.integer(args$max_lag %||% "10")
parallel_cores <- as.integer(args$parallel_cores %||% "0")
if (is.na(parallel_cores) || parallel_cores < 1L) {
  slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA_character_)))
  parallel_cores <- if (is.finite(slurm_cpus) && slurm_cpus > 0L) {
    slurm_cpus
  } else {
    max(1L, parallel::detectCores(logical = FALSE))
  }
}
track_cols <- c("source", "ploidy", "Gemcitabine", "dose_label", "site_id", "migration_track_id", "frame", "track_step", "x", "y")

select_track_cols <- function(x) {
  if (is.null(x) || nrow(x) == 0L || !all(track_cols %in% names(x))) {
    return(tibble(
      source = character(),
      ploidy = character(),
      Gemcitabine = numeric(),
      dose_label = character(),
      site_id = character(),
      migration_track_id = character(),
      frame = integer(),
      track_step = integer(),
      x = numeric(),
      y = numeric()
    ))
  }
  x |> select(all_of(track_cols))
}

prepare_observed_tracks_light <- function(analysis_dir, conditions) {
  tracks_rds <- file.path(
    analysis_dir,
    "data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"
  )
  tracks_raw <- readRDS(tracks_rds)
  tracks_raw |>
    transmute(
      source = "Observed",
      ploidy = as.character(.data$ploidy),
      Gemcitabine = as.numeric(.data$Gemcitabine),
      dose_label = order_dose_label(paste0(.data$Gemcitabine, " nM"), .data$Gemcitabine),
      site_id = as.character(.data$site_id),
      migration_track_id = as.character(.data$migration_track_id),
      frame = as.integer(.data$frame),
      x = as.numeric(.data$nucleus_x),
      y = as.numeric(.data$nucleus_y)
    ) |>
    semi_join(conditions, by = c("ploidy", "Gemcitabine")) |>
    arrange(.data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$site_id, .data$migration_track_id) |>
    mutate(track_step = row_number() - 1L) |>
    ungroup() |>
    select(all_of(track_cols))
}

split_track_groups <- function(tracks) {
  tracks |>
    mutate(.split_key = paste(.data$source, .data$ploidy, .data$Gemcitabine, sep = "\r")) |>
    group_split(.data$.split_key, .keep = FALSE)
}

parallel_bind <- function(groups, fun, cores) {
  if (length(groups) == 0L) {
    return(tibble())
  }
  bind_rows(parallel::mclapply(
    groups,
    fun,
    mc.cores = min(cores, length(groups)),
    mc.preschedule = FALSE
  ))
}

condition_files <- sort(list.files(
  condition_dir,
  pattern = "^ou_censoring_simulation_.*_seed[0-9]+[.]rds$",
  full.names = TRUE
))
if (length(condition_files) == 0L) {
  stop("No condition simulation artifacts found in ", condition_dir, call. = FALSE)
}

message("Loading ", length(condition_files), " condition simulation artifact(s)")
condition_artifacts <- lapply(condition_files, readRDS)
conditions <- bind_rows(lapply(condition_artifacts, `[[`, "condition")) |>
  distinct(.data$ploidy, .data$Gemcitabine)
fits <- bind_rows(lapply(condition_artifacts, `[[`, "fit"))

message("Preparing observed comparison tracks without object-context joins")
observed_tracks <- prepare_observed_tracks_light(analysis_dir, conditions)

simulated_uncensored <- bind_rows(lapply(condition_artifacts, `[[`, "simulated_uncensored"))
simulated_context_censored <- bind_rows(lapply(condition_artifacts, `[[`, "simulated_context_censored"))
simulated_full_censored <- bind_rows(lapply(condition_artifacts, `[[`, "simulated_full_censored"))

comparison_tracks <- bind_rows(
  observed_tracks,
  select_track_cols(simulated_uncensored),
  select_track_cols(simulated_context_censored),
  select_track_cols(simulated_full_censored)
)

track_counts <- comparison_tracks |>
  group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label) |>
  summarize(
    n_tracks = n_distinct(.data$migration_track_id),
    n_points = n(),
    .groups = "drop"
  )

message("Computing comparison metrics with ", parallel_cores, " core(s)")
track_groups <- split_track_groups(comparison_tracks)
steps <- parallel_bind(track_groups, compute_steps, parallel_cores)
metrics <- list(
  msd = parallel_bind(track_groups, function(x) summarize_msd(x, max_lag = max_lag), parallel_cores),
  velocity_autocorrelation = parallel_bind(
    split_track_groups(steps),
    function(x) summarize_velocity_autocorrelation(x, max_lag = max_lag),
    parallel_cores
  ),
  track_lengths = parallel_bind(track_groups, summarize_track_lengths, parallel_cores),
  speeds = steps |>
    mutate(speed_capped = cap_quantile(.data$speed)),
  track_counts = track_counts,
  batch_summary = bind_rows(lapply(condition_artifacts, `[[`, "batch_summary")),
  model_check_counts = bind_rows(lapply(condition_artifacts, `[[`, "model_check_counts"))
)

censoring_artifact <- readRDS(models_rds)
artifact <- list(
  simulated_tracks = NULL,
  metrics = metrics,
  fits = fits,
  conditions = conditions,
  condition_files = condition_files,
  max_lag = max_lag,
  parallel_cores = parallel_cores,
  censoring_model_comparison = censoring_artifact$model_comparison,
  created_at = Sys.time(),
  seeds = sort(unique(vapply(condition_artifacts, function(x) x$seed, numeric(1))))
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(artifact, out_rds)
message("Saved combined simulation artifact: ", out_rds)
