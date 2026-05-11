#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
source(file.path(dirname(dirname(script_file)), "R/censoring_audit_utils.R"))

usage <- paste0(
  "Usage: simulate_ou_censoring.R [options]\n\n",
  "Options:\n",
  "  --analysis_dir=/path/to/analyses/K00_GemcitabineExposure_033023\n",
  "  --models_rds=/path/to/censoring_models.rds\n",
  "  --fits_csv=/path/to/ou_tracking_fits.csv\n",
  "  --out_rds=/path/to/ou_censoring_simulation.rds\n",
  "  --conditions=2N:0,2N:12.5,4N:0,4N:12.5\n",
  "  --n_tracks_per_condition=500\n",
  "  --max_frame=39\n",
  "  --max_lag=10\n",
  "  --seed=20260509\n"
)

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
audit_dir <- audit_dir_from_script()
analysis_dir <- normalizePath(args$analysis_dir %||% analysis_dir_from_audit_dir(audit_dir), mustWork = TRUE)
fits_csv <- normalizePath(
  args$fits_csv %||% file.path(
    analysis_dir,
    "data/ou_tracking_fits_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv"
  ),
  mustWork = TRUE
)
models_rds <- normalizePath(args$models_rds %||% file.path(analysis_dir, "data/yellow_censoring_models.rds"), mustWork = TRUE)
out_rds <- normalizePath(
  args$out_rds %||% file.path(analysis_dir, "data/ou_censoring_simulation.rds"),
  mustWork = FALSE
)
conditions <- parse_condition_tokens(args$conditions %||% "2N:0,2N:12.5,4N:0,4N:12.5")
n_tracks_per_condition <- as.integer(args$n_tracks_per_condition %||% "500")
max_frame <- as.integer(args$max_frame %||% "39")
max_lag <- as.integer(args$max_lag %||% "10")
seed <- as.integer(args$seed %||% "20260509")

set.seed(seed)
message("Loading censoring models: ", models_rds)
censoring_artifact <- readRDS(models_rds)

message("Preparing observed comparison data")
observed <- prepare_observed_at_risk(analysis_dir)
fits <- best_ou_fits(fits_csv, conditions = conditions)
if (nrow(fits) == 0L) {
  stop("No OU fits matched requested conditions.", call. = FALSE)
}

start_pool <- censoring_artifact$track_start_pool |>
  mutate(ploidy = as.character(.data$ploidy), Gemcitabine = as.numeric(.data$Gemcitabine)) |>
  semi_join(fits |> select("ploidy", "Gemcitabine"), by = c("ploidy", "Gemcitabine")) |>
  filter(.data$first_frame < pmin(.data$site_last_frame, max_frame))

context_pool <- censoring_artifact$context_pool |>
  mutate(
    ploidy = as.character(.data$ploidy),
    Gemcitabine = as.numeric(.data$Gemcitabine),
    context_key = paste(.data$ploidy, .data$Gemcitabine, .data$frame, sep = "::"),
    condition_key = paste(.data$ploidy, .data$Gemcitabine, sep = "::"),
    frame_key = as.character(.data$frame)
  )
context_by_key <- split(seq_len(nrow(context_pool)), context_pool$context_key)
context_by_condition <- split(seq_len(nrow(context_pool)), context_pool$condition_key)
context_by_frame <- split(seq_len(nrow(context_pool)), context_pool$frame_key)

sample_context_row <- function(ploidy, dose, frame) {
  context_key <- paste(ploidy, dose, frame, sep = "::")
  condition_key <- paste(ploidy, dose, sep = "::")
  frame_key <- as.character(frame)
  idx <- context_by_key[[context_key]]
  if (is.null(idx)) {
    idx <- context_by_condition[[condition_key]]
  }
  if (is.null(idx)) {
    idx <- context_by_frame[[frame_key]]
  }
  if (is.null(idx)) {
    idx <- seq_len(nrow(context_pool))
  }
  context_pool[sample(idx, 1L), c("nearest_neighbor_distance", "confluency", "confluency_percent")]
}

message("Simulating OU tracks for ", nrow(fits), " conditions")
simulated_uncensored <- bind_rows(lapply(seq_len(nrow(fits)), function(fit_i) {
  fit <- fits[fit_i, , drop = FALSE]
  starts <- start_pool |>
    filter(.data$ploidy == fit$ploidy[[1]], .data$Gemcitabine == fit$Gemcitabine[[1]])
  if (nrow(starts) == 0L) {
    warning("No observed starts for ", fit$ploidy[[1]], " ", fit$Gemcitabine[[1]], " nM", call. = FALSE)
    return(NULL)
  }
  sampled_starts <- starts[sample.int(nrow(starts), n_tracks_per_condition, replace = TRUE), , drop = FALSE]

  bind_rows(lapply(seq_len(nrow(sampled_starts)), function(track_i) {
    start <- sampled_starts[track_i, , drop = FALSE]
    last_frame <- min(start$site_last_frame[[1]], max_frame)
    sim <- simulate_ou_track(
      tau = fit$tau[[1]],
      velocity_scale = fit$velocity_scale[[1]],
      obs_noise = fit$obs_noise[[1]],
      first_frame = start$first_frame[[1]],
      last_frame = last_frame,
      dt = fit$frame_interval[[1]]
    )
    sim |>
      mutate(
        source = "OU uncensored",
        ploidy = fit$ploidy[[1]],
        Gemcitabine = fit$Gemcitabine[[1]],
        dose_label = paste0(.data$Gemcitabine, " nM"),
        site_id = start$site_id[[1]],
        well = start$well[[1]],
        position = start$position[[1]],
        migration_track_id = paste("sim", fit$fit_label[[1]], track_i, sep = "::"),
        width = start$width[[1]],
        height = start$height[[1]],
        site_last_frame = start$site_last_frame[[1]],
        x = start$start_x[[1]] + .data$rel_x - first(.data$rel_x),
        y = start$start_y[[1]] + .data$rel_y - first(.data$rel_y),
        distance_to_edge = pmin(.data$x, .data$width - .data$x, .data$y, .data$height - .data$y, na.rm = FALSE),
        dose_log10 = log10(.data$Gemcitabine + 1),
        frame_window = time_window_label(.data$frame)
      )
  }))
}))

simulated_uncensored <- simulated_uncensored |>
  mutate(
    dose_label = order_dose_label(.data$dose_label, .data$Gemcitabine),
    ploidy = factor(.data$ploidy, levels = levels(observed$at_risk$ploidy))
  )

message("Attaching empirical whole-frame confluency and object-gap contexts")
context_rows <- bind_rows(lapply(seq_len(nrow(simulated_uncensored)), function(i) {
  sample_context_row(
    simulated_uncensored$ploidy[[i]],
    simulated_uncensored$Gemcitabine[[i]],
    simulated_uncensored$frame[[i]]
  )
}))

simulated_uncensored <- bind_cols(simulated_uncensored, context_rows) |>
  add_track_speeds() |>
  filter(
    is.finite(.data$model_speed),
    is.finite(.data$distance_to_edge),
    is.finite(.data$nearest_neighbor_distance),
    is.finite(.data$confluency)
  )

message("Applying context-only censoring model")
simulated_context_censored <- censor_tracks_with_model(
  simulated_uncensored,
  model = censoring_artifact$models$context,
  model_name = "context",
  scaler = censoring_artifact$scaler,
  seed = seed + 1L
)

message("Applying full censoring model")
simulated_full_censored <- censor_tracks_with_model(
  simulated_uncensored,
  model = censoring_artifact$models$full,
  model_name = "full",
  scaler = censoring_artifact$scaler,
  seed = seed + 2L
)

observed_tracks <- observed$tracks |>
  mutate(
    source = "Observed",
    migration_track_id = as.character(.data$migration_track_id),
    dose_label = order_dose_label(paste0(.data$Gemcitabine, " nM"), .data$Gemcitabine),
    track_step = ave(.data$frame, .data$migration_track_id, FUN = function(x) rank(x, ties.method = "first") - 1L)
  ) |>
  semi_join(fits |> select("ploidy", "Gemcitabine"), by = c("ploidy", "Gemcitabine")) |>
  select(all_of(c(
    "source", "ploidy", "Gemcitabine", "dose_label", "site_id",
    "migration_track_id", "frame", "track_step", "x", "y"
  )))

comparison_tracks <- bind_rows(
  observed_tracks,
  simulated_uncensored |>
    select(all_of(c("source", "ploidy", "Gemcitabine", "dose_label", "site_id", "migration_track_id", "frame", "track_step", "x", "y"))),
  simulated_context_censored |>
    select(all_of(c("source", "ploidy", "Gemcitabine", "dose_label", "site_id", "migration_track_id", "frame", "track_step", "x", "y"))),
  simulated_full_censored |>
    select(all_of(c("source", "ploidy", "Gemcitabine", "dose_label", "site_id", "migration_track_id", "frame", "track_step", "x", "y")))
)

steps <- compute_steps(comparison_tracks)
metrics <- list(
  msd = summarize_msd(comparison_tracks, max_lag = max_lag),
  velocity_autocorrelation = summarize_velocity_autocorrelation(steps, max_lag = max_lag),
  track_lengths = summarize_track_lengths(comparison_tracks),
  speeds = steps |>
    mutate(speed_capped = cap_quantile(.data$speed))
)

artifact <- list(
  simulated_tracks = comparison_tracks,
  metrics = metrics,
  fits = fits,
  conditions = conditions,
  n_tracks_per_condition = n_tracks_per_condition,
  max_frame = max_frame,
  max_lag = max_lag,
  censoring_model_comparison = censoring_artifact$model_comparison,
  created_at = Sys.time(),
  seed = seed
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(artifact, out_rds)
message("Saved simulation artifact: ", out_rds)
