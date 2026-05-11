#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
  library(tibble)
})

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
source(file.path(dirname(dirname(script_file)), "R/censoring_audit_utils.R"))

usage <- paste0(
  "Usage: 02_simulate_ou_censoring_condition.R [options]\n\n",
  "Options:\n",
  "  --analysis_dir=/path/to/analyses/K00_GemcitabineExposure_033023\n",
  "  --models_rds=/path/to/censoring_models.rds\n",
  "  --fits_csv=/path/to/ou_tracking_fits.csv\n",
  "  --out_dir=/path/to/condition_artifacts\n",
  "  --condition=2N:0\n",
  "  --file_index=1\n",
  "  --manifest=/path/to/condition_manifest.tsv\n",
  "  --batch_tracks=1000\n",
  "  --parallel_batches=0  Use 0 for SLURM_CPUS_PER_TASK or parallel::detectCores().\n",
  "  --uncensored_sample_tracks=1000\n",
  "  --include_context_model=FALSE\n",
  "  --min_segment_frames=3\n",
  "  --max_frame=39\n",
  "  --max_lag=10\n",
  "  --seed=20260509\n",
  "  --overwrite=FALSE\n"
)

safe_condition_label <- function(ploidy, dose) {
  dose_txt <- gsub("\\.", "p", format(dose, trim = TRUE, scientific = FALSE))
  paste(ploidy, paste0(dose_txt, "nM"), sep = "_")
}

simulate_condition_batch <- function(fit, starts, batch_id, max_frame) {
  bind_rows(lapply(seq_len(nrow(starts)), function(track_i) {
    start <- starts[track_i, , drop = FALSE]
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
        migration_track_id = paste("sim", fit$fit_label[[1]], batch_id, track_i, sep = "::"),
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
}

keep_min_consecutive_segments <- function(tracks, min_segment_frames = 3L) {
  if (nrow(tracks) == 0L) {
    return(tracks)
  }

  tracks |>
    arrange(.data$source, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$source, .data$site_id, .data$migration_track_id) |>
    mutate(
      segment_index = cumsum(row_number() == 1L | .data$frame != lag(.data$frame) + 1L),
      original_migration_track_id = as.character(.data$migration_track_id)
    ) |>
    group_by(.data$source, .data$site_id, .data$migration_track_id, .data$segment_index) |>
    mutate(
      segment_n_frames = n_distinct(.data$frame),
      segment_track_step = row_number() - 1L
    ) |>
    ungroup() |>
    filter(.data$segment_n_frames >= min_segment_frames) |>
    mutate(
      migration_track_id = if_else(
        .data$segment_index == 1L,
        .data$original_migration_track_id,
        paste(.data$original_migration_track_id, "segment", .data$segment_index, sep = "::")
      ),
      track_step = .data$segment_track_step
    ) |>
    select(-all_of(c("segment_index", "original_migration_track_id", "segment_n_frames", "segment_track_step")))
}

attach_condition_context <- function(tracks, context_pool) {
  if (nrow(tracks) == 0L) {
    return(bind_cols(tracks, context_pool[integer(), c("nearest_neighbor_distance", "confluency", "confluency_percent")]))
  }
  context_idx_by_frame <- split(seq_len(nrow(context_pool)), as.character(context_pool$frame))
  sampled_idx <- integer(nrow(tracks))
  fallback_idx <- seq_len(nrow(context_pool))
  row_groups <- split(seq_len(nrow(tracks)), as.character(tracks$frame))
  for (frame_key in names(row_groups)) {
    rows <- row_groups[[frame_key]]
    choices <- context_idx_by_frame[[frame_key]]
    if (is.null(choices) || length(choices) == 0L) {
      choices <- fallback_idx
    }
    sampled_idx[rows] <- sample(choices, length(rows), replace = TRUE)
  }
  bind_cols(
    tracks,
    context_pool[sampled_idx, c("nearest_neighbor_distance", "confluency", "confluency_percent")]
  )
}

append_track_sample <- function(existing, incoming, max_tracks) {
  if (max_tracks <= 0L || nrow(incoming) == 0L) {
    return(existing)
  }
  combined <- bind_rows(existing, incoming)
  keep_ids <- unique(combined$migration_track_id)
  if (length(keep_ids) > max_tracks) {
    keep_ids <- keep_ids[seq_len(max_tracks)]
  }
  combined |>
    filter(.data$migration_track_id %in% keep_ids)
}

simulate_censored_batch <- function(
  batch_id, fit, start_batches, context_pool, censoring_artifact,
  max_frame, include_context_model, min_segment_frames, seed
) {
  set.seed(seed + 100000L + batch_id)
  batch_starts <- start_batches[[batch_id]]
  proposed <- simulate_condition_batch(fit, batch_starts, batch_id, max_frame) |>
    mutate(
      dose_label = order_dose_label(.data$dose_label, .data$Gemcitabine),
      ploidy = factor(.data$ploidy, levels = sort(unique(censoring_artifact$yellow_start_pool$ploidy)))
    )

  proposed <- attach_condition_context(proposed, context_pool) |>
    add_track_speeds() |>
    filter(
      is.finite(.data$model_speed),
      is.finite(.data$distance_to_edge),
      is.finite(.data$nearest_neighbor_distance),
      is.finite(.data$confluency)
    )

  context_batch <- tibble()
  if (include_context_model) {
    context_batch <- censor_tracks_with_model(
      proposed,
      model = censoring_artifact$models$context,
      model_name = "context",
      scaler = censoring_artifact$scaler,
      seed = seed + 1000L + batch_id
    ) |>
      keep_min_consecutive_segments(min_segment_frames = min_segment_frames)
  }

  full_batch <- censor_tracks_with_model(
    proposed,
    model = censoring_artifact$models$full,
    model_name = "full",
    scaler = censoring_artifact$scaler,
    seed = seed + 2000L + batch_id
  ) |>
    keep_min_consecutive_segments(min_segment_frames = min_segment_frames)

  list(
    batch_id = batch_id,
    proposed = proposed,
    context_censored = context_batch,
    full_censored = full_batch,
    full_steps = compute_steps(full_batch)
  )
}

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
models_rds <- normalizePath(args$models_rds %||% file.path(audit_dir, "artifacts/censoring_models.rds"), mustWork = TRUE)
out_dir <- normalizePath(
  args$out_dir %||% file.path(audit_dir, "artifacts/ou_censoring_simulation_conditions"),
  mustWork = FALSE
)
condition <- args$condition %||% NULL
file_index <- as.integer(args$file_index %||% Sys.getenv("FILE_INDEX", Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_)))
manifest <- args$manifest %||% NULL
batch_tracks <- as.integer(args$batch_tracks %||% "1000")
parallel_batches <- as.integer(args$parallel_batches %||% "0")
uncensored_sample_tracks <- as.integer(args$uncensored_sample_tracks %||% "1000")
include_context_model <- as_flag(args$include_context_model, default = FALSE)
min_segment_frames <- as.integer(args$min_segment_frames %||% "3")
max_frame <- as.integer(args$max_frame %||% "39")
max_lag <- as.integer(args$max_lag %||% "10")
seed <- as.integer(args$seed %||% "20260509")
overwrite <- as_flag(args$overwrite, default = FALSE)
seed_offset <- if (is.na(file_index)) 0L else file_index
if (is.na(parallel_batches) || parallel_batches < 1L) {
  slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA_character_)))
  parallel_batches <- if (is.finite(slurm_cpus) && slurm_cpus > 0L) {
    slurm_cpus
  } else {
    max(1L, parallel::detectCores(logical = FALSE))
  }
}

if (is.null(condition) && !is.null(manifest)) {
  manifest_dt <- read.delim(manifest, stringsAsFactors = FALSE)
  if (is.na(file_index) || file_index < 1L || file_index > nrow(manifest_dt)) {
    message("file_index is outside manifest rows; exiting without work.")
    quit(status = 0)
  }
  condition <- manifest_dt$condition[[file_index]]
}
if (is.null(condition) || !nzchar(condition)) {
  stop("Provide --condition=ploidy:dose or --manifest plus --file_index.", call. = FALSE)
}
condition_df <- parse_condition_tokens(condition)
if (nrow(condition_df) != 1L) {
  stop("Condition-level simulation requires exactly one condition.", call. = FALSE)
}

set.seed(seed + seed_offset)
message("Loading censoring models: ", models_rds)
censoring_artifact <- readRDS(models_rds)
if (is.null(censoring_artifact$yellow_start_pool)) {
  censoring_artifact$yellow_start_pool <- censoring_artifact$track_start_pool
}

fits <- best_ou_fits(fits_csv, conditions = condition_df)
if (nrow(fits) != 1L) {
  stop("Expected exactly one OU fit for condition ", condition, "; found ", nrow(fits), call. = FALSE)
}
fit <- fits[1, , drop = FALSE]
condition_label <- safe_condition_label(fit$ploidy[[1]], fit$Gemcitabine[[1]])
out_rds <- file.path(out_dir, paste0("ou_censoring_simulation_", condition_label, "_seed", seed, ".rds"))
if (file.exists(out_rds) && !overwrite) {
  message("Output exists; rerun with --overwrite=TRUE: ", out_rds)
  quit(status = 0)
}

start_pool <- censoring_artifact$yellow_start_pool |>
  mutate(ploidy = as.character(.data$ploidy), Gemcitabine = as.numeric(.data$Gemcitabine)) |>
  filter(
    .data$ploidy == fit$ploidy[[1]],
    .data$Gemcitabine == fit$Gemcitabine[[1]],
    .data$first_frame < pmin(.data$site_last_frame, max_frame)
  )
if (nrow(start_pool) == 0L) {
  stop("No start-pool rows for condition ", condition, call. = FALSE)
}
n_sim_starts <- nrow(start_pool)
assigned_starts <- start_pool[sample.int(n_sim_starts), , drop = FALSE]
assigned_starts$sim_start_index <- seq_len(n_sim_starts)
start_batches <- split(
  assigned_starts,
  ceiling(seq_len(n_sim_starts) / batch_tracks)
)
n_batches <- length(start_batches)

context_pool <- censoring_artifact$context_pool |>
  mutate(ploidy = as.character(.data$ploidy), Gemcitabine = as.numeric(.data$Gemcitabine)) |>
  filter(.data$ploidy == fit$ploidy[[1]], .data$Gemcitabine == fit$Gemcitabine[[1]])
if (nrow(context_pool) == 0L) {
  stop("No context-pool rows for condition ", condition, call. = FALSE)
}

uncensored_sample <- tibble()
context_censored <- tibble()
full_censored <- tibble()
batch_summary <- vector("list", n_batches)
cumulative_full_steps <- 0L

message(
  "Simulating condition ", condition, " from ", n_sim_starts,
  " observed yellow starts in batches of up to ", batch_tracks,
  "; parallel_batches=", parallel_batches
)
batch_start <- 1L
while (batch_start <= n_batches) {
  batch_ids <- batch_start:min(n_batches, batch_start + parallel_batches - 1L)
  message("Starting batch wave: ", paste(batch_ids, collapse = ","))
  batch_results <- parallel::mclapply(
    batch_ids,
    simulate_censored_batch,
    fit = fit,
    start_batches = start_batches,
    context_pool = context_pool,
    censoring_artifact = censoring_artifact,
    max_frame = max_frame,
    include_context_model = include_context_model,
    min_segment_frames = min_segment_frames,
    seed = seed,
    mc.cores = min(parallel_batches, length(batch_ids)),
    mc.preschedule = FALSE
  )

  for (result in batch_results) {
    proposed_tracks <- n_distinct(result$proposed$migration_track_id)
    uncensored_sample <- append_track_sample(uncensored_sample, result$proposed, uncensored_sample_tracks)
    if (include_context_model) {
      context_censored <- bind_rows(context_censored, result$context_censored)
    }
    full_censored <- bind_rows(full_censored, result$full_censored)
    cumulative_full_steps <- cumulative_full_steps + nrow(result$full_steps)

    batch_summary[[result$batch_id]] <- tibble(
      condition = condition,
      batch_id = result$batch_id,
      proposed_tracks = proposed_tracks,
      proposed_rows = nrow(result$proposed),
      full_min3_tracks = n_distinct(result$full_censored$migration_track_id),
      full_min3_rows = nrow(result$full_censored),
      full_min3_steps = nrow(result$full_steps),
      cumulative_proposed_tracks = min(result$batch_id * batch_tracks, n_sim_starts),
      cumulative_full_min3_tracks = n_distinct(full_censored$migration_track_id),
      cumulative_full_min3_rows = nrow(full_censored),
      cumulative_full_min3_steps = cumulative_full_steps
    )
  }

  latest <- batch_summary[[max(batch_ids)]]
  message(
    "Completed through batch ", max(batch_ids), ": full min3 tracks=",
    latest$cumulative_full_min3_tracks, ", steps=", latest$cumulative_full_min3_steps
  )

  batch_start <- max(batch_ids) + 1L
}
batch_summary <- bind_rows(batch_summary)

observed_final_source <- censoring_artifact$final_track_summary
if (is.null(observed_final_source)) {
  observed_final_source <- censoring_artifact$track_start_pool
}
observed_final_tracks <- observed_final_source |>
  mutate(ploidy = as.character(.data$ploidy), Gemcitabine = as.numeric(.data$Gemcitabine)) |>
  filter(.data$ploidy == fit$ploidy[[1]], .data$Gemcitabine == fit$Gemcitabine[[1]])
observed_final_rows <- if ("track_length" %in% names(observed_final_tracks)) {
  sum(observed_final_tracks$track_length, na.rm = TRUE)
} else {
  NA_integer_
}
observed_final_steps <- if (is.finite(observed_final_rows)) {
  sum(pmax(observed_final_tracks$track_length - 1L, 0L), na.rm = TRUE)
} else {
  NA_integer_
}
model_check_counts <- tibble(
  condition = condition,
  observed_yellow_starts = n_sim_starts,
  simulated_proposed_starts = n_sim_starts,
  observed_final_min3_tracks = nrow(observed_final_tracks),
  simulated_final_min3_tracks = n_distinct(full_censored$migration_track_id),
  observed_final_rows = observed_final_rows,
  simulated_final_rows = nrow(full_censored),
  observed_final_steps = observed_final_steps,
  simulated_final_steps = nrow(compute_steps(full_censored))
)

artifact <- list(
  condition = condition_df,
  condition_label = condition_label,
  simulated_uncensored = uncensored_sample,
  simulated_context_censored = context_censored,
  simulated_full_censored = full_censored,
  batch_summary = batch_summary,
  model_check_counts = model_check_counts,
  fit = fit,
  n_sim_starts = n_sim_starts,
  batch_tracks = batch_tracks,
  parallel_batches = parallel_batches,
  include_context_model = include_context_model,
  uncensored_sample_tracks = uncensored_sample_tracks,
  min_segment_frames = min_segment_frames,
  max_frame = max_frame,
  max_lag = max_lag,
  seed = seed,
  created_at = Sys.time()
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(artifact, out_rds)
message("Saved condition simulation artifact: ", out_rds)
