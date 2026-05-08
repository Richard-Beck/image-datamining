#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

usage <- paste0(
  "Usage: fit_ou_tracking_manifest_row.R --manifest=/path/to/manifest.csv --file_index=1 [options]\n\n",
  "Options:\n",
  "  --repo_root=/path/to/repo\n",
  "  --job_id=1\n",
  "  --n_starts=25\n",
  "  --start_id=1\n",
  "  --seed=17\n",
  "  --fail_on_error=FALSE\n",
  "  --log_file=/path/to/task.log\n",
  "  --quiet=TRUE\n"
)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
analysis_dir_guess <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(analysis_dir_guess, "R/k00_batch_utils.R"))
source(file.path(analysis_dir_guess, "R/ou_velocity_model.R"))

select_manifest_row <- function(manifest, args) {
  if (!is.null(args$job_id)) {
    job_id <- as.integer(args$job_id)
    row <- manifest |> filter(.data$job_id == job_id)
    if (nrow(row) == 0) {
      stop("job_id not found in manifest: ", job_id, call. = FALSE)
    }
    return(row[1, , drop = FALSE])
  }

  file_index <- as.integer(args$file_index %||% Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
  if (is.na(file_index) || file_index < 1) {
    stop("--file_index must be a positive integer", call. = FALSE)
  }
  if (file_index > nrow(manifest)) {
    message(sprintf("file_index %d is beyond %d manifest rows; exiting without work.", file_index, nrow(manifest)))
    quit(save = "no", status = 0)
  }
  manifest[file_index, , drop = FALSE]
}

as_optional_numeric <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    return(NA_real_)
  }
  as.numeric(x)
}

prepare_ou_tracks <- function(tracks, row) {
  required <- c(
    "frame", "migration_track_id", "Center_of_the_object_1", "Center_of_the_object_0",
    "ploidy", "Gemcitabine"
  )
  if (!all(required %in% names(tracks))) {
    stop("Input tracks are missing columns: ", paste(setdiff(required, names(tracks)), collapse = ", "), call. = FALSE)
  }

  min_frame <- as_optional_numeric(row$min_frame[[1]])
  max_frame <- as_optional_numeric(row$max_frame[[1]])
  ploidy_filter <- row$ploidy[[1]]
  gemcitabine_filter <- row$Gemcitabine[[1]]

  out <- tracks |>
    filter(
      .data$ploidy == ploidy_filter,
      .data$Gemcitabine == gemcitabine_filter
    )

  if (is.finite(min_frame)) {
    out <- out |> filter(.data$frame >= min_frame)
  }
  if (is.finite(max_frame)) {
    out <- out |> filter(.data$frame <= max_frame)
  }

  out |>
    transmute(
      split_track_id = .data$migration_track_id,
      frame = as.numeric(.data$frame),
      x = as.numeric(.data$Center_of_the_object_1),
      y = as.numeric(.data$Center_of_the_object_0)
    )
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
repo_root <- normalizePath(args$repo_root %||% normalizePath(file.path(analysis_dir_guess, "../.."), mustWork = TRUE), mustWork = TRUE)
analysis_dir <- file.path(repo_root, "analyses/K00_GemcitabineExposure_033023")
manifest_path <- normalizePath(
  args$manifest %||% file.path(analysis_dir, "data/ou_tracking_manifest_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv"),
  mustWork = TRUE
)
fail_on_error <- as_flag(args$fail_on_error, default = FALSE)
quiet <- as_flag(args$quiet, default = FALSE)
if (!is.null(args$log_file)) {
  log_file <- normalizePath(args$log_file, mustWork = FALSE)
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  message("Log file requested: ", log_file)
  message("Task stdout/stderr should already be redirected there by the Slurm wrapper.")
}

log_msg <- function(...) {
  if (!quiet) {
    message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), paste0(..., collapse = "")))
  }
}

manifest <- read_csv(manifest_path, show_col_types = FALSE)
row <- select_manifest_row(manifest, args)
input_rds <- normalizePath(row$input_rds[[1]], mustWork = TRUE)
out_csv <- normalizePath(row$out_csv[[1]], mustWork = FALSE)
frame_interval <- as.numeric(row$frame_interval[[1]])
min_segment_frames <- as.integer(row$min_segment_frames[[1]])
max_tracks <- as.integer(row$max_tracks[[1]])
max_segments <- as.integer(row$max_segments[[1]])
n_starts <- as.integer(args$n_starts %||% if ("n_starts" %in% names(row)) row$n_starts[[1]] else "25")
start_id <- if (!is.null(args$start_id)) as.integer(args$start_id) else NA_integer_
seed <- as.integer(args$seed %||% if ("seed" %in% names(row)) row$seed[[1]] else "17")

if (is.na(frame_interval) || frame_interval <= 0) {
  stop("frame_interval must be positive for job_id=", row$job_id[[1]], call. = FALSE)
}
if (is.na(min_segment_frames) || min_segment_frames < 2) {
  stop("min_segment_frames must be >= 2 for job_id=", row$job_id[[1]], call. = FALSE)
}
if (is.na(n_starts) || n_starts < 1) {
  stop("n_starts must be a positive integer for job_id=", row$job_id[[1]], call. = FALSE)
}
if (is.na(seed)) {
  stop("seed must be an integer for job_id=", row$job_id[[1]], call. = FALSE)
}
if (!is.na(start_id) && (start_id < 1 || start_id > n_starts)) {
  stop("start_id must be between 1 and n_starts for job_id=", row$job_id[[1]], call. = FALSE)
}

log_msg("Manifest: ", manifest_path)
log_msg("Manifest rows: ", nrow(manifest))
log_msg("Selected job_id=", row$job_id[[1]], ", file_index=", args$file_index %||% Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
log_msg("Fit label: ", row$fit_label[[1]])
log_msg("Input RDS: ", input_rds)
log_msg("Output CSV: ", out_csv)
log_msg("Filters: ploidy=", row$ploidy[[1]], ", Gemcitabine=", row$Gemcitabine[[1]],
        ", min_frame=", row$min_frame[[1]], ", max_frame=", row$max_frame[[1]])
log_msg("OU fit starts: n_starts=", n_starts, ", start_id=", ifelse(is.na(start_id), "ALL", start_id), ", seed=", seed)

raw_tracks <- readRDS(input_rds)
log_msg("Input rows loaded: ", nrow(raw_tracks))
tracks <- prepare_ou_tracks(raw_tracks, row)
log_msg("Rows after manifest filters: ", nrow(tracks))
log_msg("Tracks after manifest filters: ", dplyr::n_distinct(tracks$split_track_id))
rm(raw_tracks)

if (max_tracks > 0 && nrow(tracks) > 0) {
  set.seed(seed)
  track_ids <- unique(tracks$split_track_id)
  keep_tracks <- sample(track_ids, size = min(length(track_ids), max_tracks))
  tracks <- tracks |> filter(.data$split_track_id %in% keep_tracks)
  log_msg("Rows after max_tracks sampling: ", nrow(tracks))
  log_msg("Tracks after max_tracks sampling: ", dplyr::n_distinct(tracks$split_track_id))
}

# Re-split after all manifest filters so time windows cannot leave short or non-contiguous segments.
tracks <- tracks |>
  arrange(.data$split_track_id, .data$frame) |>
  group_by(.data$split_track_id) |>
  mutate(
    starts_segment = row_number() == 1L | .data$frame != lag(.data$frame) + 1,
    post_filter_segment = cumsum(starts_segment),
    split_track_id = paste(.data$split_track_id, .data$post_filter_segment, sep = "::")
  ) |>
  ungroup() |>
  group_by(.data$split_track_id) |>
  filter(n_distinct(.data$frame) >= min_segment_frames) |>
  ungroup() |>
  select("split_track_id", "frame", "x", "y")
log_msg("Rows after post-filter contiguous segment split/min length: ", nrow(tracks))
log_msg("Segments after post-filter contiguous segment split/min length: ", dplyr::n_distinct(tracks$split_track_id))

segments <- split_ou_segments(tracks, frame_interval = frame_interval)
if (max_segments > 0 && length(segments) > max_segments) {
  set.seed(seed)
  segments <- segments[sample(seq_along(segments), max_segments)]
  attr(segments, "total_track_time") <- sum(vapply(segments, function(seg) sum(seg$dt), numeric(1)))
  log_msg("Segments after max_segments sampling: ", length(segments))
}

n_segments <- length(segments)
total_track_time <- attr(segments, "total_track_time") %||% 0
log_msg("OU segments prepared: ", n_segments)
log_msg("Total track time: ", total_track_time)
base_cols <- tibble(
  job_id = row$job_id[[1]],
  fit_label = row$fit_label[[1]],
  input_rds = input_rds,
  out_csv = out_csv,
  ploidy = row$ploidy[[1]],
  Gemcitabine = row$Gemcitabine[[1]],
  min_frame = row$min_frame[[1]],
  max_frame = row$max_frame[[1]],
  frame_interval = frame_interval,
  min_segment_frames = min_segment_frames,
  max_tracks = max_tracks,
  max_segments = max_segments
)

result <- if (n_segments == 0) {
  base_cols |>
    mutate(
      tau = NA_real_,
      velocity_scale = NA_real_,
      effective_diffusivity = NA_real_,
      obs_noise = NA_real_,
      log_likelihood = NA_real_,
      n_starts = n_starts,
      seed = seed,
      start_id = start_id,
      start_tau = NA_real_,
      start_velocity_scale = NA_real_,
      start_obs_noise = NA_real_,
      n_segments = 0L,
      total_track_time = 0,
      fit_status = "no_valid_segments"
    )
} else {
  tryCatch(
    {
      log_msg("Starting OU fit")
      fit <- fit_ou_velocity(
        segments,
        n_starts = n_starts,
        seed = seed,
        start_id = if (is.na(start_id)) NULL else start_id
      )
      log_msg("Finished OU fit rows: ", nrow(fit), "; statuses: ", paste(unique(fit$fit_status), collapse = ", "))
      bind_cols(base_cols, fit) |>
        mutate(
          n_segments = n_segments,
          total_track_time = total_track_time
        )
    },
    error = function(e) {
      if (fail_on_error) {
        stop(e)
      }
      base_cols |>
        mutate(
          tau = NA_real_,
          velocity_scale = NA_real_,
          effective_diffusivity = NA_real_,
          obs_noise = NA_real_,
          log_likelihood = NA_real_,
          n_starts = n_starts,
          seed = seed,
          start_id = start_id,
          start_tau = NA_real_,
          start_velocity_scale = NA_real_,
          start_obs_noise = NA_real_,
          n_segments = n_segments,
          total_track_time = total_track_time,
          fit_status = paste0("error: ", conditionMessage(e))
        )
    }
  )
}

append_delim_locked(result, out_csv, delim = ",", quiet = quiet)
log_msg("Appended result to ", out_csv)
log_msg("Output file exists: ", file.exists(out_csv), "; size bytes: ", if (file.exists(out_csv)) file.info(out_csv)$size else NA)
if (!quiet) {
  print(result)
}
