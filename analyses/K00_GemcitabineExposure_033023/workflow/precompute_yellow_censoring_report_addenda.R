#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
  library(readr)
  library(tibble)
  library(tidyr)
})

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
source(file.path(dirname(dirname(script_file)), "R/censoring_audit_utils.R"))

usage <- paste0(
  "Usage: precompute_yellow_censoring_report_addenda.R [options]\n\n",
  "Options:\n",
  "  --analysis_dir=/path/to/analyses/K00_GemcitabineExposure_033023\n",
  "  --simulation_rds=/path/to/ou_censoring_simulation.rds\n",
  "  --condition_dir=/path/to/ou_censoring_simulation_conditions\n",
  "  --tracks_rds=/path/to/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds\n",
  "  --out_rds=/path/to/yellow_censoring_simulation_report_addenda.rds\n",
  "  --selected_doses=0,25,200\n",
  "  --max_lag=10\n",
  "  --sample_per_group=5000\n",
  "  --mechanistic_tracks_per_condition=350\n",
  "  --mechanistic_step_size_quantile=0.99\n",
  "  --mechanistic_min_segment_frames=3\n",
  "  --skip_cosine=FALSE\n",
  "  --parallel_cores=0\n",
  "  --seed=20260511\n"
)

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
audit_dir <- audit_dir_from_script()
analysis_dir <- normalizePath(args$analysis_dir %||% analysis_dir_from_audit_dir(audit_dir), mustWork = TRUE)
simulation_rds <- normalizePath(args$simulation_rds %||% file.path(analysis_dir, "data/ou_censoring_simulation.rds"), mustWork = TRUE)
condition_dir <- normalizePath(args$condition_dir %||% file.path(analysis_dir, "data/ou_censoring_simulation_conditions"), mustWork = TRUE)
tracks_rds <- normalizePath(
  args$tracks_rds %||% file.path(analysis_dir, "data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"),
  mustWork = TRUE
)
out_rds <- normalizePath(
  args$out_rds %||% file.path(analysis_dir, "data/yellow_censoring_simulation_report_addenda.rds"),
  mustWork = FALSE
)
selected_doses <- as.numeric(strsplit(args$selected_doses %||% "0,25,200", ",", fixed = TRUE)[[1]])
max_lag <- as.integer(args$max_lag %||% "10")
sample_per_group <- as.integer(args$sample_per_group %||% "5000")
mechanistic_tracks_per_condition <- as.integer(args$mechanistic_tracks_per_condition %||% "350")
mechanistic_step_size_quantile <- as.numeric(args$mechanistic_step_size_quantile %||% "0.99")
mechanistic_min_segment_frames <- as.integer(args$mechanistic_min_segment_frames %||% "3")
skip_cosine <- as_flag(args$skip_cosine, default = FALSE)
parallel_cores <- as.integer(args$parallel_cores %||% "0")
seed <- as.integer(args$seed %||% "20260511")
if (is.na(parallel_cores) || parallel_cores < 1L) {
  slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA_character_)))
  parallel_cores <- if (is.finite(slurm_cpus) && slurm_cpus > 0L) slurm_cpus else max(1L, parallel::detectCores(logical = FALSE))
}

empty_tracks <- function() {
  tibble(
    source = character(), ploidy = character(), Gemcitabine = numeric(), dose_label = character(),
    site_id = character(), migration_track_id = character(), frame = integer(), track_step = integer(),
    x = numeric(), y = numeric()
  )
}

track_cols <- names(empty_tracks())
select_tracks <- function(x) {
  if (is.null(x) || nrow(x) == 0L || !all(track_cols %in% names(x))) {
    return(empty_tracks())
  }
  x |>
    select(all_of(track_cols)) |>
    mutate(
      source = as.character(.data$source),
      ploidy = as.character(.data$ploidy),
      Gemcitabine = as.numeric(.data$Gemcitabine),
      dose_label = as.character(.data$dose_label),
      site_id = as.character(.data$site_id),
      migration_track_id = as.character(.data$migration_track_id),
      frame = as.integer(.data$frame),
      track_step = as.integer(.data$track_step),
      x = as.numeric(.data$x),
      y = as.numeric(.data$y)
    )
}

prepare_observed_tracks <- function(tracks_rds, conditions) {
  readRDS(tracks_rds) |>
    transmute(
      source = "Observed",
      ploidy = as.character(.data$ploidy),
      Gemcitabine = as.numeric(.data$Gemcitabine),
      dose_label = paste0(.data$Gemcitabine, " nM"),
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
    mutate(dose_label = as.character(order_dose_label(.data$dose_label, .data$Gemcitabine))) |>
    select(all_of(track_cols))
}

compute_cosine_pairs <- function(tracks, max_lag) {
  if (nrow(tracks) == 0L) {
    return(tibble())
  }
  steps <- compute_steps(tracks)
  steps |>
    arrange(.data$source, .data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$site_id, .data$migration_track_id) |>
    group_modify(~ {
      df <- .x |> arrange(.data$frame)
      bind_rows(lapply(seq_len(max_lag), function(lag_frame) {
        j <- match(df$frame + lag_frame, df$frame)
        ok <- !is.na(j)
        if (!any(ok)) return(NULL)
        tibble(
          lag_frames = lag_frame,
          cosine_autocorrelation = (df$dx[ok] * df$dx[j[ok]] + df$dy[ok] * df$dy[j[ok]]) /
            pmax(df$speed[ok] * df$speed[j[ok]], 1e-9)
        )
      }))
    }) |>
    ungroup()
}

summarize_cosine_pairs <- function(x) {
  x |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$lag_frames) |>
    summarize(
      mean_cosine_autocorrelation = mean(.data$cosine_autocorrelation, na.rm = TRUE),
      sd_cosine_autocorrelation = sd(.data$cosine_autocorrelation, na.rm = TRUE),
      n_pairs = n(),
      se = .data$sd_cosine_autocorrelation / sqrt(.data$n_pairs),
      z = .data$mean_cosine_autocorrelation / .data$se,
      p_above_zero = pnorm(.data$z, lower.tail = FALSE),
      significantly_above_zero = is.finite(.data$p_above_zero) & .data$mean_cosine_autocorrelation > 0 & .data$p_above_zero < 0.05,
      .groups = "drop"
    )
}

sample_cosine_pairs <- function(x, sample_per_group, seed) {
  set.seed(seed)
  x |>
    filter(.data$Gemcitabine %in% selected_doses, .data$lag_frames == 1L) |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$lag_frames) |>
    group_modify(~ {
      n <- min(nrow(.x), sample_per_group)
      .x[sample.int(nrow(.x), n), , drop = FALSE]
    }) |>
    ungroup()
}

random_walk_comparator <- function(reference, sample_per_group, seed) {
  set.seed(seed)
  reference |>
    distinct(.data$ploidy, .data$Gemcitabine, .data$dose_label, .data$lag_frames) |>
    filter(.data$Gemcitabine %in% selected_doses, .data$lag_frames == 1L) |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$dose_label, .data$lag_frames) |>
    group_modify(~ tibble(
      source = "Random walk",
      cosine_autocorrelation = cos(runif(sample_per_group, -pi, pi))
    )) |>
    ungroup() |>
    select(source, ploidy, Gemcitabine, dose_label, lag_frames, cosine_autocorrelation)
}

split_mechanistic_segments <- function(tracks, min_segment_frames) {
  if (nrow(tracks) == 0L) {
    return(tracks |> mutate(segment_track_id = character(), segment_track_step = integer()))
  }
  tracks |>
    arrange(.data$source, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$source, .data$site_id, .data$migration_track_id) |>
    mutate(
      segment_index = cumsum(row_number() == 1L | .data$frame != lag(.data$frame) + 1L),
      original_track_id = as.character(.data$migration_track_id)
    ) |>
    group_by(.data$source, .data$site_id, .data$migration_track_id, .data$segment_index) |>
    mutate(
      segment_frames = n_distinct(.data$frame),
      segment_track_step = row_number() - 1L
    ) |>
    ungroup() |>
    filter(.data$segment_frames >= min_segment_frames) |>
    mutate(
      segment_track_id = paste(.data$original_track_id, "mechseg", .data$segment_index, sep = "::"),
      migration_track_id = .data$segment_track_id,
      track_step = .data$segment_track_step
    ) |>
    select(-"segment_index", -"original_track_id", -"segment_frames", -"segment_track_id", -"segment_track_step")
}

nearest_other_distance <- function(x, y) {
  n <- length(x)
  if (n < 2L) {
    return(rep(Inf, n))
  }
  dx <- outer(x, x, "-")
  dy <- outer(y, y, "-")
  distance <- sqrt(dx^2 + dy^2)
  diag(distance) <- Inf
  apply(distance, 1L, min)
}

mark_mechanistic_censoring <- function(tracks, disc_radius) {
  tracks |>
    mutate(
      row_id = row_number(),
      boundary_censored = .data$x - disc_radius < 0 |
        .data$x + disc_radius > .data$width |
        .data$y - disc_radius < 0 |
        .data$y + disc_radius > .data$height
    ) |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$site_id, .data$frame) |>
    mutate(
      nearest_current_object_distance = nearest_other_distance(.data$x, .data$y),
      overlap_censored = .data$nearest_current_object_distance < 2 * disc_radius
    ) |>
    ungroup() |>
    mutate(
      is_mechanically_censored = .data$boundary_censored | .data$overlap_censored,
      censor_reason = case_when(
        .data$boundary_censored & .data$overlap_censored ~ "boundary and overlap",
        .data$boundary_censored ~ "outside image boundary",
        .data$overlap_censored ~ "object overlap",
        TRUE ~ "kept"
      )
    )
}

prepare_mechanistic_steps <- function(tracks) {
  compute_steps(tracks) |>
    select("source", "ploidy", "Gemcitabine", "dose_label", "site_id", "migration_track_id", "frame", "dx", "dy", "speed")
}

mark_proactive_mechanistic_censoring <- function(tracks, reference_tracks, distance_threshold) {
  if (nrow(tracks) == 0L) {
    return(tracks |> mutate(proactive_censored = logical(), nearest_reference_object_distance = numeric()))
  }
  reference_light <- reference_tracks |>
    select("ploidy", "Gemcitabine", "site_id", "frame", "row_id", "x", "y")

  tracks |>
    select(-any_of(c("nearest_reference_object_distance", "proactive_censored"))) |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$site_id, .data$frame) |>
    group_modify(~ {
      ref <- reference_light |>
        filter(
          .data$ploidy == .y$ploidy[[1]],
          .data$Gemcitabine == .y$Gemcitabine[[1]],
          .data$site_id == .y$site_id[[1]],
          .data$frame == .y$frame[[1]]
        )
      if (nrow(ref) < 2L || nrow(.x) == 0L) {
        return(.x |> mutate(nearest_reference_object_distance = Inf, proactive_censored = FALSE))
      }
      dx <- outer(.x$x, ref$x, "-")
      dy <- outer(.x$y, ref$y, "-")
      distance <- sqrt(dx^2 + dy^2)
      same_object <- outer(.x$row_id, ref$row_id, "==")
      distance[same_object] <- Inf
      nearest_distance <- apply(distance, 1L, min)
      .x |>
        mutate(
          nearest_reference_object_distance = nearest_distance,
          proactive_censored = nearest_distance < distance_threshold
        )
    }) |>
    ungroup()
}

summarize_track_counts <- function(tracks) {
  tracks |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label) |>
    summarize(
      n_tracks = n_distinct(.data$migration_track_id),
      n_points = n(),
      .groups = "drop"
    )
}

mechanistic_worker <- function(artifact, disc_radius, max_lag, min_segment_frames, step_size_quantile, max_tracks, seed) {
  set.seed(seed)
  base <- artifact$simulated_uncensored
  if (is.null(base) || nrow(base) == 0L) {
    return(NULL)
  }
  n_available_tracks <- n_distinct(base$migration_track_id)
  keep_ids <- base |>
    distinct(.data$migration_track_id) |>
    slice_sample(n = min(max_tracks, n_available_tracks), replace = FALSE) |>
    pull("migration_track_id")
  base <- base |>
    filter(.data$migration_track_id %in% keep_ids) |>
    mutate(
      source = "OU fitted raw",
      ploidy = as.character(.data$ploidy),
      Gemcitabine = as.numeric(.data$Gemcitabine),
      dose_label = as.character(.data$dose_label),
      site_id = as.character(.data$site_id),
      migration_track_id = as.character(.data$migration_track_id),
      frame = as.integer(.data$frame),
      track_step = as.integer(.data$track_step),
      x = as.numeric(.data$x),
      y = as.numeric(.data$y),
      width = as.numeric(.data$width),
      height = as.numeric(.data$height)
    )

  marked <- mark_mechanistic_censoring(base, disc_radius = disc_radius)
  boundary <- marked |>
    filter(!.data$boundary_censored) |>
    mutate(source = "OU fitted boundary") |>
    split_mechanistic_segments(min_segment_frames = min_segment_frames)
  overlap <- marked |>
    filter(!.data$overlap_censored) |>
    mutate(source = "OU fitted overlap") |>
    split_mechanistic_segments(min_segment_frames = min_segment_frames)
  mechanical <- marked |>
    filter(!.data$is_mechanically_censored) |>
    mutate(source = "OU fitted mechanical") |>
    split_mechanistic_segments(min_segment_frames = min_segment_frames)

  mechanical_steps <- prepare_mechanistic_steps(mechanical)
  step_size_threshold <- mechanical_steps |>
    summarize(threshold = quantile(.data$speed, probs = step_size_quantile, na.rm = TRUE, names = FALSE)) |>
    pull("threshold")
  if (!is.finite(step_size_threshold)) {
    step_size_threshold <- 0
  }
  proactive_threshold <- 2 * disc_radius + step_size_threshold
  proactive_marked <- mark_proactive_mechanistic_censoring(
    mechanical,
    reference_tracks = marked,
    distance_threshold = proactive_threshold
  )
  proactive <- proactive_marked |>
    filter(!.data$proactive_censored) |>
    mutate(source = "OU fitted mechanical + proactive") |>
    split_mechanistic_segments(min_segment_frames = min_segment_frames)

  metric_tracks <- bind_rows(
    base |> select(any_of(c(track_cols, "width", "height"))) |> select(all_of(track_cols)),
    boundary |> select(all_of(track_cols)),
    overlap |> select(all_of(track_cols)),
    mechanical |> select(all_of(track_cols)),
    proactive |> select(all_of(track_cols))
  )

  metric_groups <- metric_tracks |>
    mutate(.split_key = paste(.data$source, .data$ploidy, .data$Gemcitabine, sep = "\r")) |>
    group_split(.data$.split_key, .keep = FALSE)

  list(
    msd = bind_rows(lapply(metric_groups, summarize_msd, max_lag = max_lag)),
    velocity_autocorrelation = bind_rows(lapply(
      metric_groups,
      function(x) summarize_velocity_autocorrelation(compute_steps(x), max_lag = max_lag)
    )),
    track_lengths = bind_rows(lapply(metric_groups, summarize_track_lengths)),
    speeds = compute_steps(metric_tracks),
    track_counts = summarize_track_counts(metric_tracks),
    censor_summary = marked |>
      group_by(.data$ploidy, .data$Gemcitabine, .data$dose_label, .data$censor_reason) |>
      summarize(n_points = n(), n_tracks = n_distinct(.data$migration_track_id), .groups = "drop"),
    proactive_summary = tibble(
      ploidy = base$ploidy[[1]],
      Gemcitabine = base$Gemcitabine[[1]],
      dose_label = base$dose_label[[1]],
      step_size_quantile = step_size_quantile,
      step_size_threshold = step_size_threshold,
      proactive_distance_threshold = proactive_threshold,
      mechanical_points = nrow(mechanical),
      proactive_removed_points = sum(proactive_marked$proactive_censored, na.rm = TRUE)
    )
  )
}

message("Loading simulation summary: ", simulation_rds)
simulation <- readRDS(simulation_rds)

condition_files <- sort(list.files(condition_dir, pattern = "^ou_censoring_simulation_.*_seed[0-9]+[.]rds$", full.names = TRUE))
if (length(condition_files) == 0L) {
  stop("No condition simulation artifacts found in ", condition_dir, call. = FALSE)
}
message("Loading ", length(condition_files), " condition artifacts")
condition_artifacts <- lapply(condition_files, readRDS)
conditions <- bind_rows(lapply(condition_artifacts, `[[`, "condition")) |>
  distinct(.data$ploidy, .data$Gemcitabine) |>
  mutate(ploidy = as.character(.data$ploidy), Gemcitabine = as.numeric(.data$Gemcitabine))

message("Preparing compact comparison tracks")
observed_tracks <- prepare_observed_tracks(tracks_rds, conditions)
uncensored_tracks <- bind_rows(lapply(condition_artifacts, function(x) select_tracks(x$simulated_uncensored)))
context_tracks <- bind_rows(lapply(condition_artifacts, function(x) select_tracks(x$simulated_context_censored)))
full_tracks <- bind_rows(lapply(condition_artifacts, function(x) select_tracks(x$simulated_full_censored)))
comparison_tracks <- bind_rows(observed_tracks, uncensored_tracks, context_tracks, full_tracks)

tracks_raw_for_area <- readRDS(tracks_rds)
disc_area <- if ("Size_in_pixels_0" %in% names(tracks_raw_for_area)) {
  median(as.numeric(tracks_raw_for_area$Size_in_pixels_0), na.rm = TRUE)
} else if ("track_mask_area_px" %in% names(tracks_raw_for_area)) {
  median(as.numeric(tracks_raw_for_area$track_mask_area_px), na.rm = TRUE)
} else {
  450
}
if (!is.finite(disc_area) || disc_area <= 0) {
  disc_area <- 450
}
disc_radius <- sqrt(disc_area / pi)

if (skip_cosine && file.exists(out_rds)) {
  message("Reusing cosine autocorrelation addenda from existing artifact: ", out_rds)
  existing_addenda <- readRDS(out_rds)
  cosine_autocorrelation_summary <- existing_addenda$cosine_autocorrelation_summary
  cosine_autocorrelation_distribution <- existing_addenda$cosine_autocorrelation_distribution
} else {
  message("Computing cosine autocorrelation pairs with ", parallel_cores, " core(s)")
  track_groups <- comparison_tracks |>
    mutate(.split_key = paste(.data$source, .data$ploidy, .data$Gemcitabine, sep = "\r")) |>
    group_split(.data$.split_key, .keep = FALSE)
  cosine_pairs <- bind_rows(parallel::mclapply(
    track_groups,
    compute_cosine_pairs,
    max_lag = max_lag,
    mc.cores = min(parallel_cores, length(track_groups)),
    mc.preschedule = FALSE
  ))

  cosine_autocorrelation_summary <- summarize_cosine_pairs(cosine_pairs)
  cosine_autocorrelation_distribution <- bind_rows(
    sample_cosine_pairs(cosine_pairs, sample_per_group, seed = seed + 1L),
    random_walk_comparator(cosine_pairs, sample_per_group, seed = seed + 2L)
  )
}

source_overlap <- simulation$metrics$track_counts |>
  select("source", "ploidy", "Gemcitabine", "dose_label", "n_tracks", "n_points") |>
  pivot_wider(
    names_from = "source",
    values_from = c("n_tracks", "n_points"),
    values_fill = 0
  ) |>
  mutate(
    full_tracks_vs_uncensored = .data$`n_tracks_OU + full censoring` / pmax(.data$`n_tracks_OU uncensored`, 1),
    full_points_vs_uncensored = .data$`n_points_OU + full censoring` / pmax(.data$`n_points_OU uncensored`, 1),
    observed_points_vs_full = .data$n_points_Observed / pmax(.data$`n_points_OU + full censoring`, 1)
  )

msd_similarity <- simulation$metrics$msd |>
  filter(.data$source %in% c("Observed", "OU uncensored", "OU + full censoring")) |>
  select("source", "ploidy", "Gemcitabine", "dose_label", "lag_frames", "mean_squared_displacement_px2") |>
  pivot_wider(names_from = "source", values_from = "mean_squared_displacement_px2") |>
  mutate(
    full_minus_uncensored = .data$`OU + full censoring` - .data$`OU uncensored`,
    observed_minus_full = .data$Observed - .data$`OU + full censoring`,
    full_over_uncensored = .data$`OU + full censoring` / .data$`OU uncensored`,
    observed_over_full = .data$Observed / .data$`OU + full censoring`
  )

message(
  "Computing mechanistic OU censoring summaries with disc area ",
  round(disc_area, 2), " px^2 and ", mechanistic_tracks_per_condition,
  " track(s) per condition"
)
mechanistic_results <- parallel::mclapply(
  seq_along(condition_artifacts),
  function(i) {
    mechanistic_worker(
      condition_artifacts[[i]],
      disc_radius = disc_radius,
      max_lag = max_lag,
      min_segment_frames = mechanistic_min_segment_frames,
      step_size_quantile = mechanistic_step_size_quantile,
      max_tracks = mechanistic_tracks_per_condition,
      seed = seed + 10000L + i
    )
  },
  mc.cores = min(parallel_cores, length(condition_artifacts)),
  mc.preschedule = FALSE
)
mechanistic_results <- Filter(Negate(is.null), mechanistic_results)

mechanistic <- list(
  disc_area = disc_area,
  disc_radius = disc_radius,
  tracks_per_condition = mechanistic_tracks_per_condition,
  step_size_quantile = mechanistic_step_size_quantile,
  min_segment_frames = mechanistic_min_segment_frames,
  msd = bind_rows(lapply(mechanistic_results, `[[`, "msd")),
  velocity_autocorrelation = bind_rows(lapply(mechanistic_results, `[[`, "velocity_autocorrelation")),
  track_lengths = bind_rows(lapply(mechanistic_results, `[[`, "track_lengths")),
  speeds = bind_rows(lapply(mechanistic_results, `[[`, "speeds")),
  track_counts = bind_rows(lapply(mechanistic_results, `[[`, "track_counts")),
  censor_summary = bind_rows(lapply(mechanistic_results, `[[`, "censor_summary")),
  proactive_summary = bind_rows(lapply(mechanistic_results, `[[`, "proactive_summary"))
)

artifact <- list(
  metadata = list(
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    simulation_rds = simulation_rds,
    condition_dir = condition_dir,
    tracks_rds = tracks_rds,
    selected_doses = selected_doses,
    max_lag = max_lag,
    sample_per_group = sample_per_group,
    mechanistic_tracks_per_condition = mechanistic_tracks_per_condition,
    mechanistic_step_size_quantile = mechanistic_step_size_quantile,
    mechanistic_min_segment_frames = mechanistic_min_segment_frames,
    parallel_cores = parallel_cores,
    seed = seed
  ),
  cosine_autocorrelation_summary = cosine_autocorrelation_summary,
  cosine_autocorrelation_distribution = cosine_autocorrelation_distribution,
  source_overlap = source_overlap,
  msd_similarity = msd_similarity,
  mechanistic = mechanistic
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(artifact, out_rds)
message("Saved report addenda: ", out_rds)
