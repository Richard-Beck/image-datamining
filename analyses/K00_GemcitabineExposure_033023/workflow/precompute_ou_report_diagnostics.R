#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
  library(readr)
  library(tibble)
  library(tidyr)
})

usage <- paste0(
  "Usage: precompute_ou_report_diagnostics.R [options]\n\n",
  "Options:\n",
  "  --repo_root=/path/to/repo\n",
  "  --tracks_rds=/path/to/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds\n",
  "  --fits_csv=/path/to/ou_tracking_fits_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv\n",
  "  --out_rds=/path/to/ou_report_diagnostics_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds\n",
  "  --cores=16\n",
  "  --max_lag=10\n",
  "  --track_plot_doses=0,6.25,25\n",
  "  --max_tracks_per_condition=350\n",
  "  --conditional_sample_per_group=2000\n",
  "  --seed=20260508\n",
  "  --quiet=TRUE\n"
)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
analysis_dir_guess <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(analysis_dir_guess, "R/k00_batch_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
repo_root <- normalizePath(args$repo_root %||% normalizePath(file.path(analysis_dir_guess, "../.."), mustWork = TRUE), mustWork = TRUE)
analysis_dir <- file.path(repo_root, "analyses/K00_GemcitabineExposure_033023")
data_dir <- file.path(analysis_dir, "data")

tracks_rds <- normalizePath(
  args$tracks_rds %||% file.path(data_dir, "tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"),
  mustWork = TRUE
)
fits_csv <- normalizePath(
  args$fits_csv %||% file.path(data_dir, "ou_tracking_fits_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv"),
  mustWork = TRUE
)
out_rds <- normalizePath(
  args$out_rds %||% file.path(data_dir, "ou_report_diagnostics_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"),
  mustWork = FALSE
)
cores <- as.integer(args$cores %||% "16")
max_lag <- as.integer(args$max_lag %||% "10")
track_plot_doses <- as.numeric(strsplit(args$track_plot_doses %||% "0,6.25,25", ",", fixed = TRUE)[[1]])
max_tracks_per_condition <- as.integer(args$max_tracks_per_condition %||% "350")
conditional_sample_per_group <- as.integer(args$conditional_sample_per_group %||% "2000")
seed <- as.integer(args$seed %||% "20260508")
quiet <- as_flag(args$quiet, default = FALSE)

if (is.na(cores) || cores < 1L) {
  stop("--cores must be a positive integer", call. = FALSE)
}
if (is.na(max_lag) || max_lag < 1L) {
  stop("--max_lag must be a positive integer", call. = FALSE)
}

order_dose_label <- function(label, dose) {
  dose_levels <- sort(unique(dose[is.finite(dose)]))
  factor(label, levels = paste0(dose_levels, " nM"))
}

format_ploidy <- function(x) {
  factor(x, levels = sort(unique(x)))
}

time_window_label <- function(frame) {
  case_when(
    frame >= 0 & frame < 10 ~ "frames 0-10",
    frame >= 10 & frame < 20 ~ "frames 10-20",
    frame >= 20 & frame <= 40 ~ "frames 20-40",
    TRUE ~ NA_character_
  )
}

ou_transition <- function(dt, tau, velocity_scale) {
  a <- dt / tau
  phi <- exp(-a)
  one_minus_phi <- -expm1(-a)
  one_minus_phi2 <- -expm1(-2 * a)
  f12 <- tau * one_minus_phi
  q22 <- velocity_scale^2 * one_minus_phi2
  q12 <- velocity_scale^2 * tau * one_minus_phi^2
  q11_unit <- if (a < 1e-4) {
    a^3 / 3 - a^4 / 4 + 7 * a^5 / 60
  } else {
    a - 2 * one_minus_phi + 0.5 * one_minus_phi2
  }
  q11 <- 2 * velocity_scale^2 * tau^2 * q11_unit

  q11 <- max(q11, 0)
  q22 <- max(q22, 0)
  if (q22 > 0 && q11 * q22 < q12^2) {
    q11 <- q12^2 / q22 + .Machine$double.eps
  }

  list(
    F = matrix(c(1, 0, f12, phi), nrow = 2),
    Q = matrix(c(q11, q12, q12, q22), nrow = 2)
  )
}

simulate_ou_track <- function(tau, velocity_scale, obs_noise, n_steps, first_frame, dt = 1) {
  transition <- ou_transition(dt = dt, tau = tau, velocity_scale = velocity_scale)
  noise_chol <- chol(transition$Q)

  state_x <- c(0, rnorm(1, sd = velocity_scale))
  state_y <- c(0, rnorm(1, sd = velocity_scale))
  rows <- vector("list", n_steps + 1L)
  rows[[1L]] <- tibble(
    track_step = 0L,
    frame = first_frame,
    latent_x = state_x[[1]],
    latent_y = state_y[[1]],
    x = state_x[[1]] + rnorm(1, sd = obs_noise),
    y = state_y[[1]] + rnorm(1, sd = obs_noise)
  )

  for (step in seq_len(n_steps)) {
    state_x <- drop(transition$F %*% state_x + drop(t(noise_chol) %*% rnorm(2)))
    state_y <- drop(transition$F %*% state_y + drop(t(noise_chol) %*% rnorm(2)))
    rows[[step + 1L]] <- tibble(
      track_step = step,
      frame = first_frame + step,
      latent_x = state_x[[1]],
      latent_y = state_y[[1]],
      x = state_x[[1]] + rnorm(1, sd = obs_noise),
      y = state_y[[1]] + rnorm(1, sd = obs_noise)
    )
  }

  bind_rows(rows)
}

compute_steps <- function(tracks, source, x_col = "x", y_col = "y") {
  tracks |>
    arrange(.data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$site_id, .data$migration_track_id) |>
    mutate(
      next_frame = lead(.data$frame),
      dx = lead(.data[[x_col]]) - .data[[x_col]],
      dy = lead(.data[[y_col]]) - .data[[y_col]],
      step_frames = .data$next_frame - .data$frame,
      step_distance_px = sqrt(.data$dx^2 + .data$dy^2),
      frame_window = time_window_label(.data$frame)
    ) |>
    ungroup() |>
    filter(.data$step_frames == 1L, is.finite(.data$step_distance_px)) |>
    mutate(source = source)
}

summarize_origin_msd <- function(tracks, source, max_lag) {
  tracks |>
    arrange(.data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$site_id, .data$migration_track_id) |>
    mutate(
      track_step = row_number() - 1L,
      origin_x = first(.data$x),
      origin_y = first(.data$y),
      squared_displacement_px2 = (.data$x - .data$origin_x)^2 + (.data$y - .data$origin_y)^2
    ) |>
    ungroup() |>
    filter(.data$track_step <= max_lag) |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$track_step) |>
    summarize(
      mean_squared_displacement_px2 = mean(.data$squared_displacement_px2, na.rm = TRUE),
      n_tracks = n_distinct(interaction(.data$site_id, .data$migration_track_id, drop = TRUE)),
      .groups = "drop"
    ) |>
    mutate(source = source)
}

summarize_windowed_msd <- function(tracks, source, max_lag) {
  tracks |>
    arrange(.data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$site_id, .data$migration_track_id) |>
    group_modify(~ {
      df <- .x |> arrange(.data$frame)
      bind_rows(lapply(seq_len(max_lag), function(lag_frame) {
        j <- match(df$frame + lag_frame, df$frame)
        ok <- !is.na(j)
        if (!any(ok)) {
          return(NULL)
        }
        tibble(
          frame = df$frame[ok],
          frame_window = time_window_label(df$frame[ok]),
          lag_frames = lag_frame,
          squared_displacement_px2 = (df$x[j[ok]] - df$x[ok])^2 + (df$y[j[ok]] - df$y[ok])^2
        )
      }))
    }) |>
    ungroup() |>
    filter(!is.na(.data$frame_window)) |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$frame_window, .data$lag_frames) |>
    summarize(
      mean_squared_displacement_px2 = mean(.data$squared_displacement_px2, na.rm = TRUE),
      n_displacements = n(),
      .groups = "drop"
    ) |>
    mutate(source = source)
}

compute_cos_autocorrelation <- function(steps, max_lag) {
  steps |>
    arrange(.data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$source, .data$site_id, .data$migration_track_id) |>
    group_modify(~ {
      df <- .x |> arrange(.data$frame)
      bind_rows(lapply(seq_len(max_lag), function(lag_frame) {
        j <- match(df$frame + lag_frame, df$frame)
        ok <- !is.na(j)
        if (!any(ok)) {
          return(NULL)
        }
        tibble(
          frame = df$frame[ok],
          frame_window = time_window_label(df$frame[ok]),
          lag_frames = lag_frame,
          velocity_cos_autocorrelation = (df$dx[ok] * df$dx[j[ok]] + df$dy[ok] * df$dy[j[ok]]) /
            pmax(df$step_distance_px[ok] * df$step_distance_px[j[ok]], 1e-9)
        )
      }))
    }) |>
    ungroup()
}

summarize_cos_autocorrelation <- function(autocorrelation, include_window = FALSE) {
  group_cols <- c("source", "ploidy", "Gemcitabine", "gemcitabine_label")
  if (include_window) {
    group_cols <- c(group_cols, "frame_window")
  }
  group_cols <- c(group_cols, "lag_frames")

  autocorrelation |>
    filter(if (include_window) !is.na(.data$frame_window) else TRUE) |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      mean_velocity_cos_autocorrelation = mean(.data$velocity_cos_autocorrelation, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    )
}

conditional_endpoint_sample <- function(steps, selected_doses, sample_per_group, seed) {
  set.seed(seed)
  endpoints <- steps |>
    filter(.data$Gemcitabine %in% selected_doses) |>
    arrange(.data$source, .data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$site_id, .data$migration_track_id) |>
    mutate(
      next_step_frame = lead(.data$frame),
      dx_next = lead(.data$dx),
      dy_next = lead(.data$dy),
      current_step_length_bin = cut(
        .data$step_distance_px,
        breaks = c(-Inf, 5, 10, 20, Inf),
        labels = c("0-5 px", "5-10 px", "10-20 px", ">=20 px")
      ),
      current_angle = atan2(.data$dy, .data$dx),
      next_dx_aligned = cos(-.data$current_angle) * .data$dx_next - sin(-.data$current_angle) * .data$dy_next,
      next_dy_aligned = sin(-.data$current_angle) * .data$dx_next + cos(-.data$current_angle) * .data$dy_next
    ) |>
    ungroup() |>
    filter(
      .data$next_step_frame == .data$frame + 1L,
      !is.na(.data$current_step_length_bin),
      is.finite(.data$next_dx_aligned),
      is.finite(.data$next_dy_aligned)
    )

  endpoints |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$current_step_length_bin) |>
    group_modify(~ {
      n <- min(nrow(.x), sample_per_group)
      .x[sample.int(nrow(.x), n), , drop = FALSE]
    }) |>
    ungroup() |>
    select(
      source, ploidy, Gemcitabine, gemcitabine_label, current_step_length_bin,
      next_dx_aligned, next_dy_aligned
    )
}

simulate_tracks_for_condition <- function(track_summary, fit_row, seed_i) {
  set.seed(seed_i)
  bind_rows(lapply(seq_len(nrow(track_summary)), function(i) {
    simulate_ou_track(
      tau = fit_row$tau[[1]],
      velocity_scale = fit_row$velocity_scale[[1]],
      obs_noise = fit_row$obs_noise[[1]],
      n_steps = track_summary$n_steps[[i]],
      first_frame = track_summary$first_frame[[i]],
      dt = fit_row$frame_interval[[1]]
    ) |>
      mutate(
        ploidy = track_summary$ploidy[[i]],
        Gemcitabine = track_summary$Gemcitabine[[i]],
        gemcitabine_label = track_summary$gemcitabine_label[[i]],
        site_id = track_summary$site_id[[i]],
        migration_track_id = track_summary$migration_track_id[[i]]
      )
  }))
}

message_if <- function(...) {
  if (!quiet) {
    message(...)
  }
}

message_if("Reading fits: ", fits_csv)
ou_fits_raw <- read_csv(fits_csv, show_col_types = FALSE)
for (missing_col in setdiff(
  c(
    "n_starts", "seed", "start_id", "best_start_id", "start_tau",
    "start_velocity_scale", "start_obs_noise", "delta_loglik_next_best",
    "n_converged_starts"
  ),
  names(ou_fits_raw)
)) {
  ou_fits_raw[[missing_col]] <- NA
}
ou_fits <- ou_fits_raw |>
  mutate(
    dataset = "yellow_reconstructed_area2x_nonnegative_tracks",
    ploidy = format_ploidy(.data$ploidy),
    gemcitabine_label = paste0(.data$Gemcitabine, " nM"),
    gemcitabine_label = order_dose_label(.data$gemcitabine_label, .data$Gemcitabine),
    clean_fit_status = if_else(grepl("^converged", .data$fit_status), "converged", .data$fit_status),
    log_likelihood_per_segment = .data$log_likelihood / .data$n_segments,
    log_likelihood_per_track_time = .data$log_likelihood / .data$total_track_time
  ) |>
  arrange(.data$job_id)

best_fit_summary <- ou_fits |>
  group_by(.data$dataset, .data$job_id) |>
  summarize(
    best_log_likelihood = max(.data$log_likelihood, na.rm = TRUE),
    second_best_log_likelihood = suppressWarnings(sort(.data$log_likelihood, decreasing = TRUE, na.last = NA)[2]),
    n_converged_starts = sum(.data$clean_fit_status == "converged", na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    second_best_log_likelihood = if_else(is.finite(.data$second_best_log_likelihood), .data$second_best_log_likelihood, NA_real_),
    delta_loglik_next_best = .data$best_log_likelihood - .data$second_best_log_likelihood
  )

best_fits <- ou_fits |>
  group_by(.data$dataset, .data$job_id) |>
  slice_max(order_by = .data$log_likelihood, n = 1, with_ties = FALSE) |>
  ungroup() |>
  left_join(best_fit_summary, by = c("dataset", "job_id")) |>
  mutate(
    best_start_id = .data$start_id,
    delta_loglik_next_best = coalesce(.data$delta_loglik_next_best.y, .data$delta_loglik_next_best.x),
    n_converged_starts = coalesce(.data$n_converged_starts.y, .data$n_converged_starts.x)
  ) |>
  select(-ends_with(".x"), -ends_with(".y"))

message_if("Reading tracks: ", tracks_rds)
tracks <- readRDS(tracks_rds) |>
  transmute(
    ploidy = format_ploidy(.data$ploidy),
    Gemcitabine = .data$Gemcitabine,
    gemcitabine_label = order_dose_label(paste0(.data$Gemcitabine, " nM"), .data$Gemcitabine),
    well = .data$well,
    position = .data$position,
    site_id = .data$site_id,
    migration_track_id = .data$migration_track_id,
    frame = .data$frame,
    x = .data$Center_of_the_object_1,
    y = .data$Center_of_the_object_0
  )

track_summary <- tracks |>
  group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$well, .data$position, .data$site_id, .data$migration_track_id) |>
  summarize(
    n_frames = n_distinct(.data$frame),
    first_frame = min(.data$frame, na.rm = TRUE),
    last_frame = max(.data$frame, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(n_steps = pmax(.data$n_frames - 1L, 0L))

condition_keys <- track_summary |>
  distinct(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label) |>
  arrange(.data$ploidy, .data$Gemcitabine)

observed_jobs <- lapply(seq_len(nrow(condition_keys)), function(i) {
  key <- condition_keys[i, , drop = FALSE]
  tracks |>
    filter(.data$ploidy == key$ploidy[[1]], .data$Gemcitabine == key$Gemcitabine[[1]])
})

observed_worker <- function(condition_tracks) {
  condition_steps <- compute_steps(condition_tracks, source = "Observed")
  list(
    steps = condition_steps,
    msd = summarize_origin_msd(condition_tracks, source = "Observed", max_lag = max_lag),
    windowed_msd = summarize_windowed_msd(condition_tracks, source = "Observed", max_lag = max_lag),
    autocorrelation = compute_cos_autocorrelation(condition_steps, max_lag = max_lag)
  )
}

message_if("Preparing observed summaries with ", cores, " core(s)")
observed_results <- parallel::mclapply(
  observed_jobs,
  observed_worker,
  mc.cores = min(cores, length(observed_jobs))
)
observed_steps <- bind_rows(lapply(observed_results, `[[`, "steps"))
observed_msd <- bind_rows(lapply(observed_results, `[[`, "msd"))
observed_windowed_msd <- bind_rows(lapply(observed_results, `[[`, "windowed_msd"))
observed_autocorrelation <- bind_rows(lapply(observed_results, `[[`, "autocorrelation"))

message_if("Simulating matched OU tracks with ", cores, " core(s)")
condition_jobs <- lapply(seq_len(nrow(condition_keys)), function(i) {
  key <- condition_keys[i, , drop = FALSE]
  list(
    key = key,
    tracks = track_summary |>
      filter(.data$ploidy == key$ploidy[[1]], .data$Gemcitabine == key$Gemcitabine[[1]], .data$n_steps > 0L),
    fit = best_fits |>
      filter(.data$ploidy == key$ploidy[[1]], .data$Gemcitabine == key$Gemcitabine[[1]]) |>
      slice_head(n = 1),
    seed_i = seed + i
  )
})

worker <- function(job) {
  simulate_tracks_for_condition(job$tracks, job$fit, job$seed_i)
}

simulated_tracks <- bind_rows(parallel::mclapply(
  condition_jobs,
  worker,
  mc.cores = min(cores, length(condition_jobs))
)) |>
  arrange(.data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
  group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$site_id, .data$migration_track_id) |>
  mutate(
    origin_x = first(.data$x),
    origin_y = first(.data$y),
    centered_x = .data$x - .data$origin_x,
    centered_y = .data$y - .data$origin_y
  ) |>
  ungroup()

simulated_steps <- compute_steps(simulated_tracks, source = "OU simulated")
simulated_msd <- summarize_origin_msd(simulated_tracks, source = "OU simulated", max_lag = max_lag)
simulated_windowed_msd <- summarize_windowed_msd(simulated_tracks, source = "OU simulated", max_lag = max_lag)
simulated_autocorrelation <- compute_cos_autocorrelation(simulated_steps, max_lag = max_lag)

message_if("Preparing plot tables")
sampled_track_ids <- track_summary |>
  filter(.data$Gemcitabine %in% track_plot_doses, .data$n_steps > 0L) |>
  arrange(.data$ploidy, .data$Gemcitabine, desc(.data$n_frames), .data$migration_track_id) |>
  group_by(.data$ploidy, .data$Gemcitabine) |>
  slice_head(n = max_tracks_per_condition) |>
  ungroup()

observed_sampled_tracks <- tracks |>
  semi_join(
    sampled_track_ids |> select(ploidy, Gemcitabine, site_id, migration_track_id),
    by = c("ploidy", "Gemcitabine", "site_id", "migration_track_id")
  ) |>
  arrange(.data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
  group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$site_id, .data$migration_track_id) |>
  mutate(
    track_step = row_number() - 1L,
    origin_x = first(.data$x),
    origin_y = first(.data$y),
    centered_x = .data$x - .data$origin_x,
    centered_y = .data$y - .data$origin_y,
    source = "Observed"
  ) |>
  ungroup() |>
  select(source, ploidy, Gemcitabine, gemcitabine_label, site_id, migration_track_id, track_step, centered_x, centered_y)

simulated_sampled_tracks <- simulated_tracks |>
  semi_join(
    sampled_track_ids |> select(ploidy, Gemcitabine, site_id, migration_track_id),
    by = c("ploidy", "Gemcitabine", "site_id", "migration_track_id")
  ) |>
  mutate(source = "OU simulated") |>
  select(source, ploidy, Gemcitabine, gemcitabine_label, site_id, migration_track_id, track_step, centered_x, centered_y)

centered_track_segments <- tracks |>
  filter(.data$Gemcitabine %in% track_plot_doses) |>
  arrange(.data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
  group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label, .data$site_id, .data$migration_track_id) |>
  mutate(
    track_observation_index = row_number(),
    track_observations = n(),
    centered_x = .data$x - first(.data$x),
    centered_y = .data$y - first(.data$y)
  ) |>
  ungroup() |>
  filter(.data$track_observations >= 10L, .data$track_observation_index <= 10L) |>
  select(ploidy, Gemcitabine, gemcitabine_label, site_id, migration_track_id, track_observation_index, centered_x, centered_y)

msd_by_step <- bind_rows(observed_msd, simulated_msd)
windowed_msd_by_lag <- bind_rows(observed_windowed_msd, simulated_windowed_msd)

ou_msd_curve <- best_fits |>
  select(ploidy, Gemcitabine, gemcitabine_label, tau, velocity_scale, obs_noise, frame_interval) |>
  tidyr::expand_grid(track_step = 1:max_lag) |>
  mutate(
    source = "OU expectation",
    lag_time = .data$track_step * .data$frame_interval,
    mean_squared_displacement_px2 = 4 * .data$velocity_scale^2 * .data$tau^2 *
      (.data$lag_time / .data$tau - 1 + exp(-.data$lag_time / .data$tau)) +
      4 * .data$obs_noise^2
  )

step_distribution <- bind_rows(observed_steps, simulated_steps) |>
  filter(.data$Gemcitabine %in% track_plot_doses) |>
  transmute(
    source,
    ploidy,
    Gemcitabine,
    gemcitabine_label,
    frame_window,
    step_distance_px,
    capped_step_distance_px = pmin(.data$step_distance_px, 50)
  )

cos_autocorrelation_by_lag <- summarize_cos_autocorrelation(
  bind_rows(observed_autocorrelation, simulated_autocorrelation),
  include_window = FALSE
)
windowed_cos_autocorrelation_by_lag <- summarize_cos_autocorrelation(
  bind_rows(observed_autocorrelation, simulated_autocorrelation),
  include_window = TRUE
)

track_length_distribution <- track_summary |>
  mutate(
    track_length_bin = factor(
      case_when(
        .data$n_frames >= 10L ~ ">=10",
        .data$n_frames >= 7L ~ "7-9",
        .data$n_frames >= 5L ~ "5-6",
        TRUE ~ "3-4"
      ),
      levels = c("3-4", "5-6", "7-9", ">=10")
    )
  ) |>
  count(ploidy, Gemcitabine, gemcitabine_label, track_length_bin, name = "n_tracks")

threshold_sensitivity <- bind_rows(lapply(c(3L, 4L, 5L), function(min_frames) {
  eligible_tracks <- track_summary |> filter(.data$n_frames >= min_frames)
  eligible_steps <- observed_steps |>
    semi_join(
      eligible_tracks |> select(ploidy, Gemcitabine, site_id, migration_track_id),
      by = c("ploidy", "Gemcitabine", "site_id", "migration_track_id")
    )
  eligible_msd <- tracks |>
    semi_join(
      eligible_tracks |> select(ploidy, Gemcitabine, site_id, migration_track_id),
      by = c("ploidy", "Gemcitabine", "site_id", "migration_track_id")
    ) |>
    summarize_origin_msd(source = "Observed", max_lag = max_lag) |>
    filter(.data$track_step %in% c(3L, 5L, 10L)) |>
    select(ploidy, Gemcitabine, gemcitabine_label, track_step, mean_squared_displacement_px2) |>
    pivot_wider(
      names_from = track_step,
      values_from = mean_squared_displacement_px2,
      names_prefix = "msd_step_"
    )

  eligible_tracks |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label) |>
    summarize(
      n_tracks = n(),
      median_track_frames = median(.data$n_frames, na.rm = TRUE),
      p90_track_frames = quantile(.data$n_frames, 0.90, na.rm = TRUE),
      .groups = "drop"
    ) |>
    left_join(
      eligible_steps |>
        group_by(.data$ploidy, .data$Gemcitabine, .data$gemcitabine_label) |>
        summarize(
          n_steps = n(),
          median_step_distance_px = median(.data$step_distance_px, na.rm = TRUE),
          mean_step_distance_px = mean(.data$step_distance_px, na.rm = TRUE),
          .groups = "drop"
        ),
      by = c("ploidy", "Gemcitabine", "gemcitabine_label")
    ) |>
    left_join(eligible_msd, by = c("ploidy", "Gemcitabine", "gemcitabine_label")) |>
    mutate(min_track_frames = min_frames, .before = 1L)
}))

conditional_endpoints <- conditional_endpoint_sample(
  bind_rows(observed_steps, simulated_steps),
  selected_doses = track_plot_doses,
  sample_per_group = conditional_sample_per_group,
  seed = seed + 1000L
)

diagnostics <- list(
  metadata = list(
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    tracks_rds = tracks_rds,
    fits_csv = fits_csv,
    max_lag = max_lag,
    track_plot_doses = track_plot_doses,
    max_tracks_per_condition = max_tracks_per_condition,
    conditional_sample_per_group = conditional_sample_per_group,
    seed = seed,
    cores = cores
  ),
  track_summary = track_summary,
  centered_track_segments = centered_track_segments,
  observed_sampled_tracks = observed_sampled_tracks,
  simulated_sampled_tracks = simulated_sampled_tracks,
  msd_by_step = msd_by_step,
  ou_msd_curve = ou_msd_curve,
  windowed_msd_by_lag = windowed_msd_by_lag,
  step_distribution = step_distribution,
  cos_autocorrelation_by_lag = cos_autocorrelation_by_lag,
  windowed_cos_autocorrelation_by_lag = windowed_cos_autocorrelation_by_lag,
  conditional_endpoints = conditional_endpoints,
  track_length_distribution = track_length_distribution,
  threshold_sensitivity = threshold_sensitivity
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(diagnostics, out_rds)
message_if("Wrote ", out_rds)
