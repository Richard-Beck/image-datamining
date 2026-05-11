#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Arguments must be passed as --name=value: ", arg, call. = FALSE)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[parts[[1L]]]] <- paste(parts[-1L], collapse = "=")
  }
  out
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

script_file <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]),
  mustWork = TRUE
)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

analysis_dir <- normalizePath(
  args$analysis_dir %||% file.path(dirname(script_file), "../../.."),
  mustWork = TRUE
)
source(file.path(analysis_dir, "R/ou_velocity_model.R"))

tracks_rds <- normalizePath(
  args$tracks_rds %||% file.path(
    analysis_dir,
    "data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"
  ),
  mustWork = TRUE
)
out_rds <- normalizePath(
  args$out_rds %||% file.path(
    analysis_dir,
    "dev/alternative_motility_models/artifacts/alternative_motility_model_fits.rds"
  ),
  mustWork = FALSE
)
selected_doses <- as.numeric(strsplit(args$selected_doses %||% "0,3.125,6.25,12.5,25", ",", fixed = TRUE)[[1L]])
max_lag <- as.integer(args$max_lag %||% "10")
sample_per_group <- as.integer(args$sample_per_group %||% "6000")
seed <- as.integer(args$seed %||% "20260511")
ou_n_starts <- as.integer(args$ou_n_starts %||% "25")
frame_interval <- as.numeric(args$frame_interval %||% "2")
force <- as_flag(args$force, default = FALSE)

if (file.exists(out_rds) && !force) {
  message("Reusing existing alternative-model artifact: ", out_rds)
  quit(save = "no", status = 0)
}
dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)

model_names <- c(
  "Heterogeneous persistent random walk",
  "Stop/go persistent switching",
  "Tethered local motion",
  "Latent PRW + rounded observation"
)
ou_reference_name <- "OU refit comparable"

wrap_angle <- function(theta) {
  ((theta + pi) %% (2 * pi)) - pi
}

clamp <- function(x, lo, hi) {
  pmin(pmax(x, lo), hi)
}

safe_mean <- function(x, default = NA_real_) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) default else mean(x)
}

safe_sd <- function(x, default = 1) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) default else sd(x)
}

safe_quantile <- function(x, p, default = NA_real_) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) default else quantile(x, p, na.rm = TRUE, names = FALSE)
}

fit_lognormal <- function(x, min_value = 1e-4) {
  x <- x[is.finite(x) & x > min_value]
  if (length(x) < 2L) {
    x <- c(min_value, min_value * 1.5)
  }
  lx <- log(pmax(x, min_value))
  list(meanlog = mean(lx), sdlog = max(sd(lx), 0.05))
}

draw_lognormal <- function(params) {
  as.numeric(rlnorm(1L, meanlog = params$meanlog, sdlog = params$sdlog))
}

circular_sd <- function(theta, default = pi / 2) {
  theta <- theta[is.finite(theta)]
  if (length(theta) < 2L) {
    return(default)
  }
  r <- sqrt(mean(cos(theta))^2 + mean(sin(theta))^2)
  clamp(sqrt(pmax(-2 * log(pmax(r, 1e-6)), 0)), 0.03, pi)
}

round_to_grid <- function(x, grid_size) {
  if (!is.finite(grid_size) || grid_size <= 0) {
    return(x)
  }
  round(x / grid_size) * grid_size
}

sample_table <- function(x, sample_per_group, seed) {
  set.seed(seed)
  x |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label) |>
    group_modify(~ {
      n <- min(nrow(.x), sample_per_group)
      if (n == 0L) {
        return(.x)
      }
      .x[sample.int(nrow(.x), n), , drop = FALSE]
    }) |>
    ungroup()
}

quantile_distance <- function(x, y, probs = seq(0.01, 0.99, by = 0.01)) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0L || length(y) == 0L) {
    return(NA_real_)
  }
  mean(abs(
    quantile(x, probs = probs, na.rm = TRUE, names = FALSE) -
      quantile(y, probs = probs, na.rm = TRUE, names = FALSE)
  ))
}

prepare_tracks <- function(tracks_rds, selected_doses) {
  readRDS(tracks_rds) |>
    transmute(
      source = "Observed",
      ploidy = as.character(.data$ploidy),
      Gemcitabine = as.numeric(.data$Gemcitabine),
      dose_label = paste0(.data$Gemcitabine, " nM"),
      site_id = as.character(.data$site_id),
      migration_track_id = paste(.data$site_id, .data$migration_track_id, sep = "::"),
      frame = as.integer(.data$frame),
      x = as.numeric(.data$nucleus_x),
      y = as.numeric(.data$nucleus_y)
    ) |>
    filter(
      .data$Gemcitabine %in% selected_doses,
      is.finite(.data$x),
      is.finite(.data$y),
      is.finite(.data$frame)
    ) |>
    arrange(.data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame)
}

assign_splits <- function(tracks, seed) {
  set.seed(seed)
  site_keys <- tracks |>
    distinct(.data$ploidy, .data$Gemcitabine, .data$dose_label, .data$site_id) |>
    group_by(.data$ploidy, .data$Gemcitabine, .data$dose_label) |>
    group_modify(~ {
      n_site <- nrow(.x)
      if (n_site >= 4L) {
        n_test <- max(1L, floor(0.25 * n_site))
        test_sites <- sample(.x$site_id, n_test)
        .x |> mutate(split = if_else(.data$site_id %in% test_sites, "test", "train"))
      } else {
        .x |> mutate(split = "in_sample")
      }
    }) |>
    ungroup()

  tracks |>
    left_join(site_keys, by = c("ploidy", "Gemcitabine", "dose_label", "site_id"))
}

compute_steps <- function(tracks) {
  tracks |>
    arrange(.data$source, .data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$site_id, .data$migration_track_id) |>
    mutate(
      next_frame = lead(.data$frame),
      dx = lead(.data$x) - .data$x,
      dy = lead(.data$y) - .data$y,
      step_frames = .data$next_frame - .data$frame,
      step_index = row_number(),
      speed = sqrt(.data$dx^2 + .data$dy^2),
      angle = atan2(.data$dy, .data$dx)
    ) |>
    ungroup() |>
    filter(.data$step_frames == 1L, is.finite(.data$speed))
}

compute_cosine_pairs <- function(steps, max_lag) {
  if (nrow(steps) == 0L) {
    return(tibble())
  }
  keys <- c("source", "ploidy", "Gemcitabine", "dose_label", "site_id", "migration_track_id")
  bind_rows(lapply(seq_len(max_lag), function(lag_frame) {
    endpoint <- steps |>
      transmute(
        across(all_of(keys)),
        frame = .data$frame - lag_frame,
        dx_end = .data$dx,
        dy_end = .data$dy,
        speed_end = .data$speed
      )
    steps |>
      select(all_of(keys), "frame", "dx", "dy", "speed") |>
      left_join(endpoint, by = c(keys, "frame")) |>
      filter(is.finite(.data$dx_end), is.finite(.data$dy_end), is.finite(.data$speed_end)) |>
      transmute(
        across(all_of(keys)),
        frame = .data$frame,
        lag_frames = lag_frame,
        step_distance_px = .data$speed,
        lagged_step_distance_px = .data$speed_end,
        cosine_autocorrelation = (.data$dx * .data$dx_end + .data$dy * .data$dy_end) /
          pmax(.data$speed * .data$speed_end, 1e-9)
      )
  }))
}

compute_lag_pairs <- function(tracks, max_lag) {
  keys <- c("source", "ploidy", "Gemcitabine", "dose_label", "site_id", "migration_track_id")
  base <- tracks |>
    select(all_of(keys), "frame", "x", "y")
  bind_rows(lapply(seq_len(max_lag), function(lag_frame) {
    endpoint <- base |>
      transmute(
        across(all_of(keys)),
        frame = .data$frame - lag_frame,
        x_end = .data$x,
        y_end = .data$y
      )
    base |>
      left_join(endpoint, by = c(keys, "frame")) |>
      filter(is.finite(.data$x_end), is.finite(.data$y_end)) |>
      transmute(
        across(all_of(keys)),
        frame = .data$frame,
        lag_frames = lag_frame,
        displacement_px = sqrt((.data$x_end - .data$x)^2 + (.data$y_end - .data$y)^2),
        squared_displacement_px2 = .data$displacement_px^2
      )
  }))
}

summarize_msd <- function(lag_pairs) {
  lag_pairs |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$lag_frames) |>
    summarize(
      mean_squared_displacement_px2 = mean(.data$squared_displacement_px2, na.rm = TRUE),
      median_squared_displacement_px2 = median(.data$squared_displacement_px2, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    )
}

consecutive_turns <- function(steps) {
  if (nrow(steps) == 0L) {
    return(tibble())
  }
  keys <- c("source", "ploidy", "Gemcitabine", "dose_label", "site_id", "migration_track_id")
  endpoint <- steps |>
    transmute(
      across(all_of(keys)),
      frame = .data$frame - 1L,
      dx_next = .data$dx,
      dy_next = .data$dy,
      speed_next = .data$speed,
      angle_next = .data$angle
    )
  steps |>
    select(all_of(keys), "frame", "dx", "dy", "speed", "angle") |>
    left_join(endpoint, by = c(keys, "frame")) |>
    filter(is.finite(.data$dx_next), is.finite(.data$dy_next), is.finite(.data$speed_next)) |>
    mutate(
      turn_angle = wrap_angle(.data$angle_next - .data$angle),
      cosine_autocorrelation = (.data$dx * .data$dx_next + .data$dy * .data$dy_next) /
        pmax(.data$speed * .data$speed_next, 1e-9)
    )
}

fit_heterogeneous_prw <- function(steps, turns) {
  track_stats <- steps |>
    group_by(.data$migration_track_id) |>
    summarize(median_speed = median(.data$speed, na.rm = TRUE), .groups = "drop")
  if (nrow(track_stats) >= 6L && length(unique(round(track_stats$median_speed, 4))) >= 2L) {
    km <- kmeans(log1p(track_stats$median_speed), centers = 2L, nstart = 10)
    center_order <- order(as.numeric(km$centers))
    track_stats$motility_class <- if_else(km$cluster == center_order[[2L]], "fast", "slow")
  } else {
    track_stats$motility_class <- "slow"
  }
  steps2 <- steps |> left_join(track_stats |> select("migration_track_id", "motility_class"), by = "migration_track_id")
  turns2 <- turns |> left_join(track_stats |> select("migration_track_id", "motility_class"), by = "migration_track_id")
  classes <- c("slow", "fast")
  class_params <- lapply(classes, function(cls) {
    cls_steps <- steps2 |> filter(.data$motility_class == cls)
    cls_turns <- turns2 |> filter(.data$motility_class == cls, .data$speed > 0.5, .data$speed_next > 0.5)
    list(
      class = cls,
      speed = fit_lognormal(cls_steps$speed, min_value = 0.05),
      turn_sd = circular_sd(cls_turns$turn_angle, default = circular_sd(turns$turn_angle)),
      n_tracks = sum(track_stats$motility_class == cls),
      median_speed = safe_mean(cls_steps$speed, default = safe_mean(steps$speed, default = 0.1))
    )
  })
  names(class_params) <- classes
  p_fast <- mean(track_stats$motility_class == "fast")
  if (!is.finite(p_fast)) {
    p_fast <- 0
  }
  list(p_fast = p_fast, classes = class_params)
}

fit_stop_go <- function(steps, turns, pause_threshold = 0.5) {
  state_steps <- steps |>
    mutate(step_state = if_else(.data$speed <= pause_threshold, "pause", "move"))
  state_turns <- consecutive_turns(state_steps) |>
    mutate(
      current_state = if_else(.data$speed <= pause_threshold, "pause", "move"),
      next_state = if_else(.data$speed_next <= pause_threshold, "pause", "move")
    )
  transition_probs <- state_turns |>
    count(.data$current_state, .data$next_state, name = "n") |>
    complete(current_state = c("pause", "move"), next_state = c("pause", "move"), fill = list(n = 0L)) |>
    mutate(n = .data$n + 1L) |>
    group_by(.data$current_state) |>
    mutate(probability = .data$n / sum(.data$n)) |>
    ungroup()
  move_turns <- turns |> filter(.data$speed > pause_threshold, .data$speed_next > pause_threshold)
  list(
    pause_threshold = pause_threshold,
    p_initial_pause = mean(state_steps$step_state == "pause", na.rm = TRUE),
    transition_probs = transition_probs,
    pause_step_sd = max(sqrt(mean(c(state_steps$dx[state_steps$step_state == "pause"], state_steps$dy[state_steps$step_state == "pause"])^2, na.rm = TRUE)), 0.02),
    move_speed = fit_lognormal(state_steps$speed[state_steps$step_state == "move"], min_value = pause_threshold),
    move_turn_sd = circular_sd(move_turns$turn_angle, default = circular_sd(turns$turn_angle))
  )
}

fit_tethered <- function(tracks) {
  centered <- tracks |>
    group_by(.data$migration_track_id) |>
    mutate(center_x = mean(.data$x, na.rm = TRUE), center_y = mean(.data$y, na.rm = TRUE)) |>
    ungroup()
  base <- centered |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$site_id, .data$migration_track_id) |>
    arrange(.data$frame, .by_group = TRUE) |>
    mutate(
      next_frame = lead(.data$frame),
      x0 = .data$x - .data$center_x,
      y0 = .data$y - .data$center_y,
      x1 = lead(.data$x) - .data$center_x,
      y1 = lead(.data$y) - .data$center_y,
      step_frames = .data$next_frame - .data$frame
    ) |>
    ungroup() |>
    filter(.data$step_frames == 1L, is.finite(.data$x1), is.finite(.data$y1))
  denom <- sum(base$x0^2 + base$y0^2, na.rm = TRUE)
  rho <- if (is.finite(denom) && denom > 0) {
    sum(base$x0 * base$x1 + base$y0 * base$y1, na.rm = TRUE) / denom
  } else {
    0.5
  }
  rho <- clamp(rho, 0, 0.98)
  residual <- c(base$x1 - rho * base$x0, base$y1 - rho * base$y0)
  innovation_sd <- max(sqrt(mean(residual^2, na.rm = TRUE)), 0.03)
  list(
    rho = rho,
    innovation_sd = innovation_sd,
    stationary_radius = innovation_sd / sqrt(pmax(1 - rho^2, 1e-6))
  )
}

coordinate_grid_size <- function(tracks) {
  z <- c(tracks$x, tracks$y)
  z <- z[is.finite(z)]
  candidates <- c(1, 0.5, 0.25)
  errors <- vapply(candidates, function(g) mean(abs(z / g - round(z / g)) * g, na.rm = TRUE), numeric(1))
  candidates[[which.min(errors)]]
}

fit_observation_prw <- function(steps, turns, tracks, pause_threshold = 0.5) {
  grid_size <- coordinate_grid_size(tracks)
  latent_turns <- turns |> filter(.data$speed > 1, .data$speed_next > 1)
  tiny_fraction <- mean(steps$speed <= pause_threshold, na.rm = TRUE)
  list(
    latent_speed = fit_lognormal(steps$speed[steps$speed > 1], min_value = 0.25),
    latent_turn_sd = circular_sd(latent_turns$turn_angle, default = circular_sd(turns$turn_angle)),
    grid_size = grid_size,
    obs_noise_sd = clamp(0.25 + 0.75 * tiny_fraction, 0.05, 1.25)
  )
}

tracks_to_ou_segments <- function(tracks, frame_interval = 2, min_segment_frames = 3L) {
  ou_tracks <- tracks |>
    transmute(
      split_track_id = as.character(.data$migration_track_id),
      frame = as.numeric(.data$frame),
      x = as.numeric(.data$x),
      y = as.numeric(.data$y)
    ) |>
    arrange(.data$split_track_id, .data$frame) |>
    group_by(.data$split_track_id) |>
    mutate(
      starts_segment = row_number() == 1L | .data$frame != lag(.data$frame) + 1,
      post_filter_segment = cumsum(.data$starts_segment),
      split_track_id = paste(.data$split_track_id, .data$post_filter_segment, sep = "::")
    ) |>
    ungroup() |>
    group_by(.data$split_track_id) |>
    filter(n_distinct(.data$frame) >= min_segment_frames) |>
    ungroup() |>
    select("split_track_id", "frame", "x", "y")
  split_ou_segments(ou_tracks, frame_interval = frame_interval)
}

fit_ou_group <- function(group_tracks, frame_interval, ou_n_starts, seed) {
  segments <- tracks_to_ou_segments(group_tracks, frame_interval = frame_interval)
  if (length(segments) == 0L) {
    return(list(
      fit = tibble(
        tau = NA_real_, velocity_scale = NA_real_, effective_diffusivity = NA_real_,
        obs_noise = NA_real_, log_likelihood = NA_real_, fit_status = "no_valid_segments",
        n_starts = ou_n_starts, seed = seed, start_id = NA_integer_
      ),
      n_segments = 0L,
      total_track_time = 0
    ))
  }
  fit_rows <- fit_ou_velocity(segments, n_starts = ou_n_starts, seed = seed)
  best <- fit_rows |>
    filter(is.finite(.data$log_likelihood)) |>
    arrange(desc(.data$log_likelihood)) |>
    slice(1)
  if (nrow(best) == 0L) {
    best <- fit_rows[1, , drop = FALSE]
  }
  list(
    fit = best,
    n_segments = length(segments),
    total_track_time = attr(segments, "total_track_time") %||% 0
  )
}

simulate_ou_skeleton_track <- function(track_df, ou_fit, frame_interval) {
  track_df <- track_df |> arrange(.data$frame)
  frames <- track_df$frame
  n <- length(frames)
  x <- numeric(n)
  y <- numeric(n)
  x[[1L]] <- track_df$x[[1L]]
  y[[1L]] <- track_df$y[[1L]]

  if (n >= 2L && is.finite(ou_fit$tau[[1]]) && is.finite(ou_fit$velocity_scale[[1]]) && is.finite(ou_fit$obs_noise[[1]])) {
    state_x <- c(0, rnorm(1L, sd = ou_fit$velocity_scale[[1]]))
    state_y <- c(0, rnorm(1L, sd = ou_fit$velocity_scale[[1]]))
    obs_x_prev <- state_x[[1L]] + rnorm(1L, sd = ou_fit$obs_noise[[1]])
    obs_y_prev <- state_y[[1L]] + rnorm(1L, sd = ou_fit$obs_noise[[1]])
    x[[1L]] <- track_df$x[[1L]]
    y[[1L]] <- track_df$y[[1L]]
    for (i in 2:n) {
      dt_i <- max(1e-6, (frames[[i]] - frames[[i - 1L]]) * frame_interval)
      transition <- ou_transition_r(dt = dt_i, tau = ou_fit$tau[[1]], velocity_scale = ou_fit$velocity_scale[[1]])
      noise_chol <- chol(transition$Q)
      state_x <- drop(transition$F %*% state_x + drop(t(noise_chol) %*% rnorm(2L)))
      state_y <- drop(transition$F %*% state_y + drop(t(noise_chol) %*% rnorm(2L)))
      obs_x <- state_x[[1L]] + rnorm(1L, sd = ou_fit$obs_noise[[1]])
      obs_y <- state_y[[1L]] + rnorm(1L, sd = ou_fit$obs_noise[[1]])
      x[[i]] <- x[[i - 1L]] + (obs_x - obs_x_prev)
      y[[i]] <- y[[i - 1L]] + (obs_y - obs_y_prev)
      obs_x_prev <- obs_x
      obs_y_prev <- obs_y
    }
  }

  track_df |> mutate(source = ou_reference_name, x = .env$x, y = .env$y)
}

simulate_ou_group <- function(group_tracks, ou_fit, frame_interval, seed) {
  set.seed(seed)
  track_list <- split(group_tracks, group_tracks$migration_track_id)
  bind_rows(lapply(track_list, simulate_ou_skeleton_track, ou_fit = ou_fit, frame_interval = frame_interval))
}

fit_group_models <- function(group_tracks) {
  steps <- compute_steps(group_tracks)
  turns <- consecutive_turns(steps)
  list(
    heterogeneous_prw = fit_heterogeneous_prw(steps, turns),
    stop_go = fit_stop_go(steps, turns),
    tethered = fit_tethered(group_tracks),
    observation_prw = fit_observation_prw(steps, turns, group_tracks),
    summaries = list(
      n_tracks = n_distinct(group_tracks$migration_track_id),
      n_sites = n_distinct(group_tracks$site_id),
      n_points = nrow(group_tracks),
      n_steps = nrow(steps),
      n_turns = nrow(turns),
      pause_fraction = mean(steps$speed <= 0.5, na.rm = TRUE),
      median_speed = median(steps$speed, na.rm = TRUE),
      p95_speed = safe_quantile(steps$speed, 0.95),
      mean_lag1_cosine = mean(turns$cosine_autocorrelation, na.rm = TRUE)
    )
  )
}

transition_probability <- function(transition_df, current_state) {
  row <- transition_df |> filter(.data$current_state == !!current_state)
  p_move <- row$probability[row$next_state == "move"]
  if (length(p_move) == 0L || !is.finite(p_move)) 0.5 else p_move[[1L]]
}

draw_one_increment <- function(model_name, fit, state) {
  if (model_name == "Heterogeneous persistent random walk") {
    if (is.null(state$motility_class)) {
      state$motility_class <- if (runif(1L) < fit$heterogeneous_prw$p_fast) "fast" else "slow"
    }
    params <- fit$heterogeneous_prw$classes[[state$motility_class]]
    speed <- draw_lognormal(params$speed)
    if (!is.finite(state$heading)) {
      state$heading <- runif(1L, -pi, pi)
    } else {
      state$heading <- wrap_angle(state$heading + rnorm(1L, 0, params$turn_sd))
    }
    return(list(dx = speed * cos(state$heading), dy = speed * sin(state$heading), state = state))
  }

  if (model_name == "Stop/go persistent switching") {
    if (is.null(state$move_state)) {
      state$move_state <- if (runif(1L) < fit$stop_go$p_initial_pause) "pause" else "move"
    } else {
      p_move <- transition_probability(fit$stop_go$transition_probs, state$move_state)
      state$move_state <- if (runif(1L) < p_move) "move" else "pause"
    }
    if (state$move_state == "pause") {
      return(list(dx = rnorm(1L, 0, fit$stop_go$pause_step_sd), dy = rnorm(1L, 0, fit$stop_go$pause_step_sd), state = state))
    }
    speed <- draw_lognormal(fit$stop_go$move_speed)
    if (!is.finite(state$heading)) {
      state$heading <- runif(1L, -pi, pi)
    } else {
      state$heading <- wrap_angle(state$heading + rnorm(1L, 0, fit$stop_go$move_turn_sd))
    }
    return(list(dx = speed * cos(state$heading), dy = speed * sin(state$heading), state = state))
  }

  if (model_name == "Latent PRW + rounded observation") {
    speed <- draw_lognormal(fit$observation_prw$latent_speed)
    if (!is.finite(state$heading)) {
      state$heading <- runif(1L, -pi, pi)
    } else {
      state$heading <- wrap_angle(state$heading + rnorm(1L, 0, fit$observation_prw$latent_turn_sd))
    }
    return(list(dx = speed * cos(state$heading), dy = speed * sin(state$heading), state = state))
  }

  stop("Unknown increment model: ", model_name, call. = FALSE)
}

simulate_tethered_track <- function(track_df, fit) {
  track_df <- track_df |> arrange(.data$frame)
  frames <- track_df$frame
  n <- length(frames)
  x <- numeric(n)
  y <- numeric(n)
  x[[1L]] <- track_df$x[[1L]]
  y[[1L]] <- track_df$y[[1L]]
  center_x <- x[[1L]]
  center_y <- y[[1L]]
  if (n >= 2L) {
    for (i in 2:n) {
      x_current <- x[[i - 1L]]
      y_current <- y[[i - 1L]]
      step_count <- max(1L, as.integer(frames[[i]] - frames[[i - 1L]]))
      for (step_i in seq_len(step_count)) {
        x_current <- center_x + fit$tethered$rho * (x_current - center_x) + rnorm(1L, 0, fit$tethered$innovation_sd)
        y_current <- center_y + fit$tethered$rho * (y_current - center_y) + rnorm(1L, 0, fit$tethered$innovation_sd)
      }
      x[[i]] <- x_current
      y[[i]] <- y_current
    }
  }
  track_df |> mutate(source = "Tethered local motion", x = .env$x, y = .env$y)
}

simulate_one_track <- function(track_df, model_name, fit) {
  if (model_name == "Tethered local motion") {
    return(simulate_tethered_track(track_df, fit))
  }
  track_df <- track_df |> arrange(.data$frame)
  frames <- track_df$frame
  n <- length(frames)
  x <- numeric(n)
  y <- numeric(n)
  latent_x <- numeric(n)
  latent_y <- numeric(n)
  x[[1L]] <- track_df$x[[1L]]
  y[[1L]] <- track_df$y[[1L]]
  latent_x[[1L]] <- x[[1L]]
  latent_y[[1L]] <- y[[1L]]
  state <- list(heading = NA_real_, move_state = NULL, motility_class = NULL)
  if (n >= 2L) {
    for (i in 2:n) {
      x_current <- x[[i - 1L]]
      y_current <- y[[i - 1L]]
      lx_current <- latent_x[[i - 1L]]
      ly_current <- latent_y[[i - 1L]]
      step_count <- max(1L, as.integer(frames[[i]] - frames[[i - 1L]]))
      for (step_i in seq_len(step_count)) {
        step <- draw_one_increment(model_name, fit, state)
        if (model_name == "Latent PRW + rounded observation") {
          lx_current <- lx_current + step$dx
          ly_current <- ly_current + step$dy
          x_current <- round_to_grid(lx_current + rnorm(1L, 0, fit$observation_prw$obs_noise_sd), fit$observation_prw$grid_size)
          y_current <- round_to_grid(ly_current + rnorm(1L, 0, fit$observation_prw$obs_noise_sd), fit$observation_prw$grid_size)
        } else {
          x_current <- x_current + step$dx
          y_current <- y_current + step$dy
          lx_current <- x_current
          ly_current <- y_current
        }
        state <- step$state
      }
      x[[i]] <- x_current
      y[[i]] <- y_current
      latent_x[[i]] <- lx_current
      latent_y[[i]] <- ly_current
    }
  }
  track_df |> mutate(source = model_name, x = .env$x, y = .env$y)
}

simulate_group <- function(group_tracks, model_name, fit, seed) {
  set.seed(seed)
  track_list <- split(group_tracks, group_tracks$migration_track_id)
  bind_rows(lapply(track_list, simulate_one_track, model_name = model_name, fit = fit))
}

round_modeled_observations <- function(tracks) {
  tracks |>
    mutate(
      x = if_else(as.character(.data$source) == "Observed", .data$x, round(.data$x)),
      y = if_else(as.character(.data$source) == "Observed", .data$y, round(.data$y))
    )
}

fit_summary_rows <- function(fit, group_keys, split_label) {
  base_row <- function(model, parameter, value) {
    tibble(
      ploidy = group_keys$ploidy,
      Gemcitabine = group_keys$Gemcitabine,
      dose_label = group_keys$dose_label,
      split_label = split_label,
      model = model,
      parameter = parameter,
      value = as.numeric(value)
    )
  }
  hetero <- bind_rows(
    base_row("Heterogeneous persistent random walk", "p_fast_track_class", fit$heterogeneous_prw$p_fast),
    bind_rows(lapply(names(fit$heterogeneous_prw$classes), function(cls) {
      params <- fit$heterogeneous_prw$classes[[cls]]
      base_row(
        "Heterogeneous persistent random walk",
        paste0(cls, "_", c("speed_meanlog", "speed_sdlog", "turn_sd", "n_training_tracks")),
        c(params$speed$meanlog, params$speed$sdlog, params$turn_sd, params$n_tracks)
      )
    }))
  )
  transition_probs <- fit$stop_go$transition_probs |>
    transmute(
      ploidy = group_keys$ploidy,
      Gemcitabine = group_keys$Gemcitabine,
      dose_label = group_keys$dose_label,
      split_label = split_label,
      model = "Stop/go persistent switching",
      parameter = paste0("p_", .data$current_state, "_to_", .data$next_state),
      value = .data$probability
    )
  bind_rows(
    hetero,
    base_row(
      "Stop/go persistent switching",
      c("pause_threshold_px", "p_initial_pause", "pause_step_sd", "move_speed_meanlog", "move_speed_sdlog", "move_turn_sd"),
      c(
        fit$stop_go$pause_threshold, fit$stop_go$p_initial_pause, fit$stop_go$pause_step_sd,
        fit$stop_go$move_speed$meanlog, fit$stop_go$move_speed$sdlog, fit$stop_go$move_turn_sd
      )
    ),
    transition_probs,
    base_row(
      "Tethered local motion",
      c("rho_to_local_center", "innovation_sd_px", "stationary_radius_px"),
      c(fit$tethered$rho, fit$tethered$innovation_sd, fit$tethered$stationary_radius)
    ),
    base_row(
      "Latent PRW + rounded observation",
      c("latent_speed_meanlog", "latent_speed_sdlog", "latent_turn_sd", "grid_size_px", "obs_noise_sd_px"),
      c(
        fit$observation_prw$latent_speed$meanlog, fit$observation_prw$latent_speed$sdlog,
        fit$observation_prw$latent_turn_sd, fit$observation_prw$grid_size, fit$observation_prw$obs_noise_sd
      )
    )
  )
}

ou_fit_summary_rows <- function(ou_result, group_keys, split_label) {
  fit <- ou_result$fit
  tibble(
    ploidy = group_keys$ploidy,
    Gemcitabine = group_keys$Gemcitabine,
    dose_label = group_keys$dose_label,
    split_label = split_label,
    model = ou_reference_name,
    parameter = c(
      "tau", "velocity_scale", "effective_diffusivity", "obs_noise",
      "log_likelihood", "n_segments", "total_track_time", "n_starts", "best_start_id"
    ),
    value = c(
      fit$tau[[1]], fit$velocity_scale[[1]], fit$effective_diffusivity[[1]], fit$obs_noise[[1]],
      fit$log_likelihood[[1]], ou_result$n_segments, ou_result$total_track_time,
      fit$n_starts[[1]], fit$start_id[[1]]
    )
  )
}

peak_mass_summary <- function(cosine_pairs) {
  cosine_pairs |>
    filter(.data$lag_frames == 1L) |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label) |>
    summarize(
      n_pairs = n(),
      exact_reverse = mean(abs(.data$cosine_autocorrelation + 1) < 1e-12, na.rm = TRUE),
      exact_orthogonal = mean(abs(.data$cosine_autocorrelation) < 1e-12, na.rm = TRUE),
      exact_forward = mean(abs(.data$cosine_autocorrelation - 1) < 1e-12, na.rm = TRUE),
      near_reverse = mean(.data$cosine_autocorrelation <= -0.95, na.rm = TRUE),
      near_orthogonal = mean(abs(.data$cosine_autocorrelation) <= 0.05, na.rm = TRUE),
      near_forward = mean(.data$cosine_autocorrelation >= 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}

observed_direction_step_summary <- function(cosine_pairs) {
  cosine_pairs |>
    filter(
      .data$source == "Observed",
      .data$lag_frames == 1L,
      is.finite(.data$cosine_autocorrelation)
    ) |>
    mutate(
      direction_class = case_when(
        .data$cosine_autocorrelation <= -0.95 ~ "reverse, cosine near -1",
        abs(.data$cosine_autocorrelation) <= 0.05 ~ "orthogonal, cosine near 0",
        .data$cosine_autocorrelation >= 0.95 ~ "forward, cosine near 1",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(.data$direction_class)) |>
    select(
      "source", "ploidy", "Gemcitabine", "dose_label", "site_id", "migration_track_id",
      "frame", "direction_class", "cosine_autocorrelation",
      "step_distance_px", "lagged_step_distance_px"
    ) |>
    pivot_longer(
      cols = c("step_distance_px", "lagged_step_distance_px"),
      names_to = "step_in_pair",
      values_to = "pixel_step_distance"
    ) |>
    mutate(
      step_in_pair = recode(
        .data$step_in_pair,
        step_distance_px = "first step in lag-1 pair",
        lagged_step_distance_px = "second step in lag-1 pair"
      ),
      direction_class = factor(
        .data$direction_class,
        levels = c("reverse, cosine near -1", "orthogonal, cosine near 0", "forward, cosine near 1")
      )
    )
}

score_models <- function(speed_distribution, cosine_distribution, lag_displacement_distribution, msd_summary) {
  model_sources <- setdiff(as.character(unique(speed_distribution$source)), "Observed")
  keys <- speed_distribution |> distinct(.data$ploidy, .data$Gemcitabine, .data$dose_label)
  by_group <- bind_rows(lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    obs_speed <- speed_distribution |>
      semi_join(key, by = c("ploidy", "Gemcitabine", "dose_label")) |>
      filter(.data$source == "Observed") |>
      pull(.data$speed)
    obs_cos <- cosine_distribution |>
      semi_join(key, by = c("ploidy", "Gemcitabine", "dose_label")) |>
      filter(.data$source == "Observed") |>
      pull(.data$cosine_autocorrelation)
    obs_msd <- msd_summary |>
      semi_join(key, by = c("ploidy", "Gemcitabine", "dose_label")) |>
      filter(.data$source == "Observed") |>
      transmute(lag_frames = .data$lag_frames, observed_msd = .data$mean_squared_displacement_px2)
    obs_disp <- lag_displacement_distribution |>
      semi_join(key, by = c("ploidy", "Gemcitabine", "dose_label")) |>
      filter(.data$source == "Observed")
    bind_rows(lapply(model_sources, function(model_source) {
      sim_speed <- speed_distribution |>
        semi_join(key, by = c("ploidy", "Gemcitabine", "dose_label")) |>
        filter(.data$source == model_source) |>
        pull(.data$speed)
      sim_cos <- cosine_distribution |>
        semi_join(key, by = c("ploidy", "Gemcitabine", "dose_label")) |>
        filter(.data$source == model_source) |>
        pull(.data$cosine_autocorrelation)
      sim_msd <- msd_summary |>
        semi_join(key, by = c("ploidy", "Gemcitabine", "dose_label")) |>
        filter(.data$source == model_source) |>
        transmute(lag_frames = .data$lag_frames, simulated_msd = .data$mean_squared_displacement_px2)
      msd_join <- obs_msd |> inner_join(sim_msd, by = "lag_frames")
      sim_disp <- lag_displacement_distribution |>
        semi_join(key, by = c("ploidy", "Gemcitabine", "dose_label")) |>
        filter(.data$source == model_source)
      lagged_disp_distance <- bind_rows(lapply(sort(unique(obs_disp$lag_frames)), function(lag_frame) {
        tibble(
          lag_frames = lag_frame,
          distance = quantile_distance(
            log1p(obs_disp$squared_displacement_px2[obs_disp$lag_frames == lag_frame]),
            log1p(sim_disp$squared_displacement_px2[sim_disp$lag_frames == lag_frame])
          )
        )
      }))
      tibble(
        ploidy = key$ploidy,
        Gemcitabine = key$Gemcitabine,
        dose_label = key$dose_label,
        source = model_source,
        speed_quantile_distance = quantile_distance(obs_speed, sim_speed),
        lag1_cosine_quantile_distance = quantile_distance(obs_cos, sim_cos),
        lagged_log_displacement_quantile_distance = mean(lagged_disp_distance$distance, na.rm = TRUE),
        msd_curve_log_rmse = sqrt(mean((log1p(msd_join$simulated_msd) - log1p(msd_join$observed_msd))^2, na.rm = TRUE))
      )
    }))
  }))
  aggregate <- by_group |>
    group_by(.data$source) |>
    summarize(
      mean_speed_quantile_distance = mean(.data$speed_quantile_distance, na.rm = TRUE),
      mean_lag1_cosine_quantile_distance = mean(.data$lag1_cosine_quantile_distance, na.rm = TRUE),
      mean_lagged_log_displacement_quantile_distance = mean(.data$lagged_log_displacement_quantile_distance, na.rm = TRUE),
      mean_msd_curve_log_rmse = mean(.data$msd_curve_log_rmse, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      composite_rank_score = rank(.data$mean_speed_quantile_distance, ties.method = "average") +
        rank(.data$mean_lag1_cosine_quantile_distance, ties.method = "average") +
        rank(.data$mean_lagged_log_displacement_quantile_distance, ties.method = "average")
    ) |>
    arrange(.data$composite_rank_score)
  list(by_group = by_group, aggregate = aggregate)
}

message("Loading tracks: ", tracks_rds)
observed_tracks_all <- prepare_tracks(tracks_rds, selected_doses) |>
  assign_splits(seed = seed + 11L)

group_keys <- observed_tracks_all |>
  distinct(.data$ploidy, .data$Gemcitabine, .data$dose_label) |>
  arrange(.data$ploidy, .data$Gemcitabine)

split_summary <- observed_tracks_all |>
  group_by(.data$ploidy, .data$Gemcitabine, .data$dose_label, .data$split) |>
  summarize(
    n_sites = n_distinct(.data$site_id),
    n_tracks = n_distinct(.data$migration_track_id),
    n_points = n(),
    .groups = "drop"
  )

fit_tables <- list()
simulated_groups <- list()
observed_eval_groups <- list()
ou_fit_tables <- list()

for (i in seq_len(nrow(group_keys))) {
  key <- group_keys[i, , drop = FALSE]
  group_tracks_all <- observed_tracks_all |>
    filter(.data$ploidy == key$ploidy, .data$Gemcitabine == key$Gemcitabine)
  has_holdout <- any(group_tracks_all$split == "test")
  fit_tracks <- if (has_holdout) {
    group_tracks_all |> filter(.data$split == "train")
  } else {
    group_tracks_all
  }
  eval_tracks <- if (has_holdout) {
    group_tracks_all |> filter(.data$split == "test")
  } else {
    group_tracks_all
  }
  split_label <- if (has_holdout) "site_holdout" else "in_sample"
  message("Fitting condition: ", key$ploidy, " ", key$dose_label, " (", split_label, ")")
  fit <- fit_group_models(fit_tracks)
  fit_tables[[length(fit_tables) + 1L]] <- fit_summary_rows(fit, key, split_label)
  message("Refitting OU for condition: ", key$ploidy, " ", key$dose_label)
  ou_result <- fit_ou_group(
    fit_tracks,
    frame_interval = frame_interval,
    ou_n_starts = ou_n_starts,
    seed = seed + 5000L + i
  )
  ou_fit_tables[[length(ou_fit_tables) + 1L]] <- ou_fit_summary_rows(ou_result, key, split_label)
  observed_eval_groups[[length(observed_eval_groups) + 1L]] <- eval_tracks |> mutate(source = "Observed")

  for (j in seq_along(model_names)) {
    model_name <- model_names[[j]]
    simulated_groups[[length(simulated_groups) + 1L]] <- simulate_group(
      eval_tracks,
      model_name = model_name,
      fit = fit,
      seed = seed + 1000L * i + 37L * j
    )
  }
  simulated_groups[[length(simulated_groups) + 1L]] <- simulate_ou_group(
    eval_tracks,
    ou_fit = ou_result$fit,
    frame_interval = frame_interval,
    seed = seed + 9000L + i
  )
}

observed_eval_tracks <- bind_rows(observed_eval_groups)
comparison_tracks <- bind_rows(observed_eval_tracks, bind_rows(simulated_groups)) |>
  mutate(
    source = factor(.data$source, levels = c("Observed", model_names, ou_reference_name)),
    dose_label = factor(.data$dose_label, levels = paste0(sort(selected_doses), " nM")),
    ploidy = factor(.data$ploidy, levels = sort(unique(.data$ploidy)))
  ) |>
  round_modeled_observations()

message("Computing comparison metrics.")
comparison_steps <- compute_steps(comparison_tracks)
comparison_cosine <- compute_cosine_pairs(comparison_steps, max_lag = max_lag)
comparison_lag_pairs <- compute_lag_pairs(comparison_tracks, max_lag = max_lag)
msd_summary <- summarize_msd(comparison_lag_pairs)

speed_cap <- quantile(comparison_steps$speed, 0.995, na.rm = TRUE, names = FALSE)
disp_cap <- quantile(comparison_lag_pairs$displacement_px, 0.995, na.rm = TRUE, names = FALSE)
speed_distribution <- comparison_steps |>
  mutate(speed_capped = pmin(.data$speed, speed_cap)) |>
  select("source", "ploidy", "Gemcitabine", "dose_label", "speed", "speed_capped") |>
  sample_table(sample_per_group = sample_per_group, seed = seed + 700L)

cosine_distribution <- comparison_cosine |>
  filter(.data$lag_frames == 1L) |>
  select("source", "ploidy", "Gemcitabine", "dose_label", "lag_frames", "cosine_autocorrelation") |>
  sample_table(sample_per_group = sample_per_group, seed = seed + 701L)

lag_displacement_distribution <- comparison_lag_pairs |>
  filter(.data$lag_frames %in% c(1L, 3L, 5L, 10L)) |>
  mutate(displacement_capped = pmin(.data$displacement_px, disp_cap)) |>
  select(
    "source", "ploidy", "Gemcitabine", "dose_label", "lag_frames",
    "displacement_px", "displacement_capped", "squared_displacement_px2"
  ) |>
  sample_table(sample_per_group = sample_per_group, seed = seed + 702L)

fit_parameters <- bind_rows(fit_tables, bind_rows(ou_fit_tables)) |>
  mutate(
    source = .data$model,
    dose_label = factor(.data$dose_label, levels = paste0(sort(selected_doses), " nM")),
    ploidy = factor(.data$ploidy, levels = sort(unique(as.character(.data$ploidy))))
  )

peak_mass <- peak_mass_summary(comparison_cosine)
observed_direction_step_distribution <- observed_direction_step_summary(comparison_cosine)
score_list <- score_models(speed_distribution, cosine_distribution, lag_displacement_distribution, msd_summary)

artifact <- list(
  metadata = list(
    created_at = Sys.time(),
    analysis_dir = analysis_dir,
    tracks_rds = tracks_rds,
    selected_doses = selected_doses,
    max_lag = max_lag,
    sample_per_group = sample_per_group,
    seed = seed,
    frame_interval = frame_interval,
    ou_n_starts = ou_n_starts,
    model_names = model_names,
    reference_model_names = ou_reference_name,
    comparison_design = "Models are fitted on training sites when at least four sites are available per condition; simulations are evaluated on held-out site track skeletons. Smaller groups are labeled in-sample.",
    note = paste(
      "Alternative models are parametric or finite-state generative simulators.",
      "They do not resample observed step vectors or turn angles as the movement model.",
      "All simulated model coordinates are rounded to nearest-pixel observations before predictive checks are computed."
    )
  ),
  split_summary = split_summary,
  fit_parameters = fit_parameters,
  speed_distribution = speed_distribution,
  cosine_distribution = cosine_distribution,
  lag_displacement_distribution = lag_displacement_distribution,
  msd_summary = msd_summary,
  peak_mass = peak_mass,
  observed_direction_step_distribution = observed_direction_step_distribution,
  scores_by_group = score_list$by_group,
  scores = score_list$aggregate,
  track_counts = comparison_tracks |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label) |>
    summarize(n_tracks = n_distinct(.data$migration_track_id), n_points = n(), .groups = "drop")
)

saveRDS(artifact, out_rds)
message("Wrote alternative-model artifact: ", out_rds)
