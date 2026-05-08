#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(parallel)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x)) {
    y
  } else {
    x
  }
}

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("Arguments must use --name=value form: ", arg, call. = FALSE)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    out[[parts[[1L]]]] <- paste(parts[-1L], collapse = "=")
  }
  out
}

timed_step <- function(label, expr) {
  start <- Sys.time()
  message("[", format(start, "%H:%M:%S"), "] Starting ", label)
  result <- force(expr)
  end <- Sys.time()
  message(
    "[", format(end, "%H:%M:%S"), "] Finished ", label,
    " in ", round(as.numeric(difftime(end, start, units = "secs")), 2), " sec"
  )
  result
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

n_tracks <- as.integer(args$n_tracks %||% 350L)
n_frames <- as.integer(args$n_frames %||% 41L)
tau <- as.numeric(args$tau %||% 8)
velocity_scale <- as.numeric(args$velocity_scale %||% 5)
obs_noise <- as.numeric(args$obs_noise %||% 0)
dt <- as.numeric(args$dt %||% 1)
box_width <- as.numeric(args$box_width %||% 512)
box_height <- as.numeric(args$box_height %||% 512)
disc_area <- as.numeric(args$disc_area %||% 450)
max_lag <- as.integer(args$max_lag %||% 20L)
step_size_quantile <- as.numeric(args$step_size_quantile %||% 0.99)
n_observation_boxes <- as.integer(args$n_observation_boxes %||% 16L)
n_cores <- as.integer(args$n_cores %||% min(16L, n_observation_boxes))
seed_base <- as.integer(args$seed_base %||% 20260508L)
output_png <- args$output_png %||% file.path(
  "analyses",
  "K00_GemcitabineExposure_033023",
  "dev",
  "predictive_checks",
  "track_censoring_summary.png"
)

if (n_tracks < 2L) {
  stop("--n_tracks must be at least 2", call. = FALSE)
}
if (n_observation_boxes < 1L || n_cores < 1L) {
  stop("--n_observation_boxes and --n_cores must be at least 1", call. = FALSE)
}
if (n_frames < 2L) {
  stop("--n_frames must be at least 2", call. = FALSE)
}
if (tau <= 0 || velocity_scale <= 0 || obs_noise < 0 || dt <= 0) {
  stop("--tau, --velocity_scale, and --dt must be positive; --obs_noise must be nonnegative", call. = FALSE)
}
if (box_width <= 0 || box_height <= 0 || disc_area <= 0) {
  stop("--box_width, --box_height, and --disc_area must be positive", call. = FALSE)
}
if (step_size_quantile <= 0 || step_size_quantile >= 1) {
  stop("--step_size_quantile must be between 0 and 1", call. = FALSE)
}

disc_radius <- sqrt(disc_area / pi)
if (2 * disc_radius >= min(box_width, box_height)) {
  stop("The disc diameter must be smaller than the observation box.", call. = FALSE)
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

simulate_ou_track <- function(track_id, n_frames, tau, velocity_scale, obs_noise, dt) {
  transition <- ou_transition(dt = dt, tau = tau, velocity_scale = velocity_scale)
  noise_chol <- chol(transition$Q)
  state_x <- c(0, rnorm(1, sd = velocity_scale))
  state_y <- c(0, rnorm(1, sd = velocity_scale))
  rows <- vector("list", n_frames)

  for (i in seq_len(n_frames)) {
    rows[[i]] <- tibble(
      track_id = track_id,
      frame = i - 1L,
      step = i - 1L,
      latent_x = state_x[[1L]],
      latent_y = state_y[[1L]],
      x = state_x[[1L]] + rnorm(1, sd = obs_noise),
      y = state_y[[1L]] + rnorm(1, sd = obs_noise)
    )

    if (i < n_frames) {
      state_x <- drop(transition$F %*% state_x + drop(t(noise_chol) %*% rnorm(2)))
      state_y <- drop(transition$F %*% state_y + drop(t(noise_chol) %*% rnorm(2)))
    }
  }

  bind_rows(rows)
}

simulate_tracks <- function(n_tracks, n_frames, tau, velocity_scale, obs_noise, dt) {
  bind_rows(lapply(seq_len(n_tracks), simulate_ou_track,
    n_frames = n_frames,
    tau = tau,
    velocity_scale = velocity_scale,
    obs_noise = obs_noise,
    dt = dt
  ))
}

place_tracks_in_box <- function(tracks, box_width, box_height, disc_radius) {
  starts <- tibble(
    track_id = sort(unique(tracks$track_id)),
    start_x = runif(n_distinct(tracks$track_id), disc_radius, box_width - disc_radius),
    start_y = runif(n_distinct(tracks$track_id), disc_radius, box_height - disc_radius)
  )

  tracks |>
    left_join(starts, by = "track_id") |>
    group_by(track_id) |>
    mutate(
      x_box = start_x + x - first(x),
      y_box = start_y + y - first(y)
    ) |>
    ungroup()
}

mark_censored_observations <- function(tracks, box_width, box_height, disc_radius) {
  tracks_marked <- tracks |>
    mutate(
      row_id = row_number(),
      boundary_censored = .data$x_box - disc_radius < 0 |
        .data$x_box + disc_radius > box_width |
        .data$y_box - disc_radius < 0 |
        .data$y_box + disc_radius > box_height,
      overlap_censored = FALSE
    )

  overlap_row_ids <- integer()
  for (frame_i in sort(unique(tracks$frame))) {
    frame_tracks <- tracks_marked |>
      filter(.data$frame == frame_i) |>
      arrange(.data$track_id)

    if (nrow(frame_tracks) < 2L) {
      next
    }

    dx <- outer(frame_tracks$x_box, frame_tracks$x_box, "-")
    dy <- outer(frame_tracks$y_box, frame_tracks$y_box, "-")
    overlap <- sqrt(dx^2 + dy^2) < 2 * disc_radius
    overlap[lower.tri(overlap, diag = TRUE)] <- FALSE

    overlap_row_ids <- c(
      overlap_row_ids,
      frame_tracks$row_id[row(overlap)[overlap]],
      frame_tracks$row_id[col(overlap)[overlap]]
    )
  }

  tracks_marked |>
    mutate(
      overlap_censored = .data$row_id %in% unique(overlap_row_ids),
      is_censored = .data$boundary_censored | .data$overlap_censored,
      censor_reason = case_when(
        .data$boundary_censored ~ "outside observation box",
        .data$overlap_censored ~ "disc overlap",
        TRUE ~ "observed"
      )
    ) |>
    mutate(
      censor_reason = factor(.data$censor_reason, levels = c("observed", "outside observation box", "disc overlap"))
    )
}

censor_tracks <- function(tracks) {
  tracks |>
    filter(!.data$is_censored)
}

censor_tracks_by_flag <- function(tracks, flag_col) {
  tracks |>
    filter(!.data[[flag_col]])
}

split_segmented_tracks <- function(censored_tracks) {
  censored_tracks |>
    arrange(.data$track_id, .data$frame) |>
    group_by(.data$track_id) |>
    mutate(
      new_segment = row_number() == 1L | .data$frame != lag(.data$frame) + 1,
      segment_index = cumsum(.data$new_segment),
      segment_id = paste(.data$track_id, .data$segment_index, sep = "_seg")
    ) |>
    ungroup() |>
    select(-new_segment, -segment_index)
}

calculate_step_sizes <- function(tracks, x_col = "x_box", y_col = "y_box") {
  tracks |>
    arrange(.data$segment_id, .data$frame) |>
    group_by(.data$segment_id) |>
    mutate(
      previous_frame = lag(.data$frame),
      step_size = sqrt((.data[[x_col]] - lag(.data[[x_col]]))^2 + (.data[[y_col]] - lag(.data[[y_col]]))^2)
    ) |>
    ungroup() |>
    filter(.data$frame == .data$previous_frame + 1L, !is.na(.data$step_size)) |>
    select(segment_id, track_id, frame, step_size)
}

mark_proactive_neighbor_censoring <- function(tracks, reference_tracks, distance_threshold, x_col = "x_box", y_col = "y_box") {
  if (!is.finite(distance_threshold) || distance_threshold <= 0) {
    stop("distance_threshold must be a positive finite value", call. = FALSE)
  }

  tracks_marked <- tracks |>
    mutate(
      proactive_censored = FALSE,
      nearest_current_object_distance = NA_real_
    )

  frame_results <- lapply(sort(unique(tracks_marked$frame)), function(frame_i) {
    frame_tracks <- tracks_marked |>
      filter(.data$frame == frame_i) |>
      arrange(.data$track_id, .data$segment_id)
    reference_frame_tracks <- reference_tracks |>
      filter(.data$frame == frame_i) |>
      arrange(.data$track_id, .data$row_id)

    if (nrow(frame_tracks) < 1L || nrow(reference_frame_tracks) < 2L) {
      return(frame_tracks)
    }

    dx <- outer(frame_tracks[[x_col]], reference_frame_tracks[[x_col]], "-")
    dy <- outer(frame_tracks[[y_col]], reference_frame_tracks[[y_col]], "-")
    distance <- sqrt(dx^2 + dy^2)
    same_object <- outer(frame_tracks$row_id, reference_frame_tracks$row_id, "==")
    distance[same_object] <- Inf
    nearest_distance <- apply(distance, 1, min)

    frame_tracks |>
      mutate(
        nearest_current_object_distance = nearest_distance,
        proactive_censored = nearest_distance < distance_threshold
      )
  })

  bind_rows(frame_results) |>
    arrange(.data$row_id)
}

proactively_censor_tracks <- function(tracks) {
  tracks |>
    filter(!.data$proactive_censored)
}

summarize_windowed_msd <- function(tracks, source, max_lag, x_col = "x", y_col = "y") {
  starts <- tracks |>
    transmute(
      segment_id = .data$segment_id,
      frame = .data$frame,
      x0 = .data[[x_col]],
      y0 = .data[[y_col]]
    )

  ends <- tracks |>
    transmute(
      segment_id = .data$segment_id,
      frame = .data$frame,
      x1 = .data[[x_col]],
      y1 = .data[[y_col]]
    )

  bind_rows(lapply(seq_len(max_lag), function(lag_frame) {
    starts |>
      inner_join(
        ends |> mutate(frame = .data$frame - lag_frame),
        by = c("segment_id", "frame")
      ) |>
      summarize(
        lag_frames = lag_frame,
        mean_squared_displacement = mean((.data$x1 - .data$x0)^2 + (.data$y1 - .data$y0)^2),
        n_displacements = n(),
        .groups = "drop"
      )
  })) |>
    mutate(source = source)
}

prepare_steps <- function(tracks, source, x_col = "x", y_col = "y") {
  tracks |>
    arrange(.data$segment_id, .data$frame) |>
    group_by(.data$segment_id) |>
    mutate(
      previous_frame = lag(.data$frame),
      dx = .data[[x_col]] - lag(.data[[x_col]]),
      dy = .data[[y_col]] - lag(.data[[y_col]])
    ) |>
    ungroup() |>
    filter(.data$frame == .data$previous_frame + 1L, !is.na(.data$dx), !is.na(.data$dy)) |>
    select(segment_id, frame, dx, dy) |>
    mutate(source = source)
}

summarize_step_acf <- function(steps, max_lag) {
  denominator <- steps |>
    group_by(.data$source) |>
    summarize(
      mean_step_squared = mean(.data$dx^2 + .data$dy^2),
      .groups = "drop"
    )

  starts <- steps |>
    transmute(
      source = .data$source,
      segment_id = .data$segment_id,
      frame = .data$frame,
      dx0 = .data$dx,
      dy0 = .data$dy
    )

  ends <- steps |>
    transmute(
      source = .data$source,
      segment_id = .data$segment_id,
      frame = .data$frame,
      dx1 = .data$dx,
      dy1 = .data$dy
    )

  bind_rows(lapply(seq_len(max_lag), function(lag_frame) {
    starts |>
      inner_join(
        ends |> mutate(frame = .data$frame - lag_frame),
        by = c("source", "segment_id", "frame")
      ) |>
      group_by(.data$source) |>
      summarize(
        lag_frames = lag_frame,
        mean_step_dot_product = mean(.data$dx0 * .data$dx1 + .data$dy0 * .data$dy1),
        n_step_pairs = n(),
        .groups = "drop"
      )
  })) |>
    left_join(denominator, by = "source") |>
    mutate(step_acf = .data$mean_step_dot_product / .data$mean_step_squared)
}

prefix_segment_ids <- function(tracks, observation_box) {
  tracks |>
    mutate(
      observation_box = observation_box,
      segment_id = paste(observation_box, .data$segment_id, sep = "_box")
    )
}

simulate_observation_box <- function(observation_box) {
  set.seed(seed_base + observation_box)

  tracks_uncensored <- simulate_tracks(
    n_tracks = n_tracks,
    n_frames = n_frames,
    tau = tau,
    velocity_scale = velocity_scale,
    obs_noise = obs_noise,
    dt = dt
  ) |>
    mutate(segment_id = as.character(.data$track_id)) |>
    prefix_segment_ids(observation_box = observation_box)

  tracks_in_box <- place_tracks_in_box(
    tracks = tracks_uncensored,
    box_width = box_width,
    box_height = box_height,
    disc_radius = disc_radius
  ) |>
    mark_censored_observations(
      box_width = box_width,
      box_height = box_height,
      disc_radius = disc_radius
    ) |>
    mutate(observation_box = observation_box)

  tracks_censored <- censor_tracks(
    tracks = tracks_in_box
  ) |>
    split_segmented_tracks() |>
    prefix_segment_ids(observation_box = observation_box)

  tracks_boundary_only_censored <- censor_tracks_by_flag(
    tracks = tracks_in_box,
    flag_col = "boundary_censored"
  ) |>
    split_segmented_tracks() |>
    prefix_segment_ids(observation_box = observation_box)

  tracks_overlap_only_censored <- censor_tracks_by_flag(
    tracks = tracks_in_box,
    flag_col = "overlap_censored"
  ) |>
    split_segmented_tracks() |>
    prefix_segment_ids(observation_box = observation_box)

  step_sizes_censored <- calculate_step_sizes(
    tracks_censored,
    x_col = "x_box",
    y_col = "y_box"
  )

  step_size_threshold <- unname(quantile(
    step_sizes_censored$step_size,
    probs = step_size_quantile,
    na.rm = TRUE
  ))
  proactive_distance_threshold <- 2 * disc_radius + step_size_threshold

  proactively_marked_tracks <- tracks_censored |>
    mark_proactive_neighbor_censoring(
      reference_tracks = tracks_in_box,
      distance_threshold = proactive_distance_threshold,
      x_col = "x_box",
      y_col = "y_box"
    )

  tracks_doubly_censored <- proactively_marked_tracks |>
    proactively_censor_tracks() |>
    select(-segment_id) |>
    split_segmented_tracks() |>
    prefix_segment_ids(observation_box = observation_box)

  censor_summary <- tracks_in_box |>
    group_by(.data$censor_reason) |>
    summarize(
      n_observations = n(),
      n_tracks_affected = n_distinct(.data$track_id),
      .groups = "drop"
    ) |>
    mutate(observation_box = observation_box)

  proactive_censor_summary <- proactively_marked_tracks |>
    summarize(
      observation_box = first(.data$observation_box),
      step_size_quantile = step_size_quantile,
      step_size_threshold = step_size_threshold,
      proactive_distance_threshold = proactive_distance_threshold,
      n_observations_input = n(),
      n_observations_proactively_censored = sum(.data$proactive_censored),
      n_tracks_proactively_censored = n_distinct(.data$track_id[.data$proactive_censored])
    )

  list(
    tracks_uncensored = tracks_uncensored,
    tracks_in_box = tracks_in_box,
    tracks_boundary_only_censored = tracks_boundary_only_censored,
    tracks_overlap_only_censored = tracks_overlap_only_censored,
    tracks_censored = tracks_censored,
    tracks_doubly_censored = tracks_doubly_censored,
    censor_summary = censor_summary,
    proactive_censor_summary = proactive_censor_summary
  )
}

message("Simulating ", n_observation_boxes, " observation boxes with ", min(n_cores, n_observation_boxes), " cores")
box_results <- timed_step("parallel observation-box simulation", {
  if (n_cores == 1L) {
    lapply(seq_len(n_observation_boxes), simulate_observation_box)
  } else {
    mclapply(
      seq_len(n_observation_boxes),
      simulate_observation_box,
      mc.cores = min(n_cores, n_observation_boxes)
    )
  }
})

timed_step("binding per-box outputs", {
  tracks_uncensored <<- bind_rows(lapply(box_results, `[[`, "tracks_uncensored"))
  tracks_in_box <<- bind_rows(lapply(box_results, `[[`, "tracks_in_box"))
  tracks_boundary_only_censored <<- bind_rows(lapply(box_results, `[[`, "tracks_boundary_only_censored"))
  tracks_overlap_only_censored <<- bind_rows(lapply(box_results, `[[`, "tracks_overlap_only_censored"))
  tracks_censored <<- bind_rows(lapply(box_results, `[[`, "tracks_censored"))
  tracks_doubly_censored <<- bind_rows(lapply(box_results, `[[`, "tracks_doubly_censored"))
  box_censor_summary <<- bind_rows(lapply(box_results, `[[`, "censor_summary"))
  box_proactive_censor_summary <<- bind_rows(lapply(box_results, `[[`, "proactive_censor_summary"))
})

msd_summary <- timed_step("MSD summaries", {
  bind_rows(
    summarize_windowed_msd(
      tracks_uncensored,
      source = "raw OU",
      max_lag = max_lag,
      x_col = "x",
      y_col = "y"
    ),
    summarize_windowed_msd(
      tracks_boundary_only_censored,
      source = "boundary only",
      max_lag = max_lag,
      x_col = "x_box",
      y_col = "y_box"
    ),
    summarize_windowed_msd(
      tracks_overlap_only_censored,
      source = "overlap only",
      max_lag = max_lag,
      x_col = "x_box",
      y_col = "y_box"
    ),
    summarize_windowed_msd(
      tracks_censored,
      source = "censored",
      max_lag = max_lag,
      x_col = "x_box",
      y_col = "y_box"
    ),
    summarize_windowed_msd(
      tracks_doubly_censored,
      source = "doubly censored",
      max_lag = max_lag,
      x_col = "x_box",
      y_col = "y_box"
    )
  )
})

step_summary <- timed_step("step extraction", {
  bind_rows(
    prepare_steps(tracks_uncensored, source = "raw OU", x_col = "x", y_col = "y"),
    prepare_steps(tracks_boundary_only_censored, source = "boundary only", x_col = "x_box", y_col = "y_box"),
    prepare_steps(tracks_overlap_only_censored, source = "overlap only", x_col = "x_box", y_col = "y_box"),
    prepare_steps(tracks_censored, source = "censored", x_col = "x_box", y_col = "y_box"),
    prepare_steps(tracks_doubly_censored, source = "doubly censored", x_col = "x_box", y_col = "y_box")
  )
})

acf_summary <- timed_step("ACF and lag-support summaries", {
  summarize_step_acf(step_summary, max_lag = max_lag)
})

timed_step("censor summaries", {
  censor_summary <<- box_censor_summary |>
    group_by(.data$censor_reason) |>
    summarize(
      n_observations = sum(.data$n_observations),
      n_tracks_affected = sum(.data$n_tracks_affected),
      .groups = "drop"
    )

  proactive_censor_summary <<- box_proactive_censor_summary |>
    summarize(
      step_size_quantile = first(.data$step_size_quantile),
      step_size_threshold = median(.data$step_size_threshold),
      proactive_distance_threshold = median(.data$proactive_distance_threshold),
      n_observations_input = sum(.data$n_observations_input),
      n_observations_proactively_censored = sum(.data$n_observations_proactively_censored),
      n_tracks_proactively_censored = sum(.data$n_tracks_proactively_censored),
      .groups = "drop"
    )

  class_support_summary <<- bind_rows(
    tibble(source = "raw OU", n_observations = nrow(tracks_uncensored), n_segments = n_distinct(tracks_uncensored$segment_id)),
    tibble(source = "boundary only", n_observations = nrow(tracks_boundary_only_censored), n_segments = n_distinct(tracks_boundary_only_censored$segment_id)),
    tibble(source = "overlap only", n_observations = nrow(tracks_overlap_only_censored), n_segments = n_distinct(tracks_overlap_only_censored$segment_id)),
    tibble(source = "censored", n_observations = nrow(tracks_censored), n_segments = n_distinct(tracks_censored$segment_id)),
    tibble(source = "doubly censored", n_observations = nrow(tracks_doubly_censored), n_segments = n_distinct(tracks_doubly_censored$segment_id))
  )
})

msd_subtitle <- paste(
  paste0(
    n_observation_boxes, " observation boxes; initial censoring: ",
    paste(
      paste0(
        censor_summary$censor_reason,
        "=", censor_summary$n_observations,
        " obs / ", censor_summary$n_tracks_affected, " tracks"
      ),
      collapse = "; "
    )
  ),
  paste0(
    "Proactive threshold: q", round(100 * proactive_censor_summary$step_size_quantile),
    " step=", round(proactive_censor_summary$step_size_threshold, 2), " px; ",
    "distance=", round(proactive_censor_summary$proactive_distance_threshold, 2), " px; ",
    proactive_censor_summary$n_observations_proactively_censored,
    " obs / ", proactive_censor_summary$n_tracks_proactively_censored,
    " tracks"
  ),
  sep = " | "
)

combined_plot <- timed_step("plot construction", {
  trajectory_box <- min(tracks_in_box$observation_box)
  trajectory_tracks_in_box <- tracks_in_box |>
    filter(.data$observation_box == trajectory_box)
  trajectory_tracks_censored <- tracks_censored |>
    filter(.data$observation_box == trajectory_box)
  trajectory_tracks_doubly_censored <- tracks_doubly_censored |>
    filter(.data$observation_box == trajectory_box)

  segment_starts <- tracks_censored |>
    filter(.data$observation_box == trajectory_box) |>
    group_by(.data$segment_id) |>
    slice_min(.data$frame, n = 1, with_ties = FALSE) |>
    ungroup()

  doubly_censored_starts <- tracks_doubly_censored |>
    filter(.data$observation_box == trajectory_box) |>
    group_by(.data$segment_id) |>
    slice_min(.data$frame, n = 1, with_ties = FALSE) |>
    ungroup()

  trajectory_plot <- ggplot() +
    geom_rect(
      aes(xmin = 0, xmax = box_width, ymin = 0, ymax = box_height),
      fill = NA,
      color = "grey25",
      linewidth = 0.5
    ) +
    geom_path(
      data = trajectory_tracks_in_box,
      aes(.data$x_box, .data$y_box, group = .data$track_id),
      color = "grey75",
      linewidth = 0.25,
      alpha = 0.45
    ) +
    geom_path(
      data = trajectory_tracks_censored,
      aes(.data$x_box, .data$y_box, group = .data$segment_id),
      color = "#1b6ca8",
      linewidth = 0.35,
      alpha = 0.75
    ) +
    geom_path(
      data = trajectory_tracks_doubly_censored,
      aes(.data$x_box, .data$y_box, group = .data$segment_id),
      color = "#7b3294",
      linewidth = 0.35,
      alpha = 0.8
    ) +
    geom_point(
      data = segment_starts,
      aes(.data$x_box, .data$y_box),
      shape = 21,
      fill = "#f2c14e",
      color = "grey20",
      size = pmax(0.6, disc_radius / 7),
      stroke = 0.2,
      alpha = 0.55
    ) +
    geom_point(
      data = doubly_censored_starts,
      aes(.data$x_box, .data$y_box),
      shape = 21,
      fill = "#d01c8b",
      color = "grey20",
      size = pmax(0.55, disc_radius / 8),
      stroke = 0.2,
      alpha = 0.55
    ) +
    coord_fixed(xlim = c(-disc_radius, box_width + disc_radius), ylim = c(-disc_radius, box_height + disc_radius)) +
    labs(
      title = "OU trajectories translated into a fixed observation box",
      subtitle = paste0(
        "Observation box ", trajectory_box, " of ", n_observation_boxes, "; ",
        n_tracks, " tracks; disc area = ", round(disc_area, 1),
        " px^2; radius = ", round(disc_radius, 1), " px"
      ),
      x = "x position in observation box",
      y = "y position in observation box"
    ) +
    theme_bw()

  msd_plot <- ggplot(msd_summary, aes(.data$lag_frames * dt, .data$mean_squared_displacement, color = .data$source)) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(size = .data$n_displacements), alpha = 0.85) +
    scale_size_continuous(range = c(1.5, 4.5)) +
    labs(
      title = "Raw, censored, and proactively censored MSD estimates",
      subtitle = msd_subtitle,
      x = "lag",
      y = "mean squared displacement",
      color = NULL,
      size = "displacements"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  acf_plot <- ggplot(acf_summary, aes(.data$lag_frames * dt, .data$step_acf, color = .data$source)) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 0.35) +
    geom_line(linewidth = 0.9) +
    geom_point(alpha = 0.85) +
    labs(
      title = "Step autocorrelation",
      x = "step lag",
      y = "mean dot product / mean step squared",
      color = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  support_plot <- ggplot(acf_summary, aes(.data$lag_frames * dt, .data$n_step_pairs, color = .data$source)) +
    geom_line(linewidth = 0.9) +
    geom_point(alpha = 0.85) +
    scale_y_log10() +
    labs(
      title = "ACF support by lag",
      x = "step lag",
      y = "step-pair count",
      color = NULL
    ) +
    theme_bw() +
    theme(legend.position = "bottom")

  (trajectory_plot / msd_plot / acf_plot / support_plot) +
    plot_layout(heights = c(1.4, 1, 1, 1)) +
    plot_annotation(title = "Track censoring predictive check")
})

print(censor_summary)
print(proactive_censor_summary)
print(class_support_summary)

timed_step("PNG save", {
  dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = output_png,
    plot = combined_plot,
    width = 12,
    height = 16,
    dpi = 180
  )
})
message("Saved combined plot: ", output_png)
