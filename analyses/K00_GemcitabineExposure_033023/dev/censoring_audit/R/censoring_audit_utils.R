suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(as.character(x))) {
    y
  } else {
    x
  }
}

parse_cli_args <- function(args, usage = NULL) {
  if (any(args %in% c("--help", "-h", "--help=TRUE"))) {
    if (!is.null(usage)) {
      cat(usage)
    }
    quit(save = "no", status = 0)
  }

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

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0L) {
    return(normalizePath(".", mustWork = TRUE))
  }
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

audit_dir_from_script <- function() {
  normalizePath(file.path(dirname(script_path()), ".."), mustWork = TRUE)
}

analysis_dir_from_audit_dir <- function(audit_dir) {
  normalizePath(file.path(audit_dir, "../.."), mustWork = TRUE)
}

order_dose_label <- function(label, dose) {
  dose_levels <- sort(unique(dose[is.finite(dose)]))
  factor(label, levels = paste0(dose_levels, " nM"))
}

condition_token <- function(ploidy, dose) {
  paste(ploidy, dose, sep = ":")
}

parse_condition_tokens <- function(x) {
  parts <- strsplit(x, ",", fixed = TRUE)[[1L]]
  rows <- lapply(parts, function(part) {
    fields <- strsplit(part, ":", fixed = TRUE)[[1L]]
    if (length(fields) != 2L) {
      stop("Conditions must use ploidy:dose form, e.g. 2N:0,4N:25", call. = FALSE)
    }
    tibble(ploidy = fields[[1L]], Gemcitabine = as.numeric(fields[[2L]]))
  })
  bind_rows(rows)
}

cap_quantile <- function(x, prob = 0.99) {
  cap <- quantile(x, prob, na.rm = TRUE, names = FALSE)
  pmin(x, cap)
}

time_window_label <- function(frame) {
  case_when(
    frame >= 0 & frame < 10 ~ "frames 0-10",
    frame >= 10 & frame < 20 ~ "frames 10-20",
    frame >= 20 & frame <= 40 ~ "frames 20-40",
    TRUE ~ NA_character_
  )
}

load_site_manifest <- function(analysis_dir) {
  manifest_dir <- file.path(analysis_dir, "cpsam_full_stacks/manifest_rows")
  manifest_files <- list.files(manifest_dir, pattern = "_cpsam_manifest[.]csv$", full.names = TRUE)
  if (length(manifest_files) == 0L) {
    stop("No CPSAM manifest files found in ", manifest_dir, call. = FALSE)
  }

  bind_rows(lapply(manifest_files, read_csv, show_col_types = FALSE)) |>
    distinct(.data$site_id, .data$well, .data$position, .data$n_frames, .data$height, .data$width)
}

load_object_context <- function(analysis_dir) {
  object_distance_dir <- file.path(analysis_dir, "cpsam_full_stacks/nearest_distances/object_tables")
  object_distance_files <- list.files(
    object_distance_dir,
    pattern = "_nearest_object_distances[.]tsv$",
    full.names = TRUE
  )
  if (length(object_distance_files) == 0L) {
    stop("No CPSAM object-distance files found in ", object_distance_dir, call. = FALSE)
  }

  raw <- bind_rows(lapply(object_distance_files, read_tsv, show_col_types = FALSE))

  object_distances <- raw |>
    transmute(
      site_id = .data$site_id,
      frame = .data$frame,
      cpsam_label = .data$cpsam_label,
      nearest_neighbor_distance = .data$nearest_empty_gap_px,
      nearest_mask_center_distance = .data$nearest_mask_center_distance_px,
      nearest_cpsam_label = .data$nearest_cpsam_label,
      n_objects_frame = .data$n_objects_frame,
      distance_definition = .data$distance_definition
    )

  confluency_by_frame <- raw |>
    group_by(.data$site_id, .data$frame) |>
    summarize(
      total_cpsam_area_px = sum(.data$cpsam_area_px, na.rm = TRUE),
      n_cpsam_objects = n(),
      .groups = "drop"
    )

  list(
    object_distances = object_distances,
    confluency_by_frame = confluency_by_frame
  )
}

prepare_observed_at_risk <- function(analysis_dir, tracks_rds = NULL) {
  tracks_rds <- tracks_rds %||% file.path(
    analysis_dir,
    "data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"
  )
  tracks_raw <- readRDS(tracks_rds)
  site_manifest <- load_site_manifest(analysis_dir)
  object_context <- load_object_context(analysis_dir)

  tracks <- tracks_raw |>
    left_join(site_manifest, by = c("site_id", "well", "position")) |>
    left_join(object_context$object_distances, by = c("site_id", "frame", "cpsam_label")) |>
    left_join(object_context$confluency_by_frame, by = c("site_id", "frame")) |>
    mutate(
      x = as.numeric(.data$nucleus_x),
      y = as.numeric(.data$nucleus_y),
      ploidy = factor(.data$ploidy, levels = sort(unique(.data$ploidy))),
      dose_label = paste0(.data$Gemcitabine, " nM"),
      dose_label = order_dose_label(.data$dose_label, .data$Gemcitabine),
      dose_log10 = log10(.data$Gemcitabine + 1),
      distance_to_edge = pmin(
        .data$x,
        .data$width - .data$x,
        .data$y,
        .data$height - .data$y,
        na.rm = FALSE
      ),
      image_area_px = .data$width * .data$height,
      confluency = .data$total_cpsam_area_px / .data$image_area_px,
      confluency_percent = 100 * .data$confluency,
      frame_window = time_window_label(.data$frame)
    ) |>
    arrange(.data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$migration_track_id) |>
    mutate(
      previous_frame = lag(.data$frame),
      previous_x = lag(.data$x),
      previous_y = lag(.data$y),
      step_speed = if_else(
        .data$frame == .data$previous_frame + 1,
        sqrt((.data$x - .data$previous_x)^2 + (.data$y - .data$previous_y)^2),
        NA_real_
      )
    ) |>
    ungroup()

  track_summary <- tracks |>
    group_by(
      .data$migration_track_id, .data$site_id, .data$well, .data$ploidy,
      .data$Gemcitabine, .data$dose_label, .data$position, .data$condition
    ) |>
    summarize(
      first_frame = min(.data$frame, na.rm = TRUE),
      last_frame = max(.data$frame, na.rm = TRUE),
      track_length = n_distinct(.data$frame),
      mean_speed = mean(.data$step_speed, na.rm = TRUE),
      early_track_speed = mean(.data$step_speed[.data$frame <= min(.data$frame, na.rm = TRUE) + 3], na.rm = TRUE),
      site_last_frame = max(.data$n_frames, na.rm = TRUE) - 1L,
      start_x = .data$x[which.min(.data$frame)],
      start_y = .data$y[which.min(.data$frame)],
      width = .data$width[which.min(.data$frame)],
      height = .data$height[which.min(.data$frame)],
      .groups = "drop"
    ) |>
    mutate(
      mean_speed = if_else(is.nan(.data$mean_speed), NA_real_, .data$mean_speed),
      early_track_speed = if_else(is.nan(.data$early_track_speed), NA_real_, .data$early_track_speed)
    )

  at_risk <- tracks |>
    select(all_of(c(
      "migration_track_id", "site_id", "well", "ploidy", "Gemcitabine",
      "dose_label", "position", "condition", "frame", "frame_window",
      "n_frames", "distance_to_edge", "nearest_neighbor_distance", "confluency",
      "confluency_percent", "step_speed", "dose_log10", "x", "y"
    ))) |>
    left_join(
      track_summary |>
        select(all_of(c("migration_track_id", "last_frame", "mean_speed", "early_track_speed"))),
      by = "migration_track_id"
    ) |>
    mutate(
      site_last_frame = .data$n_frames - 1L,
      end_next_frame = as.integer(.data$frame == .data$last_frame & .data$frame < .data$site_last_frame),
      model_speed = coalesce(.data$step_speed, .data$early_track_speed, .data$mean_speed)
    ) |>
    filter(
      .data$frame < .data$site_last_frame,
      is.finite(.data$distance_to_edge),
      is.finite(.data$nearest_neighbor_distance),
      is.finite(.data$confluency),
      is.finite(.data$model_speed),
      !is.na(.data$ploidy)
    )

  context_pool <- at_risk |>
    select(all_of(c(
      "ploidy", "Gemcitabine", "dose_label", "frame",
      "nearest_neighbor_distance", "confluency", "confluency_percent"
    )))

  list(
    tracks = tracks,
    track_summary = track_summary,
    at_risk = at_risk,
    context_pool = context_pool,
    site_manifest = site_manifest
  )
}

best_ou_fits <- function(fits_csv, conditions = NULL) {
  fits <- read_csv(fits_csv, show_col_types = FALSE) |>
    mutate(
      ploidy = as.character(.data$ploidy),
      Gemcitabine = as.numeric(.data$Gemcitabine),
      dose_label = paste0(.data$Gemcitabine, " nM")
    ) |>
    group_by(.data$job_id, .data$ploidy, .data$Gemcitabine) |>
    slice_max(order_by = .data$log_likelihood, n = 1L, with_ties = FALSE) |>
    ungroup()

  if (!is.null(conditions)) {
    fits <- fits |>
      semi_join(conditions, by = c("ploidy", "Gemcitabine"))
  }

  fits |>
    mutate(dose_label = order_dose_label(.data$dose_label, .data$Gemcitabine))
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
  list(F = matrix(c(1, 0, f12, phi), nrow = 2), Q = matrix(c(q11, q12, q12, q22), nrow = 2))
}

simulate_ou_track <- function(tau, velocity_scale, obs_noise, first_frame, last_frame, dt = 1) {
  n_steps <- last_frame - first_frame
  if (n_steps < 1L) {
    stop("Simulated tracks need at least two frames.", call. = FALSE)
  }

  transition <- ou_transition(dt = dt, tau = tau, velocity_scale = velocity_scale)
  noise_chol <- chol(transition$Q)
  state_x <- c(0, rnorm(1, sd = velocity_scale))
  state_y <- c(0, rnorm(1, sd = velocity_scale))
  rows <- vector("list", n_steps + 1L)

  for (i in seq_len(n_steps + 1L)) {
    rows[[i]] <- tibble(
      frame = first_frame + i - 1L,
      track_step = i - 1L,
      rel_x = state_x[[1L]] + rnorm(1, sd = obs_noise),
      rel_y = state_y[[1L]] + rnorm(1, sd = obs_noise)
    )
    if (i <= n_steps) {
      state_x <- drop(transition$F %*% state_x + drop(t(noise_chol) %*% rnorm(2)))
      state_y <- drop(transition$F %*% state_y + drop(t(noise_chol) %*% rnorm(2)))
    }
  }

  bind_rows(rows)
}

add_track_speeds <- function(tracks) {
  tracks |>
    arrange(.data$source, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$source, .data$site_id, .data$migration_track_id) |>
    mutate(
      previous_frame = lag(.data$frame),
      previous_x = lag(.data$x),
      previous_y = lag(.data$y),
      step_speed = if_else(
        .data$frame == .data$previous_frame + 1,
        sqrt((.data$x - .data$previous_x)^2 + (.data$y - .data$previous_y)^2),
        NA_real_
      ),
      mean_speed = mean(.data$step_speed, na.rm = TRUE),
      early_track_speed = mean(.data$step_speed[.data$track_step <= 3], na.rm = TRUE),
      model_speed = coalesce(.data$step_speed, .data$early_track_speed, .data$mean_speed)
    ) |>
    ungroup() |>
    mutate(
      mean_speed = if_else(is.nan(.data$mean_speed), NA_real_, .data$mean_speed),
      early_track_speed = if_else(is.nan(.data$early_track_speed), NA_real_, .data$early_track_speed)
    )
}

scale_with_stats <- function(x, center, scale) {
  (x - center) / scale
}

predict_dropout_probability <- function(model, model_name, data, scaler) {
  new_data <- data |>
    mutate(
      speed_z = scale_with_stats(.data$model_speed, scaler$model_speed[["center"]], scaler$model_speed[["scale"]]),
      distance_to_edge_z = scale_with_stats(
        .data$distance_to_edge,
        scaler$distance_to_edge[["center"]],
        scaler$distance_to_edge[["scale"]]
      ),
      nearest_neighbor_distance_z = scale_with_stats(
        .data$nearest_neighbor_distance,
        scaler$nearest_neighbor_distance[["center"]],
        scaler$nearest_neighbor_distance[["scale"]]
      ),
      confluency_z = scale_with_stats(.data$confluency, scaler$confluency[["center"]], scaler$confluency[["scale"]]),
      frame_z = scale_with_stats(.data$frame, scaler$frame[["center"]], scaler$frame[["scale"]]),
      dose_log10_z = scale_with_stats(.data$dose_log10, scaler$dose_log10[["center"]], scaler$dose_log10[["scale"]])
    )

  plogis(predict(model, newdata = new_data, type = "link", re.form = NA, allow.new.levels = TRUE))
}

censor_tracks_with_model <- function(tracks, model, model_name, scaler, seed = 1L) {
  set.seed(seed)
  prediction_rows <- tracks |>
    filter(.data$frame < .data$site_last_frame, is.finite(.data$model_speed)) |>
    mutate(random_uniform = runif(n()))
  prediction_rows$dropout_probability <- predict_dropout_probability(model, model_name, prediction_rows, scaler)

  prediction_rows <- prediction_rows |>
    arrange(.data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$site_id, .data$migration_track_id) |>
    summarize(
      censor_frame = {
        hit <- .data$frame[.data$random_uniform < .data$dropout_probability]
        if (length(hit) == 0L) Inf else min(hit)
      },
      .groups = "drop"
    )

  tracks |>
    left_join(prediction_rows, by = c("site_id", "migration_track_id")) |>
    filter(is.infinite(.data$censor_frame) | .data$frame <= .data$censor_frame) |>
    select(-"censor_frame") |>
    mutate(source = paste0("OU + ", model_name, " censoring"))
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
      speed = sqrt(.data$dx^2 + .data$dy^2),
      frame_window = time_window_label(.data$frame)
    ) |>
    ungroup() |>
    filter(.data$step_frames == 1L, is.finite(.data$speed))
}

summarize_msd <- function(tracks, max_lag = 10L) {
  tracks |>
    arrange(.data$source, .data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$site_id, .data$migration_track_id) |>
    group_modify(~ {
      df <- .x |> arrange(.data$frame)
      bind_rows(lapply(seq_len(max_lag), function(lag_frame) {
        j <- match(df$frame + lag_frame, df$frame)
        ok <- !is.na(j)
        if (!any(ok)) {
          return(NULL)
        }
        tibble(
          lag_frames = lag_frame,
          squared_displacement_px2 = (df$x[j[ok]] - df$x[ok])^2 + (df$y[j[ok]] - df$y[ok])^2
        )
      }))
    }) |>
    ungroup() |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$lag_frames) |>
    summarize(
      mean_squared_displacement_px2 = mean(.data$squared_displacement_px2, na.rm = TRUE),
      n_displacements = n(),
      .groups = "drop"
    )
}

summarize_velocity_autocorrelation <- function(steps, max_lag = 10L) {
  steps |>
    arrange(.data$source, .data$ploidy, .data$Gemcitabine, .data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$site_id, .data$migration_track_id) |>
    group_modify(~ {
      df <- .x |> arrange(.data$frame)
      bind_rows(lapply(seq_len(max_lag), function(lag_frame) {
        j <- match(df$frame + lag_frame, df$frame)
        ok <- !is.na(j)
        if (!any(ok)) {
          return(NULL)
        }
        numerator <- df$dx[ok] * df$dx[j[ok]] + df$dy[ok] * df$dy[j[ok]]
        denominator <- pmax(df$speed[ok] * df$speed[j[ok]], 1e-9)
        tibble(lag_frames = lag_frame, cosine_autocorrelation = numerator / denominator)
      }))
    }) |>
    ungroup() |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$lag_frames) |>
    summarize(
      mean_cosine_autocorrelation = mean(.data$cosine_autocorrelation, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    )
}

summarize_track_lengths <- function(tracks) {
  tracks |>
    group_by(.data$source, .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$site_id, .data$migration_track_id) |>
    summarize(track_length = n_distinct(.data$frame), .groups = "drop")
}
