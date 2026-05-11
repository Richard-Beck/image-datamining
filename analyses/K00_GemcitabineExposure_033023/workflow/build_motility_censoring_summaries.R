#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(as.character(x))) y else x
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

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
}

analysis_dir_from_script <- function() {
  normalizePath(file.path(dirname(script_path()), ".."), mustWork = TRUE)
}

ensure_columns <- function(df, required, label) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

safe_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

scale_with_stats <- function(x, stats) {
  scale <- stats[["scale"]]
  if (!is.finite(scale) || scale == 0) {
    return(rep(0, length(x)))
  }
  (x - stats[["center"]]) / scale
}

mean_ci <- function(x, w = NULL) {
  ok <- is.finite(x)
  if (!is.null(w)) {
    ok <- ok & is.finite(w) & w > 0
  }
  x <- x[ok]
  if (!is.null(w)) {
    w <- w[ok]
  }
  if (length(x) == 0L) {
    return(tibble(mean = NA_real_, se = NA_real_, n = 0L))
  }
  if (is.null(w)) {
    tibble(mean = mean(x), se = sd(x) / sqrt(length(x)), n = length(x))
  } else {
    mu <- weighted.mean(x, w)
    neff <- sum(w)^2 / sum(w^2)
    variance <- sum(w * (x - mu)^2) / sum(w)
    tibble(mean = mu, se = sqrt(variance / neff), n = length(x))
  }
}

condition_fields <- function(df) {
  df |>
    mutate(
      ploidy = factor(as.character(.data$ploidy), levels = c("2N", "4N")),
      Gemcitabine = as.numeric(.data$Gemcitabine),
      dose_label = paste0(.data$Gemcitabine, " nM"),
      dose_label = factor(.data$dose_label, levels = paste0(sort(unique(.data$Gemcitabine)), " nM"))
    )
}

frame_window_label <- function(frame) {
  case_when(
    frame >= 0 & frame < 10 ~ "frames 0-9",
    frame >= 10 & frame < 20 ~ "frames 10-19",
    frame >= 20 & frame <= 40 ~ "frames 20-40",
    TRUE ~ NA_character_
  )
}

compute_steps <- function(tracks) {
  tracks |>
    arrange(.data$site_id, .data$migration_track_id, .data$frame) |>
    group_by(
      .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$well,
      .data$site_id, .data$trackId, .data$migration_track_id
    ) |>
    mutate(
      next_frame = lead(.data$frame),
      dx = lead(.data$x) - .data$x,
      dy = lead(.data$y) - .data$y,
      step_frames = .data$next_frame - .data$frame,
      step_length_px = sqrt(.data$dx^2 + .data$dy^2),
      squared_step_px2 = .data$step_length_px^2,
      frame_window = frame_window_label(.data$frame)
    ) |>
    ungroup() |>
    filter(.data$step_frames == 1L, is.finite(.data$step_length_px))
}

compute_lag_pairs <- function(tracks, max_lag) {
  keys <- c("ploidy", "Gemcitabine", "dose_label", "well", "site_id", "trackId", "migration_track_id")
  base <- tracks |>
    select(all_of(keys), .data$frame, .data$x, .data$y)
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

compute_step_autocorrelation <- function(steps, max_lag) {
  keys <- c("ploidy", "Gemcitabine", "dose_label", "well", "site_id", "trackId", "migration_track_id")
  base <- steps |>
    select(all_of(keys), .data$frame, .data$dx, .data$dy, .data$step_length_px)
  bind_rows(lapply(seq_len(max_lag), function(lag_frame) {
    endpoint <- base |>
      transmute(
        across(all_of(keys)),
        frame = .data$frame - lag_frame,
        dx_end = .data$dx,
        dy_end = .data$dy,
        step_length_end = .data$step_length_px
      )
    base |>
      left_join(endpoint, by = c(keys, "frame")) |>
      filter(is.finite(.data$dx_end), is.finite(.data$dy_end), is.finite(.data$step_length_end)) |>
      transmute(
        across(all_of(keys)),
        frame = .data$frame,
        lag_frames = lag_frame,
        cosine_autocorrelation = (.data$dx * .data$dx_end + .data$dy * .data$dy_end) /
          pmax(.data$step_length_px * .data$step_length_end, 1e-9)
      )
  }))
}

summarize_by_condition <- function(df, value_col, metric, weight_col = NULL) {
  df |>
    group_by(.data$Gemcitabine, .data$dose_label, .data$ploidy) |>
    group_modify(~ {
      if (is.null(weight_col)) {
        mean_ci(.x[[value_col]])
      } else {
        mean_ci(.x[[value_col]], .x[[weight_col]])
      }
    }) |>
    ungroup() |>
    mutate(metric = metric, .before = 1)
}

contrast_4n_2n <- function(summary_df) {
  wide <- summary_df |>
    select(any_of(c("metric", "lag_frames", "ipcw_cap", "Gemcitabine", "dose_label", "ploidy", "mean", "se", "n"))) |>
    pivot_wider(
      names_from = .data$ploidy,
      values_from = c(.data$mean, .data$se, .data$n),
      names_sep = "_"
    )
  for (col in c("mean_2N", "mean_4N", "se_2N", "se_4N", "n_2N", "n_4N")) {
    if (!col %in% names(wide)) {
      wide[[col]] <- NA_real_
    }
  }
  wide |>
    mutate(
      difference_4N_minus_2N = .data$mean_4N - .data$mean_2N,
      ratio_4N_over_2N = .data$mean_4N / .data$mean_2N,
      se_difference = sqrt(.data$se_4N^2 + .data$se_2N^2)
    )
}

prepare_censoring_model_data <- function(at_risk, scaler = NULL) {
  out <- at_risk |>
    mutate(
      ploidy = factor(as.character(.data$ploidy), levels = c("2N", "4N")),
      Gemcitabine = as.numeric(.data$Gemcitabine),
      dose_log10 = log10(.data$Gemcitabine + 1),
      frame_window = frame_window_label(.data$frame),
      well = factor(.data$well)
    )

  if (is.null(scaler)) {
    out |>
      mutate(
        speed_z = safe_scale(.data$model_speed),
        distance_to_edge_z = safe_scale(.data$distance_to_edge),
        nearest_neighbor_distance_z = safe_scale(.data$nearest_neighbor_distance),
        confluency_z = safe_scale(.data$confluency),
        frame_z = safe_scale(.data$frame),
        dose_log10_z = safe_scale(.data$dose_log10)
      )
  } else {
    out |>
      mutate(
        speed_z = scale_with_stats(.data$model_speed, scaler$model_speed),
        distance_to_edge_z = scale_with_stats(.data$distance_to_edge, scaler$distance_to_edge),
        nearest_neighbor_distance_z = scale_with_stats(
          .data$nearest_neighbor_distance,
          scaler$nearest_neighbor_distance
        ),
        confluency_z = scale_with_stats(.data$confluency, scaler$confluency),
        frame_z = scale_with_stats(.data$frame, scaler$frame),
        dose_log10_z = scale_with_stats(.data$dose_log10, scaler$dose_log10)
      )
  }
}

predict_dropout <- function(model, newdata) {
  if (inherits(model, "merMod")) {
    return(predict(model, newdata = newdata, type = "response", re.form = NA, allow.new.levels = TRUE))
  }
  predict(model, newdata = newdata, type = "response")
}

load_dropout_models <- function(at_risk, max_rows, seed, censoring_model_rds) {
  if (is.null(censoring_model_rds) || !nzchar(censoring_model_rds) || !file.exists(censoring_model_rds)) {
    stop(
      "Required censoring model artifact not found: ", censoring_model_rds,
      "\nRun the censoring audit model job first, then rerun this script.",
      call. = FALSE
    )
  }
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("The censoring model artifact requires the lme4 package.", call. = FALSE)
  }
  external_artifact <- readRDS(censoring_model_rds)
  if (is.null(external_artifact$models$context) || is.null(external_artifact$scaler)) {
    stop("Censoring model artifact is missing models$context or scaler: ", censoring_model_rds, call. = FALSE)
  }
  message("Using censoring audit model artifact: ", censoring_model_rds)

  model_data <- prepare_censoring_model_data(
    at_risk,
    scaler = external_artifact$scaler
  ) |>
    filter(
      is.finite(.data$end_next_frame),
      is.finite(.data$speed_z),
      is.finite(.data$distance_to_edge_z),
      is.finite(.data$nearest_neighbor_distance_z),
      is.finite(.data$confluency_z),
      is.finite(.data$frame_z),
      !is.na(.data$ploidy)
    )

  if (max_rows > 0L && nrow(model_data) > max_rows) {
    set.seed(seed)
    model_data <- model_data[sample.int(nrow(model_data), max_rows), , drop = FALSE]
  }

  set.seed(seed + 1L)
  eval_flag <- runif(nrow(model_data)) < 0.2
  if (sum(eval_flag) < 100L || sum(!eval_flag) < 100L) {
    eval_flag <- rep(FALSE, nrow(model_data))
  }
  train_data <- model_data[!eval_flag, , drop = FALSE]
  eval_data <- model_data[eval_flag, , drop = FALSE]
  if (nrow(eval_data) == 0L) {
    train_data <- model_data
    eval_data <- model_data
  }

  context_model <- external_artifact$models$context
  full_model <- if (is.null(external_artifact$models$full)) {
    external_artifact$models$context
  } else {
    external_artifact$models$full
  }
  numerator_model <- glm(end_next_frame ~ frame_z, data = train_data, family = binomial())

  eval_predictions <- eval_data |>
    mutate(
      context_pred = pmin(pmax(predict_dropout(context_model, eval_data), 1e-6), 1 - 1e-6),
      full_pred = pmin(pmax(predict_dropout(full_model, eval_data), 1e-6), 1 - 1e-6)
    )

  quality <- bind_rows(
    eval_predictions |>
      summarize(
        model = "context_only",
        n_eval_rows = n(),
        event_rate = mean(.data$end_next_frame),
        mean_predicted_event_rate = mean(.data$context_pred),
        brier_score = mean((.data$end_next_frame - .data$context_pred)^2),
        log_loss = -mean(.data$end_next_frame * log(.data$context_pred) +
          (1 - .data$end_next_frame) * log(1 - .data$context_pred))
      ),
    eval_predictions |>
      summarize(
        model = "with_ploidy_dose",
        n_eval_rows = n(),
        event_rate = mean(.data$end_next_frame),
        mean_predicted_event_rate = mean(.data$full_pred),
        brier_score = mean((.data$end_next_frame - .data$full_pred)^2),
        log_loss = -mean(.data$end_next_frame * log(.data$full_pred) +
          (1 - .data$end_next_frame) * log(1 - .data$full_pred))
      )
  )

  comparison <- tibble(
    model = c("context_only", "with_ploidy_dose"),
    df = c(attr(logLik(context_model), "df"), attr(logLik(full_model), "df")),
    log_likelihood = c(as.numeric(logLik(context_model)), as.numeric(logLik(full_model))),
    aic = c(AIC(context_model), AIC(full_model)),
    n_rows = nrow(train_data),
    event_rate = mean(train_data$end_next_frame)
  ) |>
    mutate(delta_aic = .data$aic - min(.data$aic))

  list(
    context = context_model,
    full = full_model,
    numerator = numerator_model,
    scaler = external_artifact$scaler,
    source = "censoring_audit_artifact",
    censoring_model_rds = normalizePath(censoring_model_rds, mustWork = TRUE),
    comparison = comparison,
    quality = quality,
    eval_predictions = eval_predictions
  )
}

write_censoring_model_metadata <- function(models, out_dir) {
  metadata <- tibble(
    source = models$source,
    censoring_model_rds = models$censoring_model_rds,
    context_model_class = paste(class(models$context), collapse = ";"),
    full_model_class = paste(class(models$full), collapse = ";"),
    numerator_model_class = paste(class(models$numerator), collapse = ";")
  )
  write_csv(metadata, file.path(out_dir, "layer2_censoring_model_source.csv"))
}

write_context_model_quality <- function(models, out_dir) {
  eval_predictions <- models$eval_predictions |>
    mutate(
      observed_dropout = .data$end_next_frame,
      context_survival = 1 - .data$context_pred,
      speed_bin = paste0("quartile ", ntile(.data$model_speed, 4)),
      edge_bin = paste0("quartile ", ntile(.data$distance_to_edge, 4)),
      nearest_neighbor_bin = paste0("quartile ", ntile(.data$nearest_neighbor_distance, 4)),
      confluency_bin = paste0("quartile ", ntile(.data$confluency, 4)),
      frame_bin = paste0("quartile ", ntile(.data$frame, 4))
    )

  write_csv(models$quality, file.path(out_dir, "layer2_censoring_model_holdout_quality.csv"))

  calibration <- eval_predictions |>
    mutate(prediction_bin = ntile(.data$context_pred, 20)) |>
    group_by(.data$prediction_bin) |>
    summarize(
      n = n(),
      mean_predicted_dropout = mean(.data$context_pred),
      observed_dropout = mean(.data$observed_dropout),
      mean_predicted_survival = mean(.data$context_survival),
      observed_survival = mean(1 - .data$observed_dropout),
      .groups = "drop"
    )
  write_csv(calibration, file.path(out_dir, "layer2_context_model_calibration.csv"))

  covariate_checks <- eval_predictions |>
    pivot_longer(
      cols = c("speed_bin", "edge_bin", "nearest_neighbor_bin", "confluency_bin", "frame_bin"),
      names_to = "diagnostic",
      values_to = "bin"
    ) |>
    group_by(.data$diagnostic, .data$bin) |>
    summarize(
      n = n(),
      mean_predicted_dropout = mean(.data$context_pred),
      observed_dropout = mean(.data$observed_dropout),
      mean_predicted_survival = mean(.data$context_survival),
      observed_survival = mean(1 - .data$observed_dropout),
      .groups = "drop"
    )
  write_csv(covariate_checks, file.path(out_dir, "layer2_context_model_covariate_checks.csv"))

  calibration_plot <- calibration |>
    ggplot(aes(.data$mean_predicted_dropout, .data$observed_dropout)) +
    geom_abline(slope = 1, intercept = 0, color = "gray55", linewidth = 0.4) +
    geom_point(aes(size = .data$n), color = "#2b6cb0") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      x = "Mean predicted next-frame dropout",
      y = "Observed next-frame dropout",
      size = "Rows"
    ) +
    theme_bw(base_size = 11)
  ggsave(
    file.path(out_dir, "layer2_context_model_calibration.png"),
    calibration_plot,
    width = 5.5,
    height = 5,
    dpi = 180
  )

  distribution_plot <- eval_predictions |>
    mutate(outcome = if_else(.data$observed_dropout == 1L, "ended next frame", "survived next frame")) |>
    ggplot(aes(.data$context_pred, fill = .data$outcome)) +
    geom_histogram(position = "identity", alpha = 0.55, bins = 50) +
    labs(x = "Predicted next-frame dropout", y = "Held-out rows", fill = NULL) +
    theme_bw(base_size = 11)
  ggsave(
    file.path(out_dir, "layer2_context_model_prediction_distribution.png"),
    distribution_plot,
    width = 7,
    height = 4,
    dpi = 180
  )

  covariate_plot <- covariate_checks |>
    ggplot(aes(.data$bin, .data$observed_dropout, group = 1)) +
    geom_point(color = "#1b4332") +
    geom_line(color = "#1b4332") +
    geom_point(aes(y = .data$mean_predicted_dropout), color = "#b45309") +
    geom_line(aes(y = .data$mean_predicted_dropout), color = "#b45309", linetype = "dashed") +
    facet_wrap(vars(.data$diagnostic), scales = "free_x") +
    labs(x = NULL, y = "Next-frame dropout rate") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(
    file.path(out_dir, "layer2_context_model_covariate_checks.png"),
    covariate_plot,
    width = 9,
    height = 5.5,
    dpi = 180
  )
}

add_ipcw_weights <- function(pairs, at_risk, models, cap_probs, max_ipcw_weight) {
  link_rows <- pairs |>
    select(
      .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$well, .data$site_id,
      .data$trackId, .data$migration_track_id, pair_frame = .data$frame, .data$lag_frames,
      .data$displacement_px, .data$squared_displacement_px2
    ) |>
    tidyr::uncount(.data$lag_frames, .id = "link_index") |>
    mutate(frame = .data$pair_frame + .data$link_index - 1L)

  pred_keys <- link_rows |>
    distinct(.data$site_id, .data$trackId, .data$frame)

  pred_data <- at_risk |>
    inner_join(pred_keys, by = c("site_id", "trackId", "frame")) |>
    distinct(.data$site_id, .data$trackId, .data$frame, .keep_all = TRUE) |>
    prepare_censoring_model_data(scaler = models$scaler)
  pred_data <- pred_data |>
    mutate(
      denom_survival = pmax(1 - predict_dropout(models$context, pred_data), 1e-4),
      numer_survival = pmax(1 - predict_dropout(models$numerator, pred_data), 1e-4)
    ) |>
    transmute(
      site_id = .data$site_id,
      trackId = .data$trackId,
      frame = .data$frame,
      denom_survival = .data$denom_survival,
      numer_survival = .data$numer_survival
    )

  link_survival <- link_rows |>
    left_join(pred_data, by = c("site_id", "trackId", "frame")) |>
    filter(is.finite(.data$denom_survival), is.finite(.data$numer_survival)) |>
    group_by(
      .data$ploidy, .data$Gemcitabine, .data$dose_label, .data$well, .data$site_id,
      .data$trackId, .data$migration_track_id, .data$pair_frame, .data$lag_frames
    ) |>
    summarize(
      displacement_px = first(.data$displacement_px),
      squared_displacement_px2 = first(.data$squared_displacement_px2),
      ipcw_raw = prod(.data$numer_survival) / prod(.data$denom_survival),
      n_links = n(),
      .groups = "drop"
    ) |>
    filter(.data$n_links == .data$lag_frames)

  bind_rows(lapply(cap_probs, function(cap_prob) {
    percentile_cap <- quantile(link_survival$ipcw_raw, cap_prob, na.rm = TRUE, names = FALSE)
    cap <- min(percentile_cap, max_ipcw_weight)
    link_survival |>
      mutate(
        ipcw_cap = paste0("p", round(100 * cap_prob), "_max", max_ipcw_weight),
        ipcw_percentile_cap = percentile_cap,
        ipcw_applied_cap = cap,
        ipcw_weight = pmin(.data$ipcw_raw, cap)
      )
  }))
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
analysis_dir <- normalizePath(args$analysis_dir %||% analysis_dir_from_script(), mustWork = TRUE)
tracks_rds <- normalizePath(
  args$tracks_rds %||% file.path(analysis_dir, "data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"),
  mustWork = TRUE
)
at_risk_rds <- args$at_risk_rds %||% file.path(analysis_dir, "data/yellow_censoring_at_risk.rds")
censoring_model_rds <- args$censoring_model_rds %||%
  file.path(analysis_dir, "data/yellow_censoring_models.rds")
out_dir <- normalizePath(
  args$out_dir %||% file.path(analysis_dir, "data/motility_censoring_summaries"),
  mustWork = FALSE
)
max_lag <- as.integer(args$max_lag %||% "5")
acf_lag <- as.integer(args$acf_lag %||% "2")
fixed_window <- as.integer(args$fixed_window %||% "5")
max_censor_rows <- as.integer(args$max_censor_rows %||% "250000")
max_ipcw_weight <- as.numeric(args$max_ipcw_weight %||% "20")
seed <- as.integer(args$seed %||% "20260509")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading tracks: ", tracks_rds)
tracks_raw <- readRDS(tracks_rds)
ensure_columns(
  tracks_raw,
  c(
    "migration_track_id", "site_id", "well", "position", "frame", "ploidy",
    "Gemcitabine", "trackId", "nucleus_x", "nucleus_y"
  ),
  "tracks_rds"
)
tracks <- tracks_raw |>
  condition_fields() |>
  mutate(x = as.numeric(.data$nucleus_x), y = as.numeric(.data$nucleus_y)) |>
  filter(!is.na(.data$ploidy), is.finite(.data$x), is.finite(.data$y), is.finite(.data$frame))

track_qc <- tracks |>
  group_by(.data$Gemcitabine, .data$dose_label, .data$ploidy, .data$well, .data$site_id, .data$migration_track_id) |>
  summarize(
    first_frame = min(.data$frame),
    last_frame = max(.data$frame),
    track_length = n_distinct(.data$frame),
    .groups = "drop"
  )

steps <- compute_steps(tracks)
lag_pairs <- compute_lag_pairs(tracks, max_lag = max_lag)
autocorr <- compute_step_autocorrelation(steps, max_lag = acf_lag)

step_summary <- bind_rows(
  summarize_by_condition(steps, "step_length_px", "one_frame_step_length_px"),
  summarize_by_condition(steps, "squared_step_px2", "one_frame_squared_step_px2")
)
msd_summary <- lag_pairs |>
  group_by(.data$Gemcitabine, .data$dose_label, .data$ploidy, .data$lag_frames) |>
  group_modify(~ mean_ci(.x$squared_displacement_px2)) |>
  ungroup() |>
  mutate(metric = "msd_px2", .before = 1)
msd_auc <- lag_pairs |>
  group_by(.data$Gemcitabine, .data$dose_label, .data$ploidy, .data$well, .data$site_id, .data$migration_track_id) |>
  summarize(msd_auc_lag1_to_max = mean(.data$squared_displacement_px2, na.rm = TRUE), .groups = "drop") |>
  summarize_by_condition("msd_auc_lag1_to_max", paste0("msd_auc_lag1_to_", max_lag, "_px2"))
acf_summary <- autocorr |>
  group_by(.data$Gemcitabine, .data$dose_label, .data$ploidy, .data$lag_frames) |>
  group_modify(~ mean_ci(.x$cosine_autocorrelation)) |>
  ungroup() |>
  mutate(metric = "step_cosine_autocorrelation", .before = 1)
fixed_window_summary <- lag_pairs |>
  filter(.data$lag_frames == fixed_window) |>
  summarize_by_condition("displacement_px", paste0("fixed_window_", fixed_window, "frame_displacement_px"))

observed_summary <- bind_rows(step_summary, msd_auc, fixed_window_summary)
write_csv(observed_summary, file.path(out_dir, "layer1_observed_process_summary.csv"))
write_csv(msd_summary, file.path(out_dir, "layer1_msd_by_lag.csv"))
write_csv(acf_summary, file.path(out_dir, "layer1_cosine_autocorrelation.csv"))
write_csv(bind_rows(
  contrast_4n_2n(observed_summary),
  contrast_4n_2n(msd_summary),
  contrast_4n_2n(acf_summary)
), file.path(out_dir, "layer1_observed_4N_vs_2N_contrasts.csv"))
write_csv(track_qc, file.path(out_dir, "layer1_track_length_qc_by_track.csv"))
write_csv(
  track_qc |>
    group_by(.data$Gemcitabine, .data$dose_label, .data$ploidy) |>
    summarize(
      n_tracks = n(),
      median_track_length = median(.data$track_length),
      q25_track_length = quantile(.data$track_length, 0.25),
      q75_track_length = quantile(.data$track_length, 0.75),
      .groups = "drop"
    ),
  file.path(out_dir, "layer1_track_length_qc_summary.csv")
)

if (file.exists(at_risk_rds)) {
  message("Loading censoring at-risk data: ", at_risk_rds)
  at_risk <- readRDS(at_risk_rds)$at_risk |>
    as_tibble() |>
    condition_fields()
} else {
  message("No at-risk RDS found; deriving minimal at-risk table from staged tracks.")
  track_last <- track_qc |>
    group_by(.data$migration_track_id) |>
    summarize(last_frame = max(.data$last_frame), .groups = "drop")
  at_risk <- tracks |>
    left_join(track_last, by = "migration_track_id") |>
    group_by(.data$site_id) |>
    mutate(n_frames = max(.data$frame) + 1L) |>
    ungroup() |>
    left_join(
      steps |> select(.data$migration_track_id, .data$frame, model_speed = .data$step_length_px),
      by = c("migration_track_id", "frame")
    ) |>
    group_by(.data$migration_track_id) |>
    mutate(model_speed = coalesce(.data$model_speed, mean(.data$model_speed, na.rm = TRUE))) |>
    ungroup() |>
    mutate(
      distance_to_edge = NA_real_,
      nearest_neighbor_distance = NA_real_,
      confluency = NA_real_,
      end_next_frame = as.integer(.data$frame == .data$last_frame & .data$frame < max(.data$frame))
    )
}

write_csv(
  track_qc |>
    count(.data$Gemcitabine, .data$dose_label, .data$ploidy, name = "yellow_start_tracks"),
  file.path(out_dir, "layer2_yellow_start_counts.csv")
)

censoring_summary <- at_risk |>
  mutate(next_frame_yellow_survival = 1 - .data$end_next_frame) |>
  group_by(.data$Gemcitabine, .data$dose_label, .data$ploidy) |>
  summarize(
    at_risk_rows = n(),
    tracks = n_distinct(.data$migration_track_id),
    next_frame_survival = mean(.data$next_frame_yellow_survival, na.rm = TRUE),
    dropout_rate = mean(.data$end_next_frame, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(censoring_summary, file.path(out_dir, "layer2_next_frame_survival_by_condition.csv"))

bin_summary <- at_risk |>
  mutate(
    speed_bin = paste0("quartile ", ntile(.data$model_speed, 4)),
    nearest_neighbor_bin = paste0("quartile ", ntile(.data$nearest_neighbor_distance, 4)),
    edge_bin = paste0("quartile ", ntile(.data$distance_to_edge, 4)),
    frame_window = frame_window_label(.data$frame)
  ) |>
  pivot_longer(
    cols = c("speed_bin", "nearest_neighbor_bin", "edge_bin", "frame_window"),
    names_to = "diagnostic",
    values_to = "bin"
  ) |>
  group_by(.data$diagnostic, .data$bin, .data$Gemcitabine, .data$dose_label, .data$ploidy) |>
  summarize(
    at_risk_rows = n(),
    next_frame_survival = mean(1 - .data$end_next_frame, na.rm = TRUE),
    dropout_rate = mean(.data$end_next_frame, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(bin_summary, file.path(out_dir, "layer2_survival_by_bins.csv"))

models <- load_dropout_models(
  at_risk,
  max_rows = max_censor_rows,
  seed = seed,
  censoring_model_rds = censoring_model_rds
)
write_censoring_model_metadata(models, out_dir)
write_csv(models$comparison, file.path(out_dir, "layer2_censoring_model_comparison.csv"))
write_context_model_quality(models, out_dir)

ipcw_pairs <- add_ipcw_weights(
  lag_pairs, at_risk, models,
  cap_probs = c(0.95, 0.99),
  max_ipcw_weight = max_ipcw_weight
)
ipcw_msd <- ipcw_pairs |>
  group_by(.data$ipcw_cap, .data$Gemcitabine, .data$dose_label, .data$ploidy, .data$lag_frames) |>
  group_modify(~ mean_ci(.x$squared_displacement_px2, .x$ipcw_weight)) |>
  ungroup() |>
  mutate(metric = "ipcw_msd_px2", .before = 1)
ipcw_auc <- ipcw_pairs |>
  group_by(.data$ipcw_cap, .data$Gemcitabine, .data$dose_label, .data$ploidy, .data$well, .data$site_id, .data$migration_track_id) |>
  summarize(
    msd_auc_lag1_to_max = weighted.mean(.data$squared_displacement_px2, .data$ipcw_weight, na.rm = TRUE),
    mean_ipcw_weight = mean(.data$ipcw_weight, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(.data$ipcw_cap, .data$Gemcitabine, .data$dose_label, .data$ploidy) |>
  group_modify(~ mean_ci(.x$msd_auc_lag1_to_max)) |>
  ungroup() |>
  mutate(metric = paste0("ipcw_msd_auc_lag1_to_", max_lag, "_px2"), .before = 1)
write_csv(ipcw_msd, file.path(out_dir, "layer3_ipcw_msd_by_lag.csv"))
write_csv(ipcw_auc, file.path(out_dir, "layer3_ipcw_msd_auc.csv"))
write_csv(
  bind_rows(contrast_4n_2n(ipcw_msd), contrast_4n_2n(ipcw_auc)),
  file.path(out_dir, "layer3_ipcw_4N_vs_2N_contrasts.csv")
)
write_csv(
  ipcw_pairs |>
    group_by(.data$ipcw_cap) |>
    summarize(
      n_pairs = n(),
      median_weight = median(.data$ipcw_weight),
      percentile_cap = max(.data$ipcw_percentile_cap),
      applied_cap = max(.data$ipcw_applied_cap),
      p95_weight = quantile(.data$ipcw_weight, 0.95),
      p99_weight = quantile(.data$ipcw_weight, 0.99),
      max_weight = max(.data$ipcw_weight),
      .groups = "drop"
    ),
  file.path(out_dir, "layer3_ipcw_weight_qc.csv")
)

observed_plot <- observed_summary |>
  filter(.data$metric %in% c("one_frame_step_length_px", paste0("msd_auc_lag1_to_", max_lag, "_px2"))) |>
  ggplot(aes(.data$dose_label, .data$mean, color = .data$ploidy, group = .data$ploidy)) +
  geom_point(position = position_dodge(width = 0.35)) +
  geom_errorbar(aes(ymin = .data$mean - 1.96 * .data$se, ymax = .data$mean + 1.96 * .data$se),
    width = 0, position = position_dodge(width = 0.35)
  ) +
  facet_wrap(vars(.data$metric), scales = "free_y") +
  labs(x = "Gemcitabine", y = "Observed mean", color = "Ploidy") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "layer1_observed_summary.png"), observed_plot, width = 8, height = 4.5, dpi = 180)

censoring_plot <- censoring_summary |>
  ggplot(aes(.data$dose_label, .data$next_frame_survival, color = .data$ploidy, group = .data$ploidy)) +
  geom_point(position = position_dodge(width = 0.35)) +
  geom_line(position = position_dodge(width = 0.35)) +
  labs(x = "Gemcitabine", y = "Next-frame yellow survival", color = "Ploidy") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "layer2_next_frame_survival.png"), censoring_plot, width = 7, height = 4, dpi = 180)

ipcw_plot <- bind_rows(
  msd_summary |> mutate(ipcw_cap = "observed"),
  ipcw_msd
) |>
  ggplot(aes(.data$lag_frames, .data$mean, color = .data$ploidy, group = .data$ploidy)) +
  geom_point(size = 1.5) +
  geom_line() +
  facet_grid(.data$ipcw_cap ~ .data$dose_label, scales = "free_y") +
  labs(x = "Lag frames", y = "Mean squared displacement px^2", color = "Ploidy") +
  theme_bw(base_size = 10)
ggsave(file.path(out_dir, "layer3_observed_vs_ipcw_msd.png"), ipcw_plot, width = 10, height = 6.5, dpi = 180)

message("Wrote motility/censoring outputs to: ", out_dir)
