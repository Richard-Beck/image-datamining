#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(tibble)
})

source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(FALSE), value = TRUE
)[1]), mustWork = TRUE))), "R/censoring_audit_utils.R"))

usage <- paste0(
  "Usage: 01_fit_censoring_models.R [options]\n\n",
  "Options:\n",
  "  --analysis_dir=/path/to/analyses/K00_GemcitabineExposure_033023\n",
  "  --at_risk_rds=/path/to/yellow_censoring_at_risk.rds\n",
  "  --out_rds=/path/to/censoring_models.rds\n",
  "  --max_rows=0        Use 0 for all rows; otherwise sample this many at-risk rows.\n",
  "  --seed=20260509\n"
)

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
audit_dir <- audit_dir_from_script()
analysis_dir <- normalizePath(args$analysis_dir %||% analysis_dir_from_audit_dir(audit_dir), mustWork = TRUE)
out_rds <- normalizePath(
  args$out_rds %||% file.path(audit_dir, "artifacts/censoring_models.rds"),
  mustWork = FALSE
)
max_rows <- as.integer(args$max_rows %||% "0")
seed <- as.integer(args$seed %||% "20260509")
at_risk_rds <- args$at_risk_rds %||% file.path(audit_dir, "staged/yellow_censoring_at_risk.rds")

if (file.exists(at_risk_rds)) {
  at_risk_rds <- normalizePath(at_risk_rds, mustWork = TRUE)
  message("Loading precomputed raw censoring at-risk data: ", at_risk_rds)
  observed <- readRDS(at_risk_rds)
  observed$at_risk <- as_tibble(observed$at_risk)
  observed$track_summary <- as_tibble(observed$track_summary)
  observed$context_pool <- as_tibble(observed$context_pool)
  if (is.null(observed$yellow_start_pool)) {
    observed$yellow_start_pool <- observed$track_summary
  }
  observed$yellow_start_pool <- as_tibble(observed$yellow_start_pool)
} else {
  message("Preparing observed at-risk data from already staged yellow tracks: ", analysis_dir)
  observed <- prepare_observed_at_risk(analysis_dir)
  observed$yellow_start_pool <- observed$track_summary
}
model_data <- observed$at_risk |>
  mutate(
    ploidy = factor(.data$ploidy),
    well = factor(.data$well)
  )

if (max_rows > 0L && nrow(model_data) > max_rows) {
  set.seed(seed)
  model_data <- model_data[sample.int(nrow(model_data), max_rows), , drop = FALSE]
}

scaler <- list(
  model_speed = c(center = mean(model_data$model_speed), scale = sd(model_data$model_speed)),
  distance_to_edge = c(center = mean(model_data$distance_to_edge), scale = sd(model_data$distance_to_edge)),
  nearest_neighbor_distance = c(
    center = mean(model_data$nearest_neighbor_distance),
    scale = sd(model_data$nearest_neighbor_distance)
  ),
  confluency = c(center = mean(model_data$confluency), scale = sd(model_data$confluency)),
  frame = c(center = mean(model_data$frame), scale = sd(model_data$frame)),
  dose_log10 = c(center = mean(model_data$dose_log10), scale = sd(model_data$dose_log10))
)

model_data <- model_data |>
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

control <- glmerControl(optimizer = "bobyqa")

message("Fitting context-only mixed dropout model")
context_model <- glmer(
  end_next_frame ~ speed_z + distance_to_edge_z + nearest_neighbor_distance_z +
    confluency_z + frame_z + (1 | well),
  data = model_data,
  family = binomial(),
  nAGQ = 0,
  control = control
)

message("Fitting full mixed dropout model")
full_model <- glmer(
  end_next_frame ~ speed_z + distance_to_edge_z + nearest_neighbor_distance_z +
    confluency_z + frame_z + ploidy + dose_log10_z +
    speed_z:ploidy + speed_z:dose_log10_z + ploidy:dose_log10_z +
    (1 | well),
  data = model_data,
  family = binomial(),
  nAGQ = 0,
  control = control
)

model_comparison <- tibble(
  model = c("context", "full"),
  df = c(attr(logLik(context_model), "df"), attr(logLik(full_model), "df")),
  log_likelihood = c(as.numeric(logLik(context_model)), as.numeric(logLik(full_model))),
  aic = c(AIC(context_model), AIC(full_model))
) |>
  mutate(
    delta_aic = .data$aic - min(.data$aic),
    log_likelihood_gain = .data$log_likelihood - first(.data$log_likelihood)
  )

artifact <- list(
  models = list(context = context_model, full = full_model),
  scaler = scaler,
  model_comparison = model_comparison,
  model_data_summary = tibble(
    n_rows = nrow(model_data),
    n_tracks = n_distinct(model_data$migration_track_id),
    n_wells = n_distinct(model_data$well),
    event_rate = mean(model_data$end_next_frame)
  ),
  context_pool = observed$context_pool,
  yellow_start_pool = observed$yellow_start_pool,
  track_start_pool = observed$yellow_start_pool,
  final_track_summary = observed$track_summary,
  stage_summary = if (!is.null(observed$stage_summary)) observed$stage_summary else NULL,
  reason_summary = if (!is.null(observed$reason_summary)) observed$reason_summary else NULL,
  at_risk_source = if (file.exists(at_risk_rds)) "raw_trackpoint_yellow_censoring" else "staged_yellow_tracks",
  created_at = Sys.time(),
  analysis_dir = analysis_dir
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(artifact, out_rds)
message("Saved censoring model artifact: ", out_rds)
