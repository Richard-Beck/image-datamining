suppressPackageStartupMessages({
  library(dplyr)
})

REPO_ROOT <- normalizePath(getwd(), mustWork = TRUE)
ANALYSIS_DIR <- file.path(REPO_ROOT, "analyses/K00_GemcitabineExposure_033023")
FILTERED_TRACKS_RDS <- file.path(ANALYSIS_DIR, "data/tracking_data_isolated_2x_min3.rds")
N_TRACKS <- 50L
N_REPLICATES <- 25L
FRAME_INTERVAL <- 2

source(file.path(ANALYSIS_DIR, "R/ou_velocity_model.R"))

tracks <- readRDS(FILTERED_TRACKS_RDS) |>
  filter(!is.na(migration_track_id)) |>
  mutate(split_track_id = migration_track_id) |>
  transmute(
    split_track_id,
    frame,
    x = Center_of_the_object_1,
    y = Center_of_the_object_0
  )

track_ids <- unique(tracks$split_track_id)
if (length(track_ids) == 0) {
  stop("No migration tracks found in ", FILTERED_TRACKS_RDS, call. = FALSE)
}

tracks_small <- tracks |>
  filter(split_track_id %in% head(track_ids, N_TRACKS))

segments <- split_ou_segments(tracks_small, frame_interval = FRAME_INTERVAL)
if (length(segments) == 0) {
  stop("No valid OU segments generated.", call. = FALSE)
}

start <- log(initial_ou_params(segments))
test_params <- list(
  start = start,
  tau_up = start + c(log(1.5), 0, 0),
  velocity_up = start + c(0, log(1.5), 0),
  noise_up = start + c(0, 0, log(1.5))
)

validation <- bind_rows(lapply(names(test_params), function(name) {
  par <- test_params[[name]]
  r_value <- ou_loglik_r(par, segments)
  cpp_value <- ou_loglik(par, segments)
  tibble(
    parameter_set = name,
    r_loglik = r_value,
    cpp_loglik = cpp_value,
    abs_diff = abs(r_value - cpp_value),
    rel_diff = abs(r_value - cpp_value) / pmax(abs(r_value), 1)
  )
}))

r_time <- system.time({
  for (i in seq_len(N_REPLICATES)) {
    ou_loglik_r(start, segments)
  }
})

cpp_time <- system.time({
  for (i in seq_len(N_REPLICATES)) {
    ou_loglik(start, segments)
  }
})

cat("Segments:", length(segments), "\n")
cat("Total track time:", attr(segments, "total_track_time"), "\n")
cat("Replicates:", N_REPLICATES, "\n\n")

cat("Validation:\n")
print(validation)

cat("\nTiming:\n")
timing <- tibble(
  implementation = c("R", "Rcpp"),
  elapsed_seconds = c(unname(r_time[["elapsed"]]), unname(cpp_time[["elapsed"]])),
  seconds_per_eval = elapsed_seconds / N_REPLICATES
) |>
  mutate(speedup = elapsed_seconds[implementation == "R"] / elapsed_seconds)
print(timing)

max_abs_diff <- max(validation$abs_diff, na.rm = TRUE)
if (!is.finite(max_abs_diff) || max_abs_diff > 1e-6) {
  stop("Rcpp likelihood validation failed; max abs diff = ", max_abs_diff, call. = FALSE)
}
