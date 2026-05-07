suppressPackageStartupMessages({
  library(dplyr)
})

REPO_ROOT <- normalizePath(getwd(), mustWork = TRUE)
ANALYSIS_DIR <- file.path(REPO_ROOT, "analyses/K00_GemcitabineExposure_033023")
DATA_DIR <- file.path(ANALYSIS_DIR, "data")
TRACKING_RDS <- file.path(DATA_DIR, "tracking_data.rds")
OUT_FILTERED_TRACKS_RDS <- file.path(DATA_DIR, "tracking_data_isolated_2x_min3.rds")
OUT_MIGRATION_SUMMARIES_RDS <- file.path(DATA_DIR, "migration_summaries_isolated_2x_min3.rds")
ISOLATION_MULTIPLIER <- 2
MIN_CONTIGUOUS_FRAMES <- 3L
MAX_LAG_FRAMES <- 10L
CORES <- 16L
QUIET <- FALSE

source(file.path(ANALYSIS_DIR, "R/migration_stats.R"))

run_migration_filter_and_summary <- function() {
  message("Step 1/4: loading tracking data")
  tracking <- readRDS(TRACKING_RDS)

  message("Step 2/4: filtering isolated contiguous tracks")
  filtered_tracking <- filter_isolated_contiguous_tracks(
    tracking,
    isolation_multiplier = ISOLATION_MULTIPLIER,
    min_contiguous_frames = MIN_CONTIGUOUS_FRAMES
  )

  message("Step 3/4: computing retention summary")
  retention_by_frame <- summarize_retention_by_frame(tracking, filtered_tracking)

  message("Step 4/4: computing migration summaries")
  migration_summaries <- compute_migration_summaries_parallel(
    filtered_tracking,
    max_lag_frames = MAX_LAG_FRAMES,
    cores = CORES
  )

  migration_summaries <- c(
    list(
      retention_by_frame = retention_by_frame,
      isolation_multiplier = ISOLATION_MULTIPLIER,
      min_contiguous_frames = MIN_CONTIGUOUS_FRAMES,
      max_lag_frames = MAX_LAG_FRAMES
    ),
    migration_summaries
  )

  message("Writing outputs")
  dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(filtered_tracking, OUT_FILTERED_TRACKS_RDS)
  saveRDS(migration_summaries, OUT_MIGRATION_SUMMARIES_RDS)

  if (!QUIET) {
    message("Wrote ", OUT_FILTERED_TRACKS_RDS)
    message("Wrote ", OUT_MIGRATION_SUMMARIES_RDS)
    message("Input rows: ", nrow(tracking))
    message("Filtered rows: ", nrow(filtered_tracking))
    message("Filtered migration tracks: ", dplyr::n_distinct(filtered_tracking$migration_track_id))
  }

  invisible(list(
    filtered_tracking = filtered_tracking,
    migration_summaries = migration_summaries
  ))
}

run_migration_filter_and_summary()
