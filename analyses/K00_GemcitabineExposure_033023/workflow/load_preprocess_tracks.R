suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

REPO_ROOT <- normalizePath(getwd(), mustWork = TRUE)
ANALYSIS_DIR <- file.path(REPO_ROOT, "analyses/K00_GemcitabineExposure_033023")
DATA_DIR <- file.path(ANALYSIS_DIR, "data")
TRACKING_DIR <- "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Tracking_CSVs"
PLATEMAP_XLSX <- file.path(ANALYSIS_DIR, "Gemcitabine_PlateMap_20240111.xlsx")
OUT_TRACKS_RDS <- file.path(DATA_DIR, "tracking_data.rds")
OUT_SUMMARY_RDS <- file.path(DATA_DIR, "tracking_summaries.rds")
OUT_ISOLATION_SUMMARY_RDS <- file.path(DATA_DIR, "isolation_summaries.rds")
KEEP <- c(
  "frame",
  "trackId",
  "Center_of_the_object_0",
  "Center_of_the_object_1",
  "Object_Area_0",
  "Diameter_0"
)
ISOLATION_MULTIPLIERS <- c(1, 2, 3)
CORES <- 16L
QUIET <- FALSE

source(file.path(ANALYSIS_DIR, "R/preprocessing.R"))

run_preprocessing <- function() {
  message("Step 1/4: reading platemap")
  platemap <- parse_platemap_conditions(PLATEMAP_XLSX)
  message("Step 2/4: loading and preprocessing tracking CSVs")
  tracks <- load_preprocessed_tracking_directory(TRACKING_DIR, keep = KEEP, platemap = platemap, cores = CORES, bind = TRUE)
  message("Step 3/4: summarizing frames")
  summaries <- summarize_tracking_frames(tracks)
  message("Step 4/4: summarizing isolation")
  isolation_summaries <- summarize_isolation_frames(tracks, isolation_multipliers = ISOLATION_MULTIPLIERS)

  message("Writing outputs")
  dir.create(dirname(OUT_TRACKS_RDS), recursive = TRUE, showWarnings = FALSE)
  saveRDS(tracks, OUT_TRACKS_RDS)
  saveRDS(summaries, OUT_SUMMARY_RDS)
  saveRDS(isolation_summaries, OUT_ISOLATION_SUMMARY_RDS)

  if (!QUIET) {
    message("Wrote ", OUT_TRACKS_RDS)
    message("Wrote ", OUT_SUMMARY_RDS)
    message("Wrote ", OUT_ISOLATION_SUMMARY_RDS)
    message("Rows: ", nrow(tracks))
  }

  invisible(list(tracks = tracks, summaries = summaries, isolation_summaries = isolation_summaries))
}

run_preprocessing()
