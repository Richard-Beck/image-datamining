#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

usage <- paste0(
  "Usage: make_ou_tracking_manifest.R [options]\n\n",
  "Options:\n",
  "  --repo_root=/path/to/repo\n",
  "  --input_rds=/path/to/tracking_data_isolated_2x_min3.rds\n",
  "  --out_csv=/path/to/ou_tracking_fits.csv\n",
  "  --out_manifest=/path/to/ou_tracking_manifest.csv\n",
  "  --frame_interval=2\n",
  "  --n_starts=25\n",
  "  --quiet=TRUE\n"
)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
analysis_dir_guess <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(analysis_dir_guess, "R/k00_batch_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
repo_root <- normalizePath(args$repo_root %||% normalizePath(file.path(analysis_dir_guess, "../.."), mustWork = TRUE), mustWork = TRUE)
analysis_dir <- file.path(repo_root, "analyses/K00_GemcitabineExposure_033023")
input_rds <- normalizePath(args$input_rds %||% file.path(analysis_dir, "data/tracking_data_isolated_2x_min3.rds"), mustWork = TRUE)
out_csv <- normalizePath(args$out_csv %||% file.path(analysis_dir, "data/ou_tracking_fits.csv"), mustWork = FALSE)
out_manifest <- normalizePath(args$out_manifest %||% file.path(analysis_dir, "data/ou_tracking_manifest.csv"), mustWork = FALSE)
frame_interval <- as.numeric(args$frame_interval %||% "2")
n_starts <- as.integer(args$n_starts %||% "25")
quiet <- as_flag(args$quiet, default = FALSE)

if (is.na(frame_interval) || frame_interval <= 0) {
  stop("--frame_interval must be a positive number", call. = FALSE)
}
if (is.na(n_starts) || n_starts < 1) {
  stop("--n_starts must be a positive integer", call. = FALSE)
}

tracks <- readRDS(input_rds)
required <- c("ploidy", "Gemcitabine")
if (!all(required %in% names(tracks))) {
  stop("Input RDS must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
}

manifest <- tracks |>
  filter(!is.na(.data$ploidy), !is.na(.data$Gemcitabine)) |>
  distinct(.data$ploidy, .data$Gemcitabine) |>
  arrange(.data$ploidy, .data$Gemcitabine) |>
  mutate(
    job_id = row_number(),
    input_rds = input_rds,
    out_csv = out_csv,
    min_frame = NA_integer_,
    max_frame = NA_integer_,
    frame_interval = frame_interval,
    min_segment_frames = 3L,
    max_tracks = 0L,
    max_segments = 0L,
    n_starts = n_starts,
    seed = 17L,
    fit_label = paste0(.data$ploidy, "_gem_", gsub("\\.", "p", format(.data$Gemcitabine, trim = TRUE, scientific = FALSE)), "nM")
  ) |>
  select(
    job_id, fit_label, input_rds, out_csv, ploidy, Gemcitabine,
    min_frame, max_frame, frame_interval, min_segment_frames,
    max_tracks, max_segments, n_starts, seed
  )

dir.create(dirname(out_manifest), recursive = TRUE, showWarnings = FALSE)
write_csv(manifest, out_manifest)

if (!quiet) {
  message("Wrote ", out_manifest)
  message("Rows: ", nrow(manifest))
}
