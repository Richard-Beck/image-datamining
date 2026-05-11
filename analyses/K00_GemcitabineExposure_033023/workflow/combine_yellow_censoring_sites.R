#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

usage <- paste0(
  "Usage: combine_yellow_censoring_sites.R [options]\n\n",
  "Options:\n",
  "  --repo_root=/path/to/image-datamining\n",
  "  --analysis_dir=/path/to/analyses/K00_GemcitabineExposure_033023\n",
  "  --site_dir=/path/to/staged/sites\n",
  "  --out_rds=/path/to/yellow_censoring_at_risk.rds\n",
  "  --out_summary_csv=/path/to/yellow_censoring_stage_summary.csv\n",
  "  --out_reason_csv=/path/to/yellow_censoring_reason_summary.csv\n",
  "  --quiet=FALSE\n"
)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
analysis_dir_guess <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(analysis_dir_guess, "R/k00_batch_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
repo_root <- normalizePath(
  args$repo_root %||% normalizePath(file.path(analysis_dir_guess, "../.."), mustWork = TRUE),
  mustWork = TRUE
)
analysis_dir <- normalizePath(
  args$analysis_dir %||% file.path(repo_root, "analyses/K00_GemcitabineExposure_033023"),
  mustWork = TRUE
)
site_dir <- normalizePath(
  args$site_dir %||% file.path(analysis_dir, "data/yellow_censoring_sites/sites"),
  mustWork = TRUE
)
out_rds <- normalizePath(
  args$out_rds %||% file.path(analysis_dir, "data/yellow_censoring_at_risk.rds"),
  mustWork = FALSE
)
out_summary_csv <- normalizePath(
  args$out_summary_csv %||% file.path(analysis_dir, "data/yellow_censoring_stage_summary.csv"),
  mustWork = FALSE
)
out_reason_csv <- normalizePath(
  args$out_reason_csv %||% file.path(analysis_dir, "data/yellow_censoring_reason_summary.csv"),
  mustWork = FALSE
)
quiet <- as_flag(args$quiet, default = FALSE)

site_files <- sort(list.files(site_dir, pattern = "_yellow_censoring_at_risk[.]rds$", full.names = TRUE))
if (length(site_files) == 0L) {
  stop("No site censoring RDS files found in ", site_dir, call. = FALSE)
}

if (!quiet) {
  message("Combining ", length(site_files), " site artifact(s)")
}
artifacts <- lapply(site_files, readRDS)
at_risk <- rbindlist(lapply(artifacts, `[[`, "at_risk"), use.names = TRUE, fill = TRUE)
track_summary <- rbindlist(lapply(artifacts, `[[`, "track_summary"), use.names = TRUE, fill = TRUE)
context_pool <- rbindlist(lapply(artifacts, `[[`, "context_pool"), use.names = TRUE, fill = TRUE)
yellow_start_parts <- Filter(Negate(is.null), lapply(artifacts, `[[`, "yellow_start_pool"))
yellow_start_pool <- if (length(yellow_start_parts) > 0L) {
  rbindlist(yellow_start_parts, use.names = TRUE, fill = TRUE)
} else {
  data.table()
}
stage_summary <- rbindlist(lapply(artifacts, `[[`, "stage_summary"), use.names = TRUE, fill = TRUE)
reason_summary <- rbindlist(lapply(artifacts, `[[`, "reason_summary"), use.names = TRUE, fill = TRUE)
if (!"yellow_starts" %in% names(stage_summary)) {
  stage_summary[, yellow_starts := NA_integer_]
}

stage_totals <- stage_summary[, .(
  site_id = "ALL",
  raw_trackpoint_rows = sum(raw_trackpoint_rows, na.rm = TRUE),
  raw_tracks = sum(raw_tracks, na.rm = TRUE),
  rows_passing_link = sum(rows_passing_link, na.rm = TRUE),
  rows_passing_object = sum(rows_passing_object, na.rm = TRUE),
  rows_passing_area_ratio = sum(rows_passing_area_ratio, na.rm = TRUE),
  rows_passing_track_id_filter = sum(rows_passing_track_id_filter, na.rm = TRUE),
  rows_passing_point_yellow = sum(rows_passing_point_yellow, na.rm = TRUE),
  yellow_starts = sum(yellow_starts, na.rm = TRUE),
  rows_in_yellow_run_len3 = sum(rows_in_yellow_run_len3, na.rm = TRUE),
  tracks_in_yellow_run_len3 = sum(tracks_in_yellow_run_len3, na.rm = TRUE),
  at_risk_rows = sum(at_risk_rows, na.rm = TRUE),
  at_risk_tracks = sum(at_risk_tracks, na.rm = TRUE),
  event_rate = if (sum(at_risk_rows, na.rm = TRUE) > 0L) {
    weighted.mean(event_rate, at_risk_rows, na.rm = TRUE)
  } else {
    NA_real_
  },
  area_ratio = unique(area_ratio)[1],
  min_segment_frames = unique(min_segment_frames)[1],
  exclude_negative_track_ids = unique(exclude_negative_track_ids)[1]
)]
stage_summary_out <- rbindlist(list(stage_summary, stage_totals), use.names = TRUE, fill = TRUE)

reason_summary_out <- reason_summary[, .(
  n_rows = sum(n_rows, na.rm = TRUE),
  n_tracks = sum(n_tracks, na.rm = TRUE)
), by = .(censor_reason_next)]
reason_summary_out[, site_id := "ALL"]
setcolorder(reason_summary_out, c("site_id", "censor_reason_next", "n_rows", "n_tracks"))
reason_summary_out <- rbindlist(list(reason_summary, reason_summary_out), use.names = TRUE, fill = TRUE)

artifact <- list(
  at_risk = as.data.frame(at_risk),
  track_summary = as.data.frame(track_summary),
  context_pool = as.data.frame(context_pool),
  yellow_start_pool = as.data.frame(yellow_start_pool),
  stage_summary = as.data.frame(stage_summary_out),
  reason_summary = as.data.frame(reason_summary_out),
  site_files = site_files,
  created_at = Sys.time()
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(artifact, out_rds)
fwrite(stage_summary_out, out_summary_csv)
fwrite(reason_summary_out, out_reason_csv)

if (!quiet) {
  message("Wrote ", out_rds)
  message("Wrote ", out_summary_csv)
  message("Wrote ", out_reason_csv)
  print(stage_totals)
}
