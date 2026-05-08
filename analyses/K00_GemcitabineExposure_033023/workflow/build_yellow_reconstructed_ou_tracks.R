#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

usage <- paste0(
  "Usage: build_yellow_reconstructed_ou_tracks.R [options]\n\n",
  "Options:\n",
  "  --repo_root=/path/to/repo\n",
  "  --segments_dir=/path/to/yellow_track_segments_len3\n",
  "  --out_rds=/path/to/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds\n",
  "  --out_summary_csv=/path/to/summary.csv\n",
  "  --out_tracks_csv=/path/to/track_summary.csv\n",
  "  --area_ratio=2\n",
  "  --exclude_negative_track_ids=TRUE\n",
  "  --min_segment_frames=3\n",
  "  --quiet=TRUE\n"
)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
analysis_dir_guess <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(analysis_dir_guess, "R/k00_batch_utils.R"))

parse_site_id <- function(site_id) {
  tibble(
    site_id = site_id,
    well = str_extract(site_id, "^[A-H][0-9]{1,2}"),
    position = suppressWarnings(as.integer(str_match(site_id, "^[A-H][0-9]{1,2}_([0-9]+)$")[, 2]))
  )
}

default_k00_platemap <- function() {
  rows <- LETTERS[1:8]
  cols <- 1:12
  doses <- c(0, 3.125, 6.25, 12.5, 25, 50, 100, 200, 400, 800)
  expand.grid(row = rows, col = cols, stringsAsFactors = FALSE) |>
    as_tibble() |>
    mutate(
      well = paste0(.data$row, .data$col),
      ploidy = case_when(
        .data$col %in% c(1L, 12L) ~ NA_character_,
        .data$row %in% LETTERS[1:4] ~ "2N",
        .data$row %in% LETTERS[5:8] ~ "4N",
        TRUE ~ NA_character_
      ),
      Gemcitabine = doses[match(.data$col, 2:11)],
      condition = if_else(
        !is.na(.data$ploidy) & !is.na(.data$Gemcitabine),
        paste0(
          .data$ploidy,
          "_gem_",
          gsub("\\.", "p", format(.data$Gemcitabine, trim = TRUE, scientific = FALSE)),
          "nM"
        ),
        NA_character_
      )
    ) |>
    arrange(.data$row, .data$col)
}

read_yellow_segments <- function(segments_dir) {
  files <- sort(list.files(
    segments_dir,
    pattern = "_yellow_track_segments_len3_.*\\.tsv$",
    full.names = TRUE
  ))
  if (length(files) == 0L) {
    stop("No yellow segment TSV files found in: ", segments_dir, call. = FALSE)
  }

  message("Reading ", length(files), " yellow segment file(s)")
  rbindlist(lapply(files, function(path) {
    dt <- fread(path)
    dt[, source_file := basename(path)]
    dt
  }), use.names = TRUE, fill = TRUE)
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
repo_root <- normalizePath(args$repo_root %||% normalizePath(file.path(analysis_dir_guess, "../.."), mustWork = TRUE), mustWork = TRUE)
analysis_dir <- file.path(repo_root, "analyses/K00_GemcitabineExposure_033023")
segments_dir <- normalizePath(
  args$segments_dir %||% file.path(analysis_dir, "cpsam_nucleus_tracking_links/yellow_track_segments_len3"),
  mustWork = TRUE
)
out_rds <- normalizePath(
  args$out_rds %||% file.path(analysis_dir, "data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"),
  mustWork = FALSE
)
out_summary_csv <- normalizePath(
  args$out_summary_csv %||% file.path(analysis_dir, "data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3_summary.csv"),
  mustWork = FALSE
)
out_tracks_csv <- normalizePath(
  args$out_tracks_csv %||% file.path(analysis_dir, "data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3_track_summary.csv"),
  mustWork = FALSE
)
area_ratio <- suppressWarnings(as.numeric(args$area_ratio %||% "2"))
exclude_negative_track_ids <- as_flag(args$exclude_negative_track_ids, default = TRUE)
min_segment_frames <- suppressWarnings(as.integer(args$min_segment_frames %||% "3"))
quiet <- as_flag(args$quiet, default = FALSE)

if (is.na(area_ratio) || area_ratio < 0) {
  stop("--area_ratio must be a non-negative number", call. = FALSE)
}
if (is.na(min_segment_frames) || min_segment_frames < 2L) {
  stop("--min_segment_frames must be an integer >= 2", call. = FALSE)
}

yellow_segments <- read_yellow_segments(segments_dir)
required <- c(
  "site_id", "trackId", "frame", "cpsam_label", "cpsam_area_px", "track_mask_area_px",
  "nucleus_x", "nucleus_y", "labelimageId", "object_status", "touching_flag", "link_status"
)
missing_cols <- setdiff(required, names(yellow_segments))
if (length(missing_cols) > 0L) {
  stop("Yellow segment files are missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

raw_rows <- nrow(yellow_segments)
raw_segment_ids <- uniqueN(yellow_segments$segment_id)

points <- yellow_segments |>
  as_tibble() |>
  distinct(
    .data$site_id,
    .data$trackId,
    .data$frame,
    .data$cpsam_label,
    .data$nucleus_x,
    .data$nucleus_y,
    .keep_all = TRUE
  ) |>
  mutate(
    trackId_numeric = suppressWarnings(as.numeric(.data$trackId)),
    cpsam_track_area_ratio = .data$cpsam_area_px / .data$track_mask_area_px,
    passes_area_ratio = is.finite(.data$cpsam_track_area_ratio) &
      .data$cpsam_track_area_ratio >= area_ratio,
    passes_track_id_filter = !exclude_negative_track_ids |
      (is.finite(.data$trackId_numeric) & .data$trackId_numeric >= 0)
  )

deduplicated_rows <- nrow(points)

filtered_points <- points |>
  filter(
    .data$passes_area_ratio,
    .data$passes_track_id_filter,
    is.finite(.data$frame),
    is.finite(.data$nucleus_x),
    is.finite(.data$nucleus_y)
  ) |>
  arrange(.data$site_id, .data$trackId, .data$frame, .data$cpsam_label) |>
  group_by(.data$site_id, .data$trackId) |>
  mutate(
    reconstructed_run_index = cumsum(row_number() == 1L | .data$frame != lag(.data$frame) + 1L)
  ) |>
  ungroup() |>
  group_by(.data$site_id, .data$trackId, .data$reconstructed_run_index) |>
  mutate(
    reconstructed_track_frames = n_distinct(.data$frame),
    reconstructed_track_start_frame = min(.data$frame),
    reconstructed_track_end_frame = max(.data$frame)
  ) |>
  ungroup() |>
  filter(.data$reconstructed_track_frames >= min_segment_frames) |>
  mutate(
    migration_track_id = sprintf(
      "%s::%s::run%03d",
      .data$site_id,
      .data$trackId,
      .data$reconstructed_run_index
    )
  )

site_meta <- parse_site_id(unique(filtered_points$site_id))
platemap <- default_k00_platemap()

tracks <- filtered_points |>
  left_join(site_meta, by = "site_id") |>
  left_join(platemap, by = "well") |>
  mutate(
    frame = as.integer(.data$frame),
    trackId = as.character(.data$trackId),
    Center_of_the_object_1 = as.numeric(.data$nucleus_x),
    Center_of_the_object_0 = as.numeric(.data$nucleus_y),
    Object_Area_0 = as.numeric(.data$cpsam_area_px),
    Size_in_pixels_0 = as.numeric(.data$cpsam_area_px),
    yellow_area_ratio_threshold = area_ratio,
    exclude_negative_track_ids = exclude_negative_track_ids,
    min_segment_frames = min_segment_frames,
    source_dataset = "cpsam_yellow_reconstructed"
  ) |>
  select(
    "source_dataset",
    "site_id", "well", "row", "col", "position", "ploidy", "Gemcitabine", "condition",
    "migration_track_id", "trackId", "reconstructed_run_index",
    "reconstructed_track_start_frame", "reconstructed_track_end_frame", "reconstructed_track_frames",
    "frame", "Center_of_the_object_1", "Center_of_the_object_0",
    "Object_Area_0", "Size_in_pixels_0",
    "cpsam_label", "cpsam_area_px", "track_mask_area_px", "cpsam_track_area_ratio",
    "track_mask_area_source", "nucleus_x", "nucleus_y", "labelimageId",
    "object_status", "track_ids_in_cpsam_object", "touching_flag", "link_status",
    "source_file", "yellow_area_ratio_threshold", "exclude_negative_track_ids", "min_segment_frames"
  ) |>
  arrange(.data$site_id, .data$trackId, .data$reconstructed_run_index, .data$frame)

track_summary <- tracks |>
  group_by(
    .data$ploidy, .data$Gemcitabine, .data$condition, .data$site_id, .data$trackId,
    .data$migration_track_id
  ) |>
  summarize(
    n_frames = n_distinct(.data$frame),
    start_frame = min(.data$frame),
    end_frame = max(.data$frame),
    mean_cpsam_track_area_ratio = mean(.data$cpsam_track_area_ratio, na.rm = TRUE),
    .groups = "drop"
  )

summary <- tibble(
  segments_dir = segments_dir,
  out_rds = out_rds,
  area_ratio = area_ratio,
  min_segment_frames = min_segment_frames,
  raw_sliding_window_rows = raw_rows,
  raw_sliding_window_segments = raw_segment_ids,
  deduplicated_point_rows = deduplicated_rows,
  points_after_area_ratio = sum(points$passes_area_ratio, na.rm = TRUE),
  points_after_track_id_filter = sum(points$passes_area_ratio & points$passes_track_id_filter, na.rm = TRUE),
  final_rows = nrow(tracks),
  final_tracks = n_distinct(tracks$migration_track_id),
  final_sites = n_distinct(tracks$site_id),
  final_conditions = n_distinct(paste(tracks$ploidy, tracks$Gemcitabine, sep = "::")),
  median_track_frames = if (nrow(track_summary) > 0L) median(track_summary$n_frames) else NA_real_,
  max_track_frames = if (nrow(track_summary) > 0L) max(track_summary$n_frames) else NA_integer_
)

dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(tracks, out_rds)
write_csv(summary, out_summary_csv)
write_csv(track_summary, out_tracks_csv)

if (!quiet) {
  message("Wrote ", out_rds)
  message("Wrote ", out_summary_csv)
  message("Wrote ", out_tracks_csv)
  print(summary)
}
