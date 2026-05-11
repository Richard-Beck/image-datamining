#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

usage <- paste0(
  "Usage: build_yellow_censoring_site.R [options]\n\n",
  "Options:\n",
  "  --repo_root=/path/to/image-datamining\n",
  "  --analysis_dir=/path/to/analyses/K00_GemcitabineExposure_033023\n",
  "  --site_id=A10_1\n",
  "  --file_index=1\n",
  "  --manifest=/path/to/site_manifest.tsv\n",
  "  --trackpoint_links_dir=/path/to/trackpoint_links\n",
  "  --nearest_distances_dir=/path/to/object_tables\n",
  "  --out_dir=/path/to/data/yellow_censoring_sites/sites\n",
  "  --area_ratio=2\n",
  "  --exclude_negative_track_ids=TRUE\n",
  "  --min_segment_frames=3\n",
  "  --overwrite=FALSE\n"
)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
analysis_dir_guess <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(analysis_dir_guess, "R/k00_batch_utils.R"))

parse_site_id <- function(site_id) {
  data.table(
    site_id = site_id,
    well = sub("^([A-H][0-9]{1,2})_.*$", "\\1", site_id),
    position = as.integer(sub("^[A-H][0-9]{1,2}_([0-9]+)$", "\\1", site_id))
  )
}

default_k00_platemap_dt <- function() {
  rows <- LETTERS[1:8]
  cols <- 1:12
  doses <- c(0, 3.125, 6.25, 12.5, 25, 50, 100, 200, 400, 800)
  dt <- as.data.table(expand.grid(row = rows, col = cols, stringsAsFactors = FALSE))
  dt[, well := paste0(row, col)]
  dt[, ploidy := fifelse(col %in% c(1L, 12L), NA_character_,
    fifelse(row %in% LETTERS[1:4], "2N", fifelse(row %in% LETTERS[5:8], "4N", NA_character_))
  )]
  dt[, Gemcitabine := doses[match(col, 2:11)]]
  dt[, condition := fifelse(
    !is.na(ploidy) & !is.na(Gemcitabine),
    paste0(ploidy, "_gem_", gsub("\\.", "p", format(Gemcitabine, trim = TRUE, scientific = FALSE)), "nM"),
    NA_character_
  )]
  setorder(dt, row, col)
  dt
}

read_one_row_csv <- function(path) {
  if (!file.exists(path)) {
    return(data.table())
  }
  fread(path, nrows = 1L)
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
repo_root <- normalizePath(
  args$repo_root %||% normalizePath(file.path(analysis_dir_guess, "../.."), mustWork = TRUE),
  mustWork = TRUE
)
analysis_dir <- normalizePath(
  args$analysis_dir %||% file.path(repo_root, "analyses/K00_GemcitabineExposure_033023"),
  mustWork = TRUE
)
trackpoint_links_dir <- normalizePath(
  args$trackpoint_links_dir %||% file.path(analysis_dir, "cpsam_nucleus_tracking_links/trackpoint_links"),
  mustWork = TRUE
)
nearest_distances_dir <- normalizePath(
  args$nearest_distances_dir %||% file.path(analysis_dir, "cpsam_full_stacks/nearest_distances/object_tables"),
  mustWork = TRUE
)
out_dir <- normalizePath(
  args$out_dir %||% file.path(analysis_dir, "data/yellow_censoring_sites/sites"),
  mustWork = FALSE
)
manifest_path <- args$manifest %||% file.path(analysis_dir, "data/yellow_censoring_sites/site_manifest.tsv")
file_index <- as.integer(args$file_index %||% Sys.getenv("FILE_INDEX", Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_)))
site_id <- args$site_id %||% NULL
area_ratio <- as.numeric(args$area_ratio %||% "2")
exclude_negative_track_ids <- as_flag(args$exclude_negative_track_ids, default = TRUE)
min_segment_frames <- as.integer(args$min_segment_frames %||% "3")
overwrite <- as_flag(args$overwrite, default = FALSE)

if (is.null(site_id) && file.exists(manifest_path)) {
  if (is.na(file_index) || file_index < 1L) {
    stop("Provide --site_id or --file_index with --manifest.", call. = FALSE)
  }
  manifest <- fread(manifest_path)
  if (file_index > nrow(manifest)) {
    message("file_index is beyond manifest rows; exiting without work.")
    quit(status = 0)
  }
  site_id <- manifest$site_id[[file_index]]
}
if (is.null(site_id) || is.na(site_id) || !nzchar(site_id)) {
  stop("No site_id supplied.", call. = FALSE)
}

out_rds <- file.path(out_dir, paste0(site_id, "_yellow_censoring_at_risk.rds"))
if (file.exists(out_rds) && !overwrite) {
  message("Output exists; rerun with --overwrite=TRUE: ", out_rds)
  quit(status = 0)
}

trackpoint_path <- file.path(
  trackpoint_links_dir,
  paste0(site_id, "_nucleus_trackpoint_cpsam_links_framesall_excluding_frame0.tsv")
)
distance_path <- file.path(nearest_distances_dir, paste0(site_id, "_nearest_object_distances.tsv"))
manifest_row_path <- file.path(analysis_dir, "cpsam_full_stacks/manifest_rows", paste0(site_id, "_cpsam_manifest.csv"))

if (!file.exists(trackpoint_path)) {
  stop("Trackpoint link file not found: ", trackpoint_path, call. = FALSE)
}
if (!file.exists(distance_path)) {
  stop("Nearest-distance file not found: ", distance_path, call. = FALSE)
}

message("Reading trackpoint links: ", trackpoint_path)
dt <- fread(trackpoint_path)
dt[, site_id := site_id]
required <- c(
  "frame", "trackId", "nucleus_x", "nucleus_y", "cpsam_label", "link_status",
  "object_status", "touching_flag", "cpsam_area_px", "track_mask_area_px"
)
missing_cols <- setdiff(required, names(dt))
if (length(missing_cols) > 0L) {
  stop("Trackpoint links are missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

setDT(dt)
dt <- unique(
  dt,
  by = c("site_id", "trackId", "frame", "cpsam_label", "nucleus_x", "nucleus_y")
)
dt[, row_id := .I]
dt[, `:=`(
  frame = as.integer(frame),
  trackId = as.character(trackId),
  cpsam_label = as.integer(cpsam_label),
  nucleus_x = as.numeric(nucleus_x),
  nucleus_y = as.numeric(nucleus_y),
  track_mask_area_px = as.numeric(track_mask_area_px),
  cpsam_area_px = as.numeric(cpsam_area_px),
  touching_flag = as.logical(touching_flag),
  trackId_numeric = suppressWarnings(as.numeric(trackId)),
  cpsam_track_area_ratio = cpsam_area_px / track_mask_area_px
)]
dt[, `:=`(
  passes_link = link_status == "linked_one_cpsam_mask" & !is.na(cpsam_label),
  passes_object = object_status == "single_track" & fifelse(is.na(touching_flag), FALSE, touching_flag == FALSE),
  passes_area_ratio = is.finite(cpsam_track_area_ratio) & cpsam_track_area_ratio >= area_ratio,
  passes_track_id_filter = !exclude_negative_track_ids | (is.finite(trackId_numeric) & trackId_numeric >= 0),
  passes_numeric = is.finite(frame) & is.finite(nucleus_x) & is.finite(nucleus_y)
)]
dt[, passes_point_yellow := passes_link & passes_object & passes_area_ratio & passes_track_id_filter & passes_numeric]
dt[, in_yellow_run_len3 := FALSE]

eligible <- dt[passes_point_yellow == TRUE]
if (nrow(eligible) > 0L) {
  setorder(eligible, site_id, trackId, frame, cpsam_label)
  eligible[, run_index := cumsum(seq_len(.N) == 1L | frame != shift(frame) + 1L), by = .(site_id, trackId)]
  eligible[, run_frames := uniqueN(frame), by = .(site_id, trackId, run_index)]
  final_row_ids <- eligible[run_frames >= min_segment_frames, row_id]
  dt[row_id %in% final_row_ids, in_yellow_run_len3 := TRUE]
}

message("Reading nearest distances: ", distance_path)
distances_raw <- fread(distance_path)
distances <- distances_raw[, .(
  site_id,
  frame = as.integer(frame),
  cpsam_label = as.integer(cpsam_label),
  nearest_neighbor_distance = as.numeric(nearest_empty_gap_px),
  nearest_mask_center_distance = as.numeric(nearest_mask_center_distance_px),
  nearest_cpsam_label = as.integer(nearest_cpsam_label),
  n_objects_frame = as.integer(n_objects_frame),
  distance_definition
)]
confluency_by_frame <- distances_raw[, .(
  total_cpsam_area_px = sum(cpsam_area_px, na.rm = TRUE),
  n_cpsam_objects = .N
), by = .(site_id, frame)]

setkey(distances, site_id, frame, cpsam_label)
setkey(dt, site_id, frame, cpsam_label)
dt <- distances[dt]
setkey(confluency_by_frame, site_id, frame)
setkey(dt, site_id, frame)
dt <- confluency_by_frame[dt]

site_meta <- parse_site_id(site_id)
manifest_row <- read_one_row_csv(manifest_row_path)
if (nrow(manifest_row) == 0L) {
  manifest_row <- data.table(site_id = site_id, n_frames = NA_integer_, height = NA_integer_, width = NA_integer_)
}
site_meta <- merge(site_meta, manifest_row[, .(site_id, n_frames, height, width)], by = "site_id", all.x = TRUE)
platemap <- default_k00_platemap_dt()
site_meta <- merge(site_meta, platemap[, .(well, row, col, ploidy, Gemcitabine, condition)], by = "well", all.x = TRUE)
dt <- merge(dt, site_meta, by = "site_id", all.x = TRUE, sort = FALSE)

dt[, `:=`(
  x = as.numeric(nucleus_x),
  y = as.numeric(nucleus_y),
  image_area_px = as.numeric(width) * as.numeric(height),
  site_last_frame = as.integer(n_frames) - 1L,
  dose_label = paste0(Gemcitabine, " nM"),
  dose_log10 = log10(Gemcitabine + 1)
)]
dt[, `:=`(
  confluency = total_cpsam_area_px / image_area_px,
  confluency_percent = 100 * total_cpsam_area_px / image_area_px,
  distance_to_edge = pmin(x, width - x, y, height - y, na.rm = FALSE)
)]

setorder(dt, site_id, trackId, frame)
dt[, `:=`(
  previous_frame = shift(frame),
  previous_x = shift(x),
  previous_y = shift(y)
), by = .(site_id, trackId)]
dt[, step_speed := fifelse(
  frame == previous_frame + 1L,
  sqrt((x - previous_x)^2 + (y - previous_y)^2),
  NA_real_
)]
dt[, mean_speed := mean(step_speed, na.rm = TRUE), by = .(site_id, trackId)]
dt[, early_track_speed := mean(step_speed[frame <= min(frame, na.rm = TRUE) + 3L], na.rm = TRUE), by = .(site_id, trackId)]
dt[is.nan(mean_speed), mean_speed := NA_real_]
dt[is.nan(early_track_speed), early_track_speed := NA_real_]
dt[, model_speed := fcoalesce(step_speed, early_track_speed, mean_speed)]
dt[, current_yellow_object := passes_point_yellow]

next_yellow <- dt[passes_point_yellow == TRUE, .(
  site_id,
  trackId,
  frame = frame - 1L,
  next_yellow_object = TRUE
)]
next_yellow <- unique(next_yellow, by = c("site_id", "trackId", "frame"))
setkey(next_yellow, site_id, trackId, frame)
setkey(dt, site_id, trackId, frame)
dt <- next_yellow[dt]
dt[is.na(next_yellow_object), next_yellow_object := FALSE]

next_raw <- dt[, .(
  site_id,
  trackId,
  frame = frame - 1L,
  raw_next_exists = TRUE,
  next_passes_link = passes_link,
  next_passes_object = passes_object,
  next_passes_area_ratio = passes_area_ratio,
  next_passes_track_id_filter = passes_track_id_filter,
  next_passes_numeric = passes_numeric,
  next_passes_point_yellow = passes_point_yellow,
  next_in_yellow_run_len3 = in_yellow_run_len3
)]
next_raw <- unique(next_raw, by = c("site_id", "trackId", "frame"))
setkey(next_raw, site_id, trackId, frame)
setkey(dt, site_id, trackId, frame)
dt <- next_raw[dt]
dt[is.na(raw_next_exists), raw_next_exists := FALSE]

dt[, censor_reason_next := fcase(
  next_yellow_object == TRUE, "observed_next_yellow",
  raw_next_exists == FALSE, "no_raw_next_detection",
  next_passes_link == FALSE, "next_link_failure",
  next_passes_object == FALSE, "next_object_not_single_or_touching",
  next_passes_area_ratio == FALSE, "next_area_ratio_failure",
  next_passes_track_id_filter == FALSE, "next_track_id_filter",
  next_passes_numeric == FALSE, "next_nonfinite_detection",
  next_passes_point_yellow == FALSE, "next_not_point_yellow",
  default = "unknown"
)]
dt[, end_next_frame := as.integer(!next_yellow_object)]

previous_yellow <- dt[passes_point_yellow == TRUE, .(
  site_id,
  trackId,
  frame = frame + 1L,
  has_previous_yellow_parent = TRUE
)]
previous_yellow <- unique(previous_yellow, by = c("site_id", "trackId", "frame"))
setkey(previous_yellow, site_id, trackId, frame)
setkey(dt, site_id, trackId, frame)
dt <- previous_yellow[dt]
dt[is.na(has_previous_yellow_parent), has_previous_yellow_parent := FALSE]
dt[, is_yellow_start := passes_point_yellow == TRUE & !has_previous_yellow_parent]

track_summary <- dt[, .(
  first_frame = min(frame, na.rm = TRUE),
  last_frame = max(frame, na.rm = TRUE),
  track_length = uniqueN(frame),
  mean_speed = mean(step_speed, na.rm = TRUE),
  early_track_speed = mean(step_speed[frame <= min(frame, na.rm = TRUE) + 3L], na.rm = TRUE),
  site_last_frame = max(site_last_frame, na.rm = TRUE),
  start_x = x[which.min(frame)],
  start_y = y[which.min(frame)],
  width = width[which.min(frame)],
  height = height[which.min(frame)],
  well = well[which.min(frame)],
  position = position[which.min(frame)],
  ploidy = ploidy[which.min(frame)],
  Gemcitabine = Gemcitabine[which.min(frame)],
  dose_label = dose_label[which.min(frame)],
  condition = condition[which.min(frame)]
), by = .(migration_track_id = paste(site_id, trackId, sep = "::raw::"), site_id, trackId)]
track_summary[is.nan(mean_speed), mean_speed := NA_real_]
track_summary[is.nan(early_track_speed), early_track_speed := NA_real_]

at_risk <- dt[
  frame < site_last_frame &
    passes_point_yellow == TRUE &
    is.finite(distance_to_edge) &
    is.finite(nearest_neighbor_distance) &
    is.finite(confluency) &
    is.finite(model_speed) &
    !is.na(ploidy),
  .(
    migration_track_id = paste(site_id, trackId, sep = "::raw::"),
    site_id, well, ploidy, Gemcitabine, dose_label, position, condition,
    frame, n_frames, distance_to_edge, nearest_neighbor_distance,
    confluency, confluency_percent, step_speed, dose_log10, x, y,
    last_frame = max(frame, na.rm = TRUE),
    mean_speed, early_track_speed, site_last_frame, model_speed,
    end_next_frame, censor_reason_next,
    current_yellow_object, next_yellow_object,
    in_yellow_run_len3, passes_point_yellow, passes_link, passes_object,
    passes_area_ratio, passes_track_id_filter, passes_numeric,
    cpsam_track_area_ratio, trackId, cpsam_label, object_status, link_status,
    touching_flag, n_objects_frame
  )
]

yellow_start_pool <- dt[
  is_yellow_start == TRUE &
    is.finite(x) &
    is.finite(y) &
    !is.na(ploidy),
  .(
    yellow_start_id = paste(site_id, trackId, frame, sep = "::yellow_start::"),
    migration_track_id = paste(site_id, trackId, sep = "::raw::"),
    site_id, trackId, well, position, ploidy, Gemcitabine, dose_label, condition,
    first_frame = frame,
    site_last_frame, start_x = x, start_y = y, width, height,
    mean_speed, early_track_speed, dose_log10,
    distance_to_edge, nearest_neighbor_distance, confluency, confluency_percent,
    cpsam_label
  )
]

context_pool <- at_risk[, .(
  ploidy, Gemcitabine, dose_label, frame,
  nearest_neighbor_distance, confluency, confluency_percent
)]

stage_summary <- data.table(
  site_id = site_id,
  raw_trackpoint_rows = nrow(dt),
  raw_tracks = uniqueN(dt$trackId),
  rows_passing_link = dt[passes_link == TRUE, .N],
  rows_passing_object = dt[passes_link == TRUE & passes_object == TRUE, .N],
  rows_passing_area_ratio = dt[passes_link == TRUE & passes_object == TRUE & passes_area_ratio == TRUE, .N],
  rows_passing_track_id_filter = dt[
    passes_link == TRUE & passes_object == TRUE & passes_area_ratio == TRUE & passes_track_id_filter == TRUE,
    .N
  ],
  rows_passing_point_yellow = dt[passes_point_yellow == TRUE, .N],
  yellow_starts = nrow(yellow_start_pool),
  rows_in_yellow_run_len3 = dt[in_yellow_run_len3 == TRUE, .N],
  tracks_in_yellow_run_len3 = uniqueN(dt[in_yellow_run_len3 == TRUE, trackId]),
  at_risk_rows = nrow(at_risk),
  at_risk_tracks = uniqueN(at_risk$migration_track_id),
  event_rate = if (nrow(at_risk) > 0L) mean(at_risk$end_next_frame) else NA_real_,
  area_ratio = area_ratio,
  min_segment_frames = min_segment_frames,
  exclude_negative_track_ids = exclude_negative_track_ids
)

reason_summary <- dt[frame < site_last_frame, .(
  n_rows = .N,
  n_tracks = uniqueN(trackId)
), by = .(site_id, censor_reason_next)]

artifact <- list(
  at_risk = at_risk,
  context_pool = context_pool,
  yellow_start_pool = yellow_start_pool,
  track_summary = track_summary,
  stage_summary = stage_summary,
  reason_summary = reason_summary,
  created_at = Sys.time(),
  inputs = list(
    trackpoint_path = trackpoint_path,
    distance_path = distance_path,
    manifest_row_path = manifest_row_path
  )
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(artifact, out_rds)
message("Saved site censoring artifact: ", out_rds)
