#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

usage <- paste0(
  "Usage: render_ou_fit_sanity_overlays.R [options]\n\n",
  "Render one deterministic random OU sanity-check track GIF. Each array job\n",
  "recomputes the same selected track set from --seed and uses --file_index to\n",
  "choose the one track it should render. No shared manifest is written.\n\n",
  "Options:\n",
  "  --repo_root=/path/to/repo\n",
  "  --input_rds=/path/to/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds\n",
  "  --out_dir=/path/to/ou_fit_sanity\n",
  "  --rendered_csv=/path/to/rendered_overlays.csv\n",
  "  --images_dir=/path/to/Images_40Frames\n",
  "  --file_index=1\n",
  "  --tracks_per_group=16\n",
  "  --crop_buffer_px=40\n",
  "  --seed=17\n",
  "  --fps=2\n",
  "  --quiet=TRUE\n"
)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
analysis_dir_guess <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(analysis_dir_guess, "R/k00_batch_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
repo_root <- normalizePath(args$repo_root %||% normalizePath(file.path(analysis_dir_guess, "../.."), mustWork = TRUE), mustWork = TRUE)
analysis_dir <- file.path(repo_root, "analyses/K00_GemcitabineExposure_033023")
data_dir <- file.path(analysis_dir, "data")

input_rds <- normalizePath(
  args$input_rds %||% file.path(data_dir, "tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"),
  mustWork = TRUE
)
out_dir <- normalizePath(args$out_dir %||% file.path(analysis_dir, "overlay_checks/ou_fit_sanity"), mustWork = FALSE)
rendered_csv <- normalizePath(args$rendered_csv %||% file.path(out_dir, "rendered_overlays.csv"), mustWork = FALSE)
images_dir <- normalizePath(
  args$images_dir %||% Sys.getenv(
    "K00_IMAGES_DIR",
    "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Images_40Frames"
  ),
  mustWork = TRUE
)

file_index <- as.integer(args$file_index %||% Sys.getenv("FILE_INDEX", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
tracks_per_group <- as.integer(args$tracks_per_group %||% Sys.getenv("TRACKS_PER_GROUP", "16"))
crop_buffer_px <- as.integer(args$crop_buffer_px %||% Sys.getenv("CROP_BUFFER_PX", "40"))
seed <- as.integer(args$seed %||% Sys.getenv("SEED", "17"))
fps <- as.integer(args$fps %||% Sys.getenv("FPS", "2"))
quiet <- as_flag(args$quiet, default = FALSE)

if (is.na(file_index) || file_index < 1L) {
  stop("--file_index must be a positive integer", call. = FALSE)
}
if (is.na(tracks_per_group) || tracks_per_group < 1L) {
  stop("--tracks_per_group must be a positive integer", call. = FALSE)
}
if (is.na(crop_buffer_px) || crop_buffer_px < 0L) {
  stop("--crop_buffer_px must be a non-negative integer", call. = FALSE)
}
if (is.na(seed)) {
  stop("--seed must be an integer", call. = FALSE)
}
if (is.na(fps) || fps < 1L) {
  stop("--fps must be a positive integer", call. = FALSE)
}

source(file.path(analysis_dir, "R/overlays.R"))

dose_label <- function(x) {
  gsub("\\.", "p", format(x, trim = TRUE, scientific = FALSE))
}

safe_track_id <- function(x) {
  gsub("[^A-Za-z0-9_=-]+", "_", x)
}

track_summary <- function(tracks) {
  tracks |>
    arrange(.data$migration_track_id, .data$frame) |>
    group_by(.data$migration_track_id, .data$ploidy, .data$Gemcitabine, .data$well, .data$position) |>
    summarize(
      n_frames = n_distinct(.data$frame),
      first_frame = min(.data$frame, na.rm = TRUE),
      last_frame = max(.data$frame, na.rm = TRUE),
      net_displacement_px = sqrt(
        (last(.data$Center_of_the_object_1) - first(.data$Center_of_the_object_1))^2 +
          (last(.data$Center_of_the_object_0) - first(.data$Center_of_the_object_0))^2
      ),
      path_length_px = sum(sqrt(
        diff(.data$Center_of_the_object_1)^2 + diff(.data$Center_of_the_object_0)^2
      ), na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(is.finite(.data$path_length_px), .data$n_frames >= 1L)
}

selected_tracks_from_seed <- function(summary, n_per_group, seed) {
  set.seed(seed)
  summary |>
    group_by(.data$Gemcitabine, .data$ploidy) |>
    group_modify(function(.x, .y) {
      sample_n <- min(n_per_group, nrow(.x))
      .x[sample.int(nrow(.x), size = sample_n), , drop = FALSE]
    }) |>
    ungroup() |>
    arrange(.data$Gemcitabine, .data$ploidy, .data$well, .data$position, .data$migration_track_id) |>
    group_by(.data$Gemcitabine, .data$ploidy) |>
    mutate(
      group_track_index = row_number(),
      condition_label = paste0("dose_", dose_label(.data$Gemcitabine), "_", .data$ploidy)
    ) |>
    ungroup() |>
    mutate(selected_track_index = row_number()) |>
    select(
      "selected_track_index", "condition_label", "group_track_index",
      "migration_track_id", "ploidy", "Gemcitabine", "well", "position",
      "n_frames", "first_frame", "last_frame",
      "path_length_px", "net_displacement_px"
    )
}

global_crop_size <- function(tracks, selected, buffer_px) {
  keys <- c("migration_track_id", "ploidy", "Gemcitabine", "well", "position")
  selected_rows <- tracks |>
    semi_join(selected |> select(all_of(keys)), by = keys)

  crop_sizes <- selected_rows |>
    group_by(.data$migration_track_id, .data$ploidy, .data$Gemcitabine, .data$well, .data$position) |>
    summarize(
      crop_width_needed = ceiling(max(.data$Center_of_the_object_0, na.rm = TRUE) - min(.data$Center_of_the_object_0, na.rm = TRUE) + 1 + 2 * buffer_px),
      crop_height_needed = ceiling(max(.data$Center_of_the_object_1, na.rm = TRUE) - min(.data$Center_of_the_object_1, na.rm = TRUE) + 1 + 2 * buffer_px),
      .groups = "drop"
    )

  c(
    width = as.integer(max(crop_sizes$crop_width_needed, na.rm = TRUE)),
    height = as.integer(max(crop_sizes$crop_height_needed, na.rm = TRUE))
  )
}

image_dimensions <- function(well, position) {
  image_path <- file.path(images_dir, paste0(well, "_", position, ".tiff"))
  if (!file.exists(image_path)) {
    stop("Image path does not exist: ", image_path, call. = FALSE)
  }
  if (!requireNamespace("tiff", quietly = TRUE)) {
    stop("The tiff package is required to inspect registered frame images.", call. = FALSE)
  }
  first_img <- normalize_overlay_image(tiff::readTIFF(image_path, as.is = TRUE, all = FALSE))
  c(width = dim(first_img)[2], height = dim(first_img)[1])
}

centered_crop_bounds <- function(track_rows, frames, crop_width, crop_height, image_width, image_height) {
  rows <- track_rows |>
    filter(.data$frame %in% frames)
  x_center <- mean(range(rows$Center_of_the_object_0, na.rm = TRUE))
  y_center <- mean(range(rows$Center_of_the_object_1, na.rm = TRUE))
  width <- min(as.integer(crop_width), as.integer(image_width))
  height <- min(as.integer(crop_height), as.integer(image_height))

  xmin <- floor(x_center - (width - 1L) / 2)
  ymin <- floor(y_center - (height - 1L) / 2)
  xmin <- max(1L, min(xmin, image_width - width + 1L))
  ymin <- max(1L, min(ymin, image_height - height + 1L))
  c(
    xmin = xmin,
    xmax = xmin + width - 1L,
    ymin = ymin,
    ymax = ymin + height - 1L
  )
}

tracks <- readRDS(input_rds)
required <- c(
  "migration_track_id", "ploidy", "Gemcitabine", "well", "position", "frame",
  "Center_of_the_object_0", "Center_of_the_object_1"
)
missing_cols <- setdiff(required, names(tracks))
if (length(missing_cols) > 0L) {
  stop("Input tracks are missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

summary <- track_summary(tracks)
selected <- selected_tracks_from_seed(summary, n_per_group = tracks_per_group, seed = seed)
if (nrow(selected) == 0L) {
  stop("No tracks were available for overlay rendering.", call. = FALSE)
}
if (file_index > nrow(selected)) {
  if (!quiet) {
    message("--file_index=", file_index, " is beyond the selected track count of ", nrow(selected), "; exiting without work.")
  }
  quit(save = "no", status = 0)
}

crop_size <- global_crop_size(tracks, selected, buffer_px = crop_buffer_px)
selected_row <- selected[file_index, , drop = FALSE]
track_rows <- tracks |>
  filter(
    .data$migration_track_id == selected_row$migration_track_id[[1]],
    .data$ploidy == selected_row$ploidy[[1]],
    .data$Gemcitabine == selected_row$Gemcitabine[[1]],
    .data$well == selected_row$well[[1]],
    .data$position == selected_row$position[[1]]
  )
frames <- sort(unique(as.integer(track_rows$frame)))
dims <- image_dimensions(selected_row$well[[1]], selected_row$position[[1]])
bounds <- centered_crop_bounds(
  track_rows,
  frames = frames,
  crop_width = crop_size[["width"]],
  crop_height = crop_size[["height"]],
  image_width = dims[["width"]],
  image_height = dims[["height"]]
)

condition_dir <- file.path(out_dir, selected_row$condition_label[[1]])
prefix <- sprintf(
  "%s_track%03d_%s",
  selected_row$condition_label[[1]],
  selected_row$group_track_index[[1]],
  safe_track_id(selected_row$migration_track_id[[1]])
)
gif_path <- file.path(condition_dir, paste0(prefix, ".gif"))

if (!quiet) {
  message("Selected track count: ", nrow(selected))
  message("Rendering file_index: ", file_index)
  message("Condition: ", selected_row$condition_label[[1]])
  message("Track: ", selected_row$migration_track_id[[1]])
  message("Frames: ", paste(frames, collapse = ";"))
  message("Rendered frame count: ", length(frames))
  message("Global crop request: ", crop_size[["width"]], "x", crop_size[["height"]])
  message("Render crop bounds: ", paste(names(bounds), bounds, sep = "=", collapse = ", "))
  message("Output GIF: ", gif_path)
  message("Provenance CSV: ", rendered_csv)
}

rendered <- tracking_overlay_png_sequence(
  tracks = track_rows,
  images_dir = images_dir,
  out_dir = condition_dir,
  output_prefix = prefix,
  frames = frames,
  crop = "tracks",
  crop_bounds = bounds,
  crop_buffer_px = crop_buffer_px,
  arrow_mode = "current_to_next"
)
png_sequence_gif(rendered$frame_paths, out_gif = gif_path, fps = fps, cleanup_frames = TRUE)

out_row <- selected_row |>
  mutate(
    rendered_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    seed = seed,
    file_index = file_index,
    selected_track_count = nrow(selected),
    tracks_per_group = tracks_per_group,
    track_length_mode = "all_track_frames",
    track_length = length(frames),
    rendered_frames = paste(frames, collapse = ";"),
    crop_buffer_px = crop_buffer_px,
    requested_crop_width = crop_size[["width"]],
    requested_crop_height = crop_size[["height"]],
    image_width = dims[["width"]],
    image_height = dims[["height"]],
    xmin = bounds[["xmin"]],
    xmax = bounds[["xmax"]],
    ymin = bounds[["ymin"]],
    ymax = bounds[["ymax"]],
    rendered_crop_width = bounds[["xmax"]] - bounds[["xmin"]] + 1L,
    rendered_crop_height = bounds[["ymax"]] - bounds[["ymin"]] + 1L,
    fps = fps,
    gif_path = gif_path
  )

append_delim_locked(out_row, rendered_csv, delim = ",", quiet = quiet)
if (!quiet) {
  message("Rendered ", gif_path)
  message("Appended ", rendered_csv)
}
