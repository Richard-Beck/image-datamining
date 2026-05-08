#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

usage <- paste0(
  "Usage: combine_ou_fit_sanity_overlay_grids.R [options]\n\n",
  "Combine one gemcitabine dose worth of per-track GIFs into a 4x8 GIF.\n",
  "For the selected dose, 16 2N GIFs are arranged as a 4x4 grid, 16 4N GIFs\n",
  "are arranged as a 4x4 grid, and those grids are appended side-by-side.\n\n",
  "Options:\n",
  "  --repo_root=/path/to/repo\n",
  "  --out_dir=/path/to/ou_fit_sanity/slurm_seed17_random_tracks_all_frames\n",
  "  --rendered_csv=/path/to/rendered_overlays.csv\n",
  "  --combined_dir=/path/to/combined\n",
  "  --file_index=1\n",
  "  --gemcitabine=0\n",
  "  --tracks_per_ploidy=16\n",
  "  --grid_cols=4\n",
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

out_dir <- normalizePath(args$out_dir %||% file.path(analysis_dir, "overlay_checks/ou_fit_sanity/slurm_seed17_random_tracks_all_frames"), mustWork = FALSE)
rendered_csv <- normalizePath(args$rendered_csv %||% file.path(out_dir, "rendered_overlays.csv"), mustWork = TRUE)
combined_dir <- normalizePath(args$combined_dir %||% file.path(out_dir, "combined"), mustWork = FALSE)
file_index <- as.integer(args$file_index %||% Sys.getenv("FILE_INDEX", Sys.getenv("SLURM_ARRAY_TASK_ID", "1")))
gemcitabine_arg <- args$gemcitabine %||% Sys.getenv("GEMCITABINE", "")
tracks_per_ploidy <- as.integer(args$tracks_per_ploidy %||% Sys.getenv("TRACKS_PER_PLOIDY", "16"))
grid_cols <- as.integer(args$grid_cols %||% Sys.getenv("GRID_COLS", "4"))
fps <- as.integer(args$fps %||% Sys.getenv("FPS", "2"))
quiet <- as_flag(args$quiet, default = FALSE)

if (is.na(file_index) || file_index < 1L) {
  stop("--file_index must be a positive integer", call. = FALSE)
}
if (is.na(tracks_per_ploidy) || tracks_per_ploidy < 1L) {
  stop("--tracks_per_ploidy must be a positive integer", call. = FALSE)
}
if (is.na(grid_cols) || grid_cols < 1L) {
  stop("--grid_cols must be a positive integer", call. = FALSE)
}
if (tracks_per_ploidy %% grid_cols != 0L) {
  stop("--tracks_per_ploidy must be divisible by --grid_cols", call. = FALSE)
}
if (is.na(fps) || fps < 1L) {
  stop("--fps must be a positive integer", call. = FALSE)
}
if (!requireNamespace("magick", quietly = TRUE)) {
  stop("The magick package is required to combine GIF animations.", call. = FALSE)
}

dose_label <- function(x) {
  gsub("\\.", "p", format(x, trim = TRUE, scientific = FALSE))
}

latest_unique_tracks <- function(rendered) {
  rendered |>
    mutate(
      rendered_at_sort = suppressWarnings(as.POSIXct(.data$rendered_at, format = "%Y-%m-%dT%H:%M:%S%z"))
    ) |>
    arrange(.data$Gemcitabine, .data$ploidy, .data$group_track_index, .data$rendered_at_sort) |>
    group_by(.data$Gemcitabine, .data$ploidy, .data$group_track_index) |>
    slice_tail(n = 1) |>
    ungroup()
}

read_gif_frames <- function(path) {
  frames <- magick::image_coalesce(magick::image_read(path))
  if (length(frames) == 0L) {
    stop("GIF has no frames: ", path, call. = FALSE)
  }
  frames
}

append_images <- function(images, stack = FALSE) {
  magick::image_append(magick::image_join(images), stack = stack)
}

add_header <- function(frames, label, fill, color = "white") {
  info <- magick::image_info(frames[[1]])
  header_height <- max(36L, ceiling(info$height[[1]] * 0.045))
  lapply(frames, function(frame) {
    header <- magick::image_blank(width = info$width[[1]], height = header_height, color = fill) |>
      magick::image_annotate(
        text = label,
        gravity = "center",
        color = color,
        size = max(20L, ceiling(header_height * 0.55)),
        weight = 700
      )
    append_images(list(header, frame), stack = TRUE)
  })
}

image_page_size <- function(img) {
  info <- magick::image_info(img[1])
  c(width = info$width[[1]], height = info$height[[1]])
}

normalize_tile_frames <- function(frames, width, height) {
  magick::image_extent(frames, geometry = paste0(width, "x", height), gravity = "center")
}

make_ploidy_grid_frames <- function(paths, output_n_frames, grid_cols) {
  gif_frames <- lapply(paths, read_gif_frames)
  sizes <- do.call(rbind, lapply(gif_frames, image_page_size))
  tile_width <- max(sizes[, "width"])
  tile_height <- max(sizes[, "height"])
  gif_frames <- lapply(gif_frames, normalize_tile_frames, width = tile_width, height = tile_height)

  grid_rows <- length(paths) / grid_cols
  grid_frames <- vector("list", output_n_frames)
  for (frame_idx in seq_len(output_n_frames)) {
    row_imgs <- vector("list", grid_rows)
    for (row_idx in seq_len(grid_rows)) {
      tile_idx <- ((row_idx - 1L) * grid_cols + 1L):(row_idx * grid_cols)
      tiles <- lapply(gif_frames[tile_idx], function(x) {
        x[((frame_idx - 1L) %% length(x)) + 1L]
      })
      row_imgs[[row_idx]] <- append_images(tiles, stack = FALSE)
    }
    grid_frames[[frame_idx]] <- append_images(row_imgs, stack = TRUE)
  }
  list(frames = grid_frames, input_frame_counts = vapply(gif_frames, length, integer(1)))
}

rendered <- read_csv(rendered_csv, show_col_types = FALSE)
required <- c(
  "Gemcitabine", "ploidy", "group_track_index", "track_length", "fps", "gif_path", "rendered_at"
)
missing_cols <- setdiff(required, names(rendered))
if (length(missing_cols) > 0L) {
  stop("rendered_csv is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

rendered <- latest_unique_tracks(rendered)
complete_dose_counts <- rendered |>
  filter(.data$ploidy %in% c("2N", "4N")) |>
  count(.data$Gemcitabine, .data$ploidy, name = "n_gifs") |>
  tidyr::pivot_wider(names_from = "ploidy", values_from = "n_gifs", values_fill = 0L) |>
  filter(.data$`2N` >= tracks_per_ploidy, .data$`4N` >= tracks_per_ploidy) |>
  arrange(.data$Gemcitabine)
doses <- complete_dose_counts$Gemcitabine
if (nzchar(gemcitabine_arg)) {
  selected_dose <- suppressWarnings(as.numeric(gemcitabine_arg))
  if (!is.finite(selected_dose)) {
    stop("--gemcitabine must be numeric when provided.", call. = FALSE)
  }
  if (!selected_dose %in% doses) {
    stop("--gemcitabine=", selected_dose, " is not present in ", rendered_csv, call. = FALSE)
  }
} else {
  if (file_index > length(doses)) {
    if (!quiet) {
      message("--file_index=", file_index, " is beyond the dose count of ", length(doses), "; exiting without work.")
    }
    quit(save = "no", status = 0)
  }
  selected_dose <- doses[[file_index]]
}

dose_rows <- rendered |>
  filter(.data$Gemcitabine == selected_dose, .data$ploidy %in% c("2N", "4N")) |>
  arrange(.data$ploidy, .data$group_track_index)

if (!all(c("2N", "4N") %in% unique(dose_rows$ploidy))) {
  stop("Selected dose does not have both 2N and 4N rendered GIFs: ", selected_dose, call. = FALSE)
}

output_n_frames <- max(dose_rows$track_length, na.rm = TRUE)
if (!is.finite(output_n_frames) || output_n_frames < 1L) {
  stop("Selected dose has no finite track_length values.", call. = FALSE)
}
output_n_frames <- as.integer(output_n_frames)

ploidy_grids <- list()
ploidy_frame_counts <- list()
for (ploidy_value in c("2N", "4N")) {
  one <- dose_rows |>
    filter(.data$ploidy == ploidy_value) |>
    arrange(.data$group_track_index) |>
    slice_head(n = tracks_per_ploidy)

  if (nrow(one) != tracks_per_ploidy) {
    stop("Expected ", tracks_per_ploidy, " ", ploidy_value, " GIFs for dose ", selected_dose, ", found ", nrow(one), ".", call. = FALSE)
  }

  missing_gifs <- one$gif_path[!file.exists(one$gif_path)]
  if (length(missing_gifs) > 0L) {
    stop("Missing GIF(s): ", paste(missing_gifs, collapse = ", "), call. = FALSE)
  }

  grid_result <- make_ploidy_grid_frames(
    paths = one$gif_path,
    output_n_frames = output_n_frames,
    grid_cols = grid_cols
  )
  ploidy_grids[[ploidy_value]] <- grid_result$frames
  ploidy_frame_counts[[ploidy_value]] <- grid_result$input_frame_counts
}
ploidy_grids[["2N"]] <- add_header(ploidy_grids[["2N"]], "2N cells", fill = "#2458a6")
ploidy_grids[["4N"]] <- add_header(ploidy_grids[["4N"]], "4N cells", fill = "#8f2f23")

combined_frames <- vector("list", output_n_frames)
for (frame_idx in seq_len(output_n_frames)) {
  frame_2n <- ploidy_grids[["2N"]][[frame_idx]]
  frame_4n <- ploidy_grids[["4N"]][[frame_idx]]
  combined_frames[[frame_idx]] <- append_images(
    list(frame_2n, frame_4n),
    stack = FALSE
  )
}

dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
out_gif <- file.path(combined_dir, paste0("dose_", dose_label(selected_dose), "_2N_4N_grid4x8.gif"))
anim <- magick::image_animate(magick::image_join(combined_frames), fps = fps)
magick::image_write(anim, out_gif)

out_row <- tibble(
  combined_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  file_index = file_index,
  Gemcitabine = selected_dose,
  tracks_per_ploidy = tracks_per_ploidy,
  grid_cols_per_ploidy = grid_cols,
  output_grid_rows = tracks_per_ploidy / grid_cols,
  output_grid_cols = grid_cols * 2L,
  min_frames_2N = min(ploidy_frame_counts[["2N"]]),
  max_frames_2N = max(ploidy_frame_counts[["2N"]]),
  min_frames_4N = min(ploidy_frame_counts[["4N"]]),
  max_frames_4N = max(ploidy_frame_counts[["4N"]]),
  n_frames = output_n_frames,
  fps = fps,
  rendered_csv = rendered_csv,
  combined_gif_path = out_gif
)
write_csv(out_row, file.path(combined_dir, paste0("dose_", dose_label(selected_dose), "_combined_overlays.csv")))

if (!quiet) {
  message("Dose: ", selected_dose)
  message("Wrote ", out_gif)
}
