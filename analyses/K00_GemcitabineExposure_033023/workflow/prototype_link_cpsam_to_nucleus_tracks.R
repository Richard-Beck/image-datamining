#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(png)
  library(tiff)
})

RAW_TIFF_DIR <- "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Images_40Frames"
TRACKING_DIR <- "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Tracking_CSVs"

usage <- paste0(
  "Usage: prototype_link_cpsam_to_nucleus_tracks.R [options]\n\n",
  "Prototype CPSAM-mask to nucleus-tracking linkage and overlays for one site.\n\n",
  "Options:\n",
  "  --repo_root=/path/to/image-datamining\n",
  "  --site_id=A2_1\n",
  "  --frames=all\n",
  "  --touch_connectivity=8\n",
  "  --file_index=1\n",
  "  --manifest=/path/to/cpsam_nucleus_tracking_link_manifest.tsv\n",
  "  --mask_tiff=/path/to/site_cpsam_masks.tiff\n",
  "  --tracks_csv=/path/to/site_CSV-Table.tiff.csv\n",
  "  --raw_tiff=/path/to/site.tiff\n",
  "  --out_dir=/path/to/output_dir\n",
  "  --overwrite=FALSE\n",
  "  --text_cex=0.95\n"
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
analysis_dir <- file.path(repo_root, "analyses/K00_GemcitabineExposure_033023")

parse_frames <- function(x) {
  x <- gsub("\\s+", "", x)
  if (identical(tolower(x), "all")) return(NULL)
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    parts <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq(parts[[1]], parts[[2]]))
  }
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

read_tiff_stack <- function(path) {
  x <- tiff::readTIFF(path, all = TRUE, as.is = TRUE)
  if (!is.list(x)) x <- list(x)
  x
}

normalize_raw_frame <- function(x) {
  if (length(dim(x)) == 3L) {
    x <- x[, , seq_len(min(3L, dim(x)[3L])), drop = FALSE]
    if (dim(x)[3L] == 1L) x <- x[, , 1L]
  }

  if (length(dim(x)) == 2L) {
    vals <- as.numeric(x)
    qs <- stats::quantile(vals, c(0.01, 0.995), na.rm = TRUE)
    if (!is.finite(qs[[1]]) || !is.finite(qs[[2]]) || qs[[2]] <= qs[[1]]) {
      qs <- range(vals, na.rm = TRUE)
    }
    y <- pmin(pmax((x - qs[[1]]) / (qs[[2]] - qs[[1]]), 0), 1)
    y[is.na(y)] <- 0
    return(array(rep(y, 3L), dim = c(nrow(y), ncol(y), 3L)))
  }

  y <- x
  if (max(y, na.rm = TRUE) > 1) y <- y / max(y, na.rm = TRUE)
  y[is.na(y)] <- 0
  pmin(pmax(y, 0), 1)
}

mask_boundary <- function(mask) {
  foreground <- mask > 0L
  nr <- nrow(mask)
  nc <- ncol(mask)
  boundary <- matrix(FALSE, nr, nc)
  if (nr > 1L) {
    diff_v <- mask[-1L, , drop = FALSE] != mask[-nr, , drop = FALSE]
    boundary[-1L, ] <- boundary[-1L, ] | diff_v
    boundary[-nr, ] <- boundary[-nr, ] | diff_v
  }
  if (nc > 1L) {
    diff_h <- mask[, -1L, drop = FALSE] != mask[, -nc, drop = FALSE]
    boundary[, -1L] <- boundary[, -1L] | diff_h
    boundary[, -nc] <- boundary[, -nc] | diff_h
  }
  boundary & foreground
}

touching_labels <- function(mask, connectivity = 8L) {
  labels <- sort(unique(as.integer(mask[mask > 0L])))
  shifts <- list(c(0L, 1L), c(1L, 0L))
  if (connectivity == 8L) {
    shifts <- c(shifts, list(c(1L, 1L), c(1L, -1L)))
  }

  touch_labels <- integer()
  touch_edges <- 0L
  nr <- nrow(mask)
  nc <- ncol(mask)

  for (shift in shifts) {
    dy <- shift[[1]]
    dx <- shift[[2]]
    y1 <- max(1L, 1L + dy):min(nr, nr + dy)
    x1 <- max(1L, 1L + dx):min(nc, nc + dx)
    y0 <- y1 - dy
    x0 <- x1 - dx

    a <- mask[y0, x0, drop = FALSE]
    b <- mask[y1, x1, drop = FALSE]
    touching <- a > 0L & b > 0L & a != b
    if (any(touching)) {
      touch_edges <- touch_edges + sum(touching)
      touch_labels <- c(touch_labels, as.integer(a[touching]), as.integer(b[touching]))
    }
  }
  list(
    labels = labels,
    touch_labels = sort(unique(touch_labels)),
    touch_edge_count = touch_edges
  )
}

point_label <- function(x, y, mask, width, height) {
  xi <- round(x)
  yi <- round(y)
  if (!is.finite(xi) || !is.finite(yi) || xi < 1L || xi > width || yi < 1L || yi > height) {
    return(NA_integer_)
  }
  lab <- as.integer(mask[yi, xi])
  if (lab > 0L) lab else NA_integer_
}

mask_areas <- function(mask) {
  labels <- as.integer(mask[mask > 0L])
  if (length(labels) < 1L) {
    return(data.table(cpsam_label = integer(), cpsam_area_px = integer()))
  }
  out <- as.data.table(as.data.frame(table(labels), stringsAsFactors = FALSE))
  setnames(out, c("labels", "Freq"), c("cpsam_label", "cpsam_area_px"))
  out[, `:=`(
    cpsam_label = as.integer(as.character(cpsam_label)),
    cpsam_area_px = as.integer(cpsam_area_px)
  )]
  out
}

build_yellow_track_segments <- function(trackpoint_links, segment_length = 3L) {
  if (nrow(trackpoint_links) < segment_length) {
    return(data.table())
  }
  yellow_points <- copy(trackpoint_links[
    object_status == "single_track" &
      link_status == "linked_one_cpsam_mask" &
      !is.na(cpsam_label)
  ])
  if (nrow(yellow_points) < segment_length) {
    return(data.table())
  }
  setorder(yellow_points, trackId, frame)
  yellow_points[, prev_frame := shift(frame), by = trackId]
  yellow_points[, new_run := is.na(prev_frame) | frame != prev_frame + 1L, by = trackId]
  yellow_points[, run_id := cumsum(new_run), by = trackId]

  segments <- yellow_points[, {
    n <- .N
    if (n < segment_length) {
      NULL
    } else {
      rbindlist(lapply(seq_len(n - segment_length + 1L), function(start_i) {
        rows <- .SD[start_i:(start_i + segment_length - 1L)]
        rows[, `:=`(
          segment_start_frame = min(frame),
          segment_end_frame = max(frame),
          segment_position = seq_len(.N)
        )]
        rows
      }))
    }
  }, by = .(trackId, run_id)]
  if (nrow(segments) < 1L) {
    return(data.table())
  }
  segments[, segment_id := sprintf(
    "%s_%s_%03d_%03d",
    site_id,
    trackId,
    segment_start_frame,
    segment_end_frame
  )]
  segments[, site_id := site_id]
  segments[, .(
    site_id,
    segment_id,
    trackId,
    segment_start_frame,
    segment_end_frame,
    segment_position,
    frame,
    cpsam_label,
    cpsam_area_px,
    track_mask_area_px,
    track_mask_area_source,
    nucleus_x,
    nucleus_y,
    labelimageId,
    object_status,
    track_ids_in_cpsam_object,
    touching_flag,
    link_status
  )]
}

draw_track_labels <- function(track_points, text_cex, label_box_alpha = 0.85) {
  if (nrow(track_points) < 1L) return(invisible(NULL))
  labels <- as.character(track_points$trackId)
  pad_x <- 3.5 * text_cex
  pad_y <- 2.5 * text_cex
  widths <- strwidth(labels, units = "user", cex = text_cex, font = 2)
  heights <- strheight(labels, units = "user", cex = text_cex, font = 2)
  rect(
    xleft = track_points$x - widths / 2 - pad_x,
    ybottom = track_points$y - heights / 2 - pad_y,
    xright = track_points$x + widths / 2 + pad_x,
    ytop = track_points$y + heights / 2 + pad_y,
    col = adjustcolor("black", alpha.f = label_box_alpha),
    border = NA
  )
  text(track_points$x, track_points$y, labels = labels, col = "white", cex = text_cex, font = 2)
  invisible(NULL)
}

draw_overlay_frame <- function(raw, mask, object_summary, track_points, frame, text_cex) {
  rgb_frame <- normalize_raw_frame(raw)
  height <- dim(rgb_frame)[1L]
  width <- dim(rgb_frame)[2L]

  palette <- list(
    touching = c(1.00, 0.16, 0.00),
    multi_track = c(1.00, 0.00, 1.00),
    no_track = c(0.00, 0.85, 1.00),
    single_track = c(1.00, 0.88, 0.00)
  )

  boundary <- mask_boundary(mask)
  for (st in names(palette)) {
    labels <- object_summary[object_summary[["status"]] == st, cpsam_label]
    if (length(labels) < 1L) next
    status_boundary <- boundary & mask %in% labels
    color <- palette[[st]]
    for (ch in seq_len(3L)) {
      rgb_frame[, , ch][status_boundary] <- color[[ch]]
    }
  }

  tmp_png <- tempfile(fileext = ".png")
  png(tmp_png, width = width, height = height, units = "px", bg = "black", type = "cairo")
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, width), ylim = c(height, 0), asp = 1)
  rasterImage(as.raster(rgb_frame), 0, height, width, 0, interpolate = FALSE)

  legend_x <- 12
  legend_y <- 18
  legend_labels <- c(
    "cyan: no track point",
    "magenta: >1 track points",
    "red: touching object",
    "yellow: one track point"
  )
  legend_cols <- c("cyan", "magenta", "red", "yellow")
  rect(6, 5, 300, 112, col = adjustcolor("black", alpha.f = 0.68), border = NA)
  text(legend_x, legend_y, labels = paste0("frame ", frame), adj = c(0, 0.5), col = "white", cex = 1.0, font = 2)
  for (i in seq_along(legend_labels)) {
    y <- legend_y + 18 * i
    segments(legend_x, y, legend_x + 24, y, col = legend_cols[[i]], lwd = 3)
    text(legend_x + 32, y, labels = legend_labels[[i]], adj = c(0, 0.5), col = "white", cex = 0.78)
  }

  draw_track_labels(track_points, text_cex = text_cex)
  dev.off()

  out <- png::readPNG(tmp_png)
  unlink(tmp_png)
  if (dim(out)[3L] > 3L) out <- out[, , 1:3, drop = FALSE]
  out
}

site_id <- args$site_id %||% "A2_1"
frames_requested <- parse_frames(args$frames %||% "all")
file_index <- as.integer(args$file_index %||% Sys.getenv("FILE_INDEX", Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_)))
touch_connectivity <- as.integer(args$touch_connectivity %||% "8")
if (!touch_connectivity %in% c(4L, 8L)) {
  stop("--touch_connectivity must be 4 or 8", call. = FALSE)
}
overwrite <- as_flag(args$overwrite, default = FALSE)
text_cex <- as.numeric(args$text_cex %||% "0.95")

manifest_path <- normalizePath(
  args$manifest %||% file.path(analysis_dir, "cpsam_nucleus_tracking_links/site_manifest.tsv"),
  mustWork = FALSE
)

manifest_row <- NULL
if (is.null(args$mask_tiff) && file.exists(manifest_path)) {
  if (is.na(file_index) || file_index < 1L) {
    stop("Provide --file_index or SLURM_ARRAY_TASK_ID when using --manifest.", call. = FALSE)
  }
  manifest <- fread(manifest_path)
  if (file_index > nrow(manifest)) {
    message("file_index is beyond manifest rows; exiting without work.")
    quit(status = 0)
  }
  manifest_row <- manifest[file_index]
  site_id <- manifest_row$site_id[[1]]
}

mask_tiff <- normalizePath(
  args$mask_tiff %||% manifest_row$mask_tiff %||% file.path(analysis_dir, "cpsam_full_stacks/masks", paste0(site_id, "_cpsam_masks.tiff")),
  mustWork = TRUE
)
tracks_csv <- normalizePath(
  args$tracks_csv %||% manifest_row$tracks_csv %||% file.path(TRACKING_DIR, paste0(site_id, "_CSV-Table.tiff.csv")),
  mustWork = TRUE
)
raw_tiff <- normalizePath(
  args$raw_tiff %||% manifest_row$raw_tiff %||% file.path(RAW_TIFF_DIR, paste0(site_id, ".tiff")),
  mustWork = TRUE
)
out_dir <- normalizePath(
  args$out_dir %||% file.path(analysis_dir, "cpsam_nucleus_tracking_links"),
  mustWork = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

frame_tag <- if (is.null(frames_requested)) "all_excluding_frame0" else paste0(min(frames_requested), "-", max(frames_requested))
out_tiff <- file.path(out_dir, "overlays", paste0(site_id, "_cpsam_nucleus_track_overlay_frames", frame_tag, ".tiff"))
out_objects <- file.path(out_dir, "object_links", paste0(site_id, "_cpsam_nucleus_object_links_frames", frame_tag, ".tsv"))
out_points <- file.path(out_dir, "trackpoint_links", paste0(site_id, "_nucleus_trackpoint_cpsam_links_frames", frame_tag, ".tsv"))
out_segments <- file.path(out_dir, "yellow_track_segments_len3", paste0(site_id, "_yellow_track_segments_len3_frames", frame_tag, ".tsv"))
dir.create(dirname(out_tiff), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_objects), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_points), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_segments), recursive = TRUE, showWarnings = FALSE)

if (file.exists(out_tiff) && !overwrite) {
  stop("Output exists; rerun with --overwrite=TRUE: ", out_tiff, call. = FALSE)
}

message("Raw TIFF: ", raw_tiff)
message("Mask TIFF: ", mask_tiff)
message("Tracking CSV: ", tracks_csv)
message("Output dir: ", out_dir)

raw_stack <- read_tiff_stack(raw_tiff)
mask_stack <- read_tiff_stack(mask_tiff)
tracks <- fread(tracks_csv)

if (is.null(frames_requested)) {
  frames_requested <- seq.int(1L, min(length(raw_stack), length(mask_stack)) - 1L)
}

required_cols <- c("frame", "trackId", "Center_of_the_object_0", "Center_of_the_object_1")
missing_cols <- setdiff(required_cols, names(tracks))
if (length(missing_cols) > 0L) {
  stop("Tracking CSV is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}
tracks <- tracks[
  frame %in% frames_requested &
    is.finite(Center_of_the_object_0) &
    is.finite(Center_of_the_object_1),
  .(
    frame = as.integer(frame),
    trackId = as.character(trackId),
    labelimageId = if ("labelimageId" %in% names(tracks)) as.integer(labelimageId) else NA_integer_,
    x = as.numeric(Center_of_the_object_0),
    y = as.numeric(Center_of_the_object_1),
    nucleus_x = as.numeric(Center_of_the_object_0),
    nucleus_y = as.numeric(Center_of_the_object_1),
    track_mask_area_px = if ("Size_in_pixels_0" %in% names(tracks)) {
      as.numeric(Size_in_pixels_0)
    } else if ("Object_Area_0" %in% names(tracks)) {
      as.numeric(Object_Area_0)
    } else {
      NA_real_
    },
    track_mask_area_source = if ("Size_in_pixels_0" %in% names(tracks)) {
      "Size_in_pixels_0"
    } else if ("Object_Area_0" %in% names(tracks)) {
      "Object_Area_0"
    } else {
      NA_character_
    }
  )
]

rendered <- vector("list", length(frames_requested))
object_links <- list()
trackpoint_links <- list()

for (idx in seq_along(frames_requested)) {
  frame_value <- frames_requested[[idx]]
  stack_index <- frame_value + 1L
  if (stack_index < 1L || stack_index > length(raw_stack) || stack_index > length(mask_stack)) {
    stop("Frame ", frame_value, " maps outside the available TIFF stack pages.", call. = FALSE)
  }

  mask <- mask_stack[[stack_index]]
  if (length(dim(mask)) > 2L) mask <- mask[, , 1L]
  mask <- round(mask)
  storage.mode(mask) <- "integer"
  touch_info <- touching_labels(mask, connectivity = touch_connectivity)
  area_dt <- mask_areas(mask)
  height <- nrow(mask)
  width <- ncol(mask)

  tracks_frame <- copy(tracks[frame == frame_value])
  if (nrow(tracks_frame) > 0L) {
    labels_at_points <- vapply(
      seq_len(nrow(tracks_frame)),
      function(i) point_label(
        tracks_frame$x[[i]],
        tracks_frame$y[[i]],
        mask,
        width,
        height
      ),
      integer(1)
    )
    tracks_frame[, cpsam_label := labels_at_points]
    tracks_frame[, n_candidate_labels := fifelse(is.na(cpsam_label), 0L, 1L)]
    tracks_frame[, candidate_labels := fifelse(is.na(cpsam_label), "", as.character(cpsam_label))]
  } else {
    tracks_frame[, `:=`(
      candidate_labels = character(),
      n_candidate_labels = integer(),
      cpsam_label = integer()
    )]
  }

  object_dt <- data.table(frame = frame_value, cpsam_label = touch_info$labels)
  object_dt <- merge(object_dt, area_dt, by = "cpsam_label", all.x = TRUE)
  object_dt[, touching_flag := cpsam_label %in% touch_info$touch_labels]
  if (nrow(tracks_frame) > 0L) {
    point_counts <- tracks_frame[!is.na(cpsam_label), .(
      n_track_points = .N,
      track_ids = paste(sort(unique(trackId)), collapse = ",")
    ), by = .(cpsam_label)]
    object_dt <- merge(object_dt, point_counts, by = "cpsam_label", all.x = TRUE)
  } else {
    object_dt[, `:=`(n_track_points = NA_integer_, track_ids = NA_character_)]
  }
  object_dt[is.na(n_track_points), n_track_points := 0L]
  object_dt[is.na(track_ids), track_ids := ""]
  object_dt[, status := fifelse(
    touching_flag,
    "touching",
    fifelse(n_track_points > 1L, "multi_track", fifelse(n_track_points == 0L, "no_track", "single_track"))
  )]
  object_dt[, touch_connectivity := touch_connectivity]
  object_dt[, touch_edge_count_frame := touch_info$touch_edge_count]
  setcolorder(
    object_dt,
    c(
      "frame", "cpsam_label", "status", "n_track_points", "track_ids",
      "cpsam_area_px", "touching_flag", "touch_connectivity", "touch_edge_count_frame"
    )
  )

  tracks_frame[, link_status := fifelse(
    n_candidate_labels == 0L,
    "no_cpsam_mask_at_track_point",
    "linked_one_cpsam_mask"
  )]
  trackpoint_links[[idx]] <- tracks_frame
  object_links[[idx]] <- object_dt
  tracks_frame <- merge(
    tracks_frame,
    object_dt[, .(
      cpsam_label,
      object_status = status,
      cpsam_area_px,
      track_ids_in_cpsam_object = track_ids,
      touching_flag
    )],
    by = "cpsam_label",
    all.x = TRUE,
    sort = FALSE
  )
  trackpoint_links[[idx]] <- tracks_frame

  rendered[[idx]] <- draw_overlay_frame(
    raw_stack[[stack_index]],
    mask,
    object_dt,
    tracks_frame,
    frame = frame_value,
    text_cex = text_cex
  )
  message("Rendered frame ", frame_value, ": ", nrow(object_dt), " CPSAM objects, ", nrow(tracks_frame), " track points")
}

object_links <- rbindlist(object_links, use.names = TRUE, fill = TRUE)
trackpoint_links <- rbindlist(trackpoint_links, use.names = TRUE, fill = TRUE)

invisible(tiff::writeTIFF(rendered, out_tiff, bits.per.sample = 8L, compression = "LZW"))
fwrite(object_links, out_objects, sep = "\t")
fwrite(trackpoint_links, out_points, sep = "\t")
segments <- build_yellow_track_segments(trackpoint_links, segment_length = 3L)
fwrite(segments, out_segments, sep = "\t")

message("Wrote overlay TIFF: ", out_tiff)
message("Wrote object links: ", out_objects)
message("Wrote track-point links: ", out_points)
message("Wrote yellow length-3 track segments: ", out_segments)
