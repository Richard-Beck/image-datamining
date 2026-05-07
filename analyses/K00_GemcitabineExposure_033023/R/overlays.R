# Reusable tracking overlays for registered K00 image frames.

required_overlay_columns <- function(
  frame_col,
  track_col,
  x_col,
  y_col,
  well_col,
  position_col
) {
  unique(c(frame_col, track_col, x_col, y_col, well_col, position_col))
}

normalize_overlay_image <- function(x) {
  if (length(dim(x)) == 2) {
    qs <- stats::quantile(as.numeric(x), c(0.01, 0.995), na.rm = TRUE)
    if (!is.finite(qs[[1]]) || !is.finite(qs[[2]]) || qs[[2]] <= qs[[1]]) {
      qs <- range(x, na.rm = TRUE)
    }
    x <- pmin(pmax((x - qs[[1]]) / (qs[[2]] - qs[[1]]), 0), 1)
    x[is.na(x)] <- 0
    return(array(rep(x, 3), dim = c(nrow(x), ncol(x), 3)))
  }

  x <- x[, , seq_len(min(3L, dim(x)[3])), drop = FALSE]
  if (max(x, na.rm = TRUE) > 1) {
    x <- x / 255
  }
  x[is.na(x)] <- 0
  x
}

extract_overlay_frame <- function(frames, frame, frame_base = 0L) {
  frame_index <- as.integer(frame) - as.integer(frame_base) + 1L
  if (!is.finite(frame_index) || frame_index < 1L) {
    stop("Requested frame is outside the image stack: ", frame, call. = FALSE)
  }

  if (is.list(frames)) {
    if (frame_index > length(frames)) {
      stop("Requested frame is outside the image stack: ", frame, call. = FALSE)
    }
    return(frames[[frame_index]])
  }
  if (length(dim(frames)) == 4L) {
    if (frame_index > dim(frames)[4]) {
      stop("Requested frame is outside the image stack: ", frame, call. = FALSE)
    }
    return(frames[, , , frame_index])
  }
  if (length(dim(frames)) == 3L && dim(frames)[3] > 3L) {
    if (frame_index > dim(frames)[3]) {
      stop("Requested frame is outside the image stack: ", frame, call. = FALSE)
    }
    return(frames[, , frame_index])
  }
  frames
}

read_overlay_frame <- function(image_path, frame, frame_base = 0L) {
  if (!requireNamespace("tiff", quietly = TRUE)) {
    stop("The tiff package is required to read registered frame images.", call. = FALSE)
  }
  extract_overlay_frame(tiff::readTIFF(image_path, as.is = TRUE, all = TRUE), frame, frame_base = frame_base)
}

draw_overlay_disk <- function(img, x, y, radius = 4, color = c(1, 0.05, 0.1)) {
  h <- dim(img)[1]
  w <- dim(img)[2]
  xi <- round(x)
  yi <- round(y)
  if (!is.finite(xi) || !is.finite(yi) ||
      xi + radius < 1 || xi - radius > w || yi + radius < 1 || yi - radius > h) {
    return(img)
  }

  xs <- max(1, xi - radius):min(w, xi + radius)
  ys <- max(1, yi - radius):min(h, yi + radius)
  grid <- expand.grid(y = ys, x = xs)
  grid <- grid[(grid$x - xi)^2 + (grid$y - yi)^2 <= radius^2, , drop = FALSE]
  for (ch in seq_len(3)) {
    img[cbind(grid$y, grid$x, ch)] <- color[[ch]]
  }
  img
}

draw_overlay_line <- function(img, x0, y0, x1, y1, color = c(1, 0.82, 0), width = 1L) {
  if (!all(is.finite(c(x0, y0, x1, y1)))) {
    return(img)
  }

  h <- dim(img)[1]
  w <- dim(img)[2]
  n <- ceiling(max(abs(x1 - x0), abs(y1 - y0), 1))
  line_points <- unique(data.frame(
    x = round(seq(x0, x1, length.out = n + 1L)),
    y = round(seq(y0, y1, length.out = n + 1L))
  ))
  line_points <- line_points[
    line_points$x >= 1 & line_points$x <= w &
      line_points$y >= 1 & line_points$y <= h,
    ,
    drop = FALSE
  ]
  if (nrow(line_points) == 0L) {
    return(img)
  }

  offsets <- expand.grid(dx = -width:width, dy = -width:width)
  offsets <- offsets[offsets$dx^2 + offsets$dy^2 <= width^2, , drop = FALSE]
  pixels <- merge(line_points, offsets, by = NULL)
  pixels$x <- pixels$x + pixels$dx
  pixels$y <- pixels$y + pixels$dy
  pixels <- unique(pixels[
    pixels$x >= 1 & pixels$x <= w & pixels$y >= 1 & pixels$y <= h,
    c("x", "y"),
    drop = FALSE
  ])
  for (ch in seq_len(3)) {
    img[cbind(pixels$y, pixels$x, ch)] <- color[[ch]]
  }
  img
}

draw_overlay_arrow <- function(img, x0, y0, x1, y1, color = c(1, 0.82, 0), width = 1L) {
  img <- draw_overlay_line(img, x0, y0, x1, y1, color = color, width = width)
  if (!all(is.finite(c(x0, y0, x1, y1)))) {
    return(img)
  }

  angle <- atan2(y1 - y0, x1 - x0)
  for (theta in c(angle + pi - pi / 7, angle + pi + pi / 7)) {
    img <- draw_overlay_line(
      img,
      x1,
      y1,
      x1 + 8 * cos(theta),
      y1 + 8 * sin(theta),
      color = color,
      width = width
    )
  }
  img
}

overlay_crop_bounds <- function(
  tracks,
  image_width,
  image_height,
  crop = c("none", "tracks"),
  crop_bounds = NULL,
  crop_buffer_px = 20,
  x_col = "Center_of_the_object_0",
  y_col = "Center_of_the_object_1"
) {
  crop <- match.arg(crop)
  if (!is.null(crop_bounds)) {
    bounds <- unlist(crop_bounds)
    if (!all(c("xmin", "xmax", "ymin", "ymax") %in% names(bounds))) {
      stop("crop_bounds must include xmin, xmax, ymin, and ymax.", call. = FALSE)
    }
  } else if (identical(crop, "tracks")) {
    bounds <- c(
      xmin = floor(min(tracks[[x_col]], na.rm = TRUE) - crop_buffer_px),
      xmax = ceiling(max(tracks[[x_col]], na.rm = TRUE) + crop_buffer_px),
      ymin = floor(min(tracks[[y_col]], na.rm = TRUE) - crop_buffer_px),
      ymax = ceiling(max(tracks[[y_col]], na.rm = TRUE) + crop_buffer_px)
    )
  } else {
    bounds <- c(xmin = 1, xmax = image_width, ymin = 1, ymax = image_height)
  }

  bounds["xmin"] <- max(1, floor(bounds[["xmin"]]))
  bounds["xmax"] <- min(image_width, ceiling(bounds[["xmax"]]))
  bounds["ymin"] <- max(1, floor(bounds[["ymin"]]))
  bounds["ymax"] <- min(image_height, ceiling(bounds[["ymax"]]))
  if (bounds[["xmin"]] > bounds[["xmax"]] || bounds[["ymin"]] > bounds[["ymax"]]) {
    stop("Invalid crop bounds after clipping to image dimensions.", call. = FALSE)
  }
  bounds
}

frame_to_elapsed_label <- function(frame, frame_interval_hours = 2, frame_base = 0L) {
  elapsed_hours <- (as.integer(frame) - as.integer(frame_base)) * frame_interval_hours
  days <- elapsed_hours %/% 24
  hours <- elapsed_hours %% 24
  sprintf("%02dd%02dh00m", days, hours)
}

# Write a PNG overlay of tracking points and track-step arrows on a registered frame.
#
# `tracks` should already be filtered to one well/position site. Points are drawn
# only for `display_frame`; arrows are drawn for every consecutive frame-to-frame
# step present in `tracks`.
tracking_overlay_png <- function(
  tracks,
  images_dir,
  out_png,
  display_frame = NULL,
  image_path = NULL,
  image_array = NULL,
  frame_col = "frame",
  track_col = "migration_track_id",
  fallback_track_col = "trackId",
  x_col = "Center_of_the_object_0",
  y_col = "Center_of_the_object_1",
  well_col = "well",
  position_col = "position",
  frame_base = 0L,
  crop = c("none", "tracks"),
  crop_bounds = NULL,
  crop_buffer_px = 20,
  point_radius = 4,
  point_color = c(1, 0.05, 0.1),
  point_outline_color = c(0, 0, 0),
  arrow_color = c(1, 0.82, 0),
  arrow_width = 1L,
  arrow_mode = c("all", "current_to_next", "none")
) {
  crop <- match.arg(crop)
  arrow_mode <- match.arg(arrow_mode)
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("The png package is required to write overlay PNGs.", call. = FALSE)
  }

  if (!track_col %in% names(tracks) && fallback_track_col %in% names(tracks)) {
    track_col <- fallback_track_col
  }
  required <- required_overlay_columns(frame_col, track_col, x_col, y_col, well_col, position_col)
  missing_cols <- setdiff(required, names(tracks))
  if (length(missing_cols) > 0L) {
    stop("tracks is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  tracks <- as.data.frame(tracks)
  tracks <- tracks[is.finite(tracks[[frame_col]]) & is.finite(tracks[[x_col]]) & is.finite(tracks[[y_col]]), , drop = FALSE]
  if (nrow(tracks) == 0L) {
    stop("tracks has no finite frame/x/y rows to overlay.", call. = FALSE)
  }

  site_values <- unique(paste(tracks[[well_col]], tracks[[position_col]], sep = "_"))
  if (length(site_values) != 1L) {
    stop("tracks must correspond to exactly one well/position site.", call. = FALSE)
  }
  if (is.null(display_frame)) {
    display_frame <- min(tracks[[frame_col]], na.rm = TRUE)
  }
  display_frame <- as.integer(display_frame)

  well <- unique(tracks[[well_col]])[[1]]
  position <- unique(tracks[[position_col]])[[1]]
  if (is.null(image_path)) {
    image_path <- file.path(images_dir, paste0(well, "_", position, ".tiff"))
  }
  if (!file.exists(image_path)) {
    stop("Image path does not exist: ", image_path, call. = FALSE)
  }

  if (is.null(image_array)) {
    img <- normalize_overlay_image(read_overlay_frame(image_path, display_frame, frame_base = frame_base))
  } else {
    img <- normalize_overlay_image(image_array)
  }
  bounds <- overlay_crop_bounds(
    tracks,
    image_width = dim(img)[2],
    image_height = dim(img)[1],
    crop = crop,
    crop_bounds = crop_bounds,
    crop_buffer_px = crop_buffer_px,
    x_col = x_col,
    y_col = y_col
  )
  img <- img[bounds[["ymin"]]:bounds[["ymax"]], bounds[["xmin"]]:bounds[["xmax"]], , drop = FALSE]

  tracks[[".overlay_x"]] <- tracks[[x_col]] - bounds[["xmin"]] + 1
  tracks[[".overlay_y"]] <- tracks[[y_col]] - bounds[["ymin"]] + 1
  tracks <- tracks[order(tracks[[track_col]], tracks[[frame_col]]), , drop = FALSE]

  if (!identical(arrow_mode, "none")) {
    split_tracks <- split(tracks, tracks[[track_col]], drop = TRUE)
    for (one_track in split_tracks) {
      if (nrow(one_track) < 2L) {
        next
      }
      one_track <- one_track[order(one_track[[frame_col]]), , drop = FALSE]
      for (i in seq_len(nrow(one_track) - 1L)) {
        is_consecutive_step <- one_track[[frame_col]][[i + 1L]] == one_track[[frame_col]][[i]] + 1L
        is_visible_step <- identical(arrow_mode, "all") || one_track[[frame_col]][[i]] == display_frame
        if (is_consecutive_step && is_visible_step) {
          img <- draw_overlay_arrow(
            img,
            one_track[[".overlay_x"]][[i]],
            one_track[[".overlay_y"]][[i]],
            one_track[[".overlay_x"]][[i + 1L]],
            one_track[[".overlay_y"]][[i + 1L]],
            color = arrow_color,
            width = arrow_width
          )
        }
      }
    }
  }

  current <- tracks[tracks[[frame_col]] == display_frame, , drop = FALSE]
  if (nrow(current) > 0L) {
    for (i in seq_len(nrow(current))) {
      img <- draw_overlay_disk(img, current[[".overlay_x"]][[i]], current[[".overlay_y"]][[i]],
        radius = point_radius + 2L,
        color = point_outline_color
      )
      img <- draw_overlay_disk(img, current[[".overlay_x"]][[i]], current[[".overlay_y"]][[i]],
        radius = point_radius,
        color = point_color
      )
    }
  }

  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
  png::writePNG(img, out_png)

  invisible(list(
    output_png = out_png,
    image_path = image_path,
    well = well,
    position = position,
    display_frame = display_frame,
    n_tracks_rows = nrow(tracks),
    n_display_points = nrow(current),
    crop_bounds = bounds
  ))
}

# Render one overlay PNG per frame using stable crop bounds across the sequence.
tracking_overlay_png_sequence <- function(
  tracks,
  images_dir,
  out_dir,
  frames = NULL,
  image_path = NULL,
  output_prefix = NULL,
  frame_col = "frame",
  track_col = "migration_track_id",
  fallback_track_col = "trackId",
  x_col = "Center_of_the_object_0",
  y_col = "Center_of_the_object_1",
  well_col = "well",
  position_col = "position",
  frame_base = 0L,
  crop = c("tracks", "none"),
  crop_bounds = NULL,
  crop_buffer_px = 30,
  arrow_mode = c("current_to_next", "all", "none"),
  point_radius = 4,
  point_color = c(1, 0.05, 0.1),
  point_outline_color = c(0, 0, 0),
  arrow_color = c(1, 0.82, 0),
  arrow_width = 1L
) {
  crop <- match.arg(crop)
  arrow_mode <- match.arg(arrow_mode)
  if (!track_col %in% names(tracks) && fallback_track_col %in% names(tracks)) {
    track_col <- fallback_track_col
  }
  required <- required_overlay_columns(frame_col, track_col, x_col, y_col, well_col, position_col)
  missing_cols <- setdiff(required, names(tracks))
  if (length(missing_cols) > 0L) {
    stop("tracks is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  tracks <- as.data.frame(tracks)
  tracks <- tracks[is.finite(tracks[[frame_col]]) & is.finite(tracks[[x_col]]) & is.finite(tracks[[y_col]]), , drop = FALSE]
  if (nrow(tracks) == 0L) {
    stop("tracks has no finite frame/x/y rows to overlay.", call. = FALSE)
  }

  site_values <- unique(paste(tracks[[well_col]], tracks[[position_col]], sep = "_"))
  if (length(site_values) != 1L) {
    stop("tracks must correspond to exactly one well/position site.", call. = FALSE)
  }
  well <- unique(tracks[[well_col]])[[1]]
  position <- unique(tracks[[position_col]])[[1]]
  if (is.null(image_path)) {
    image_path <- file.path(images_dir, paste0(well, "_", position, ".tiff"))
  }
  if (!file.exists(image_path)) {
    stop("Image path does not exist: ", image_path, call. = FALSE)
  }

  if (is.null(frames)) {
    frames <- sort(unique(as.integer(tracks[[frame_col]])))
  } else {
    frames <- sort(unique(as.integer(frames)))
  }
  if (length(frames) == 0L) {
    stop("No frames were requested for rendering.", call. = FALSE)
  }

  if (!requireNamespace("tiff", quietly = TRUE)) {
    stop("The tiff package is required to read registered frame images.", call. = FALSE)
  }
  image_stack <- tiff::readTIFF(image_path, as.is = TRUE, all = TRUE)
  first_img <- normalize_overlay_image(extract_overlay_frame(image_stack, frames[[1]], frame_base = frame_base))
  sequence_crop_bounds <- overlay_crop_bounds(
    tracks,
    image_width = dim(first_img)[2],
    image_height = dim(first_img)[1],
    crop = crop,
    crop_bounds = crop_bounds,
    crop_buffer_px = crop_buffer_px,
    x_col = x_col,
    y_col = y_col
  )

  if (is.null(output_prefix)) {
    output_prefix <- paste0(well, "_", position, "_tracking_overlay")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  rendered <- vector("list", length(frames))
  for (i in seq_along(frames)) {
    frame <- frames[[i]]
    out_png <- file.path(out_dir, sprintf("%s_frame%03d.png", output_prefix, frame))
    rendered[[i]] <- tracking_overlay_png(
      tracks = tracks,
      images_dir = images_dir,
      out_png = out_png,
      display_frame = frame,
      image_path = image_path,
      image_array = extract_overlay_frame(image_stack, frame, frame_base = frame_base),
      frame_col = frame_col,
      track_col = track_col,
      fallback_track_col = fallback_track_col,
      x_col = x_col,
      y_col = y_col,
      well_col = well_col,
      position_col = position_col,
      frame_base = frame_base,
      crop = crop,
      crop_bounds = sequence_crop_bounds,
      crop_buffer_px = crop_buffer_px,
      point_radius = point_radius,
      point_color = point_color,
      point_outline_color = point_outline_color,
      arrow_color = arrow_color,
      arrow_width = arrow_width,
      arrow_mode = arrow_mode
    )
  }

  frame_paths <- vapply(rendered, function(x) x$output_png, character(1))
  invisible(list(
    frame_paths = frame_paths,
    frames = frames,
    image_path = image_path,
    well = well,
    position = position,
    crop_bounds = sequence_crop_bounds
  ))
}

# Convert ordered PNG frames into a GIF.
png_sequence_gif <- function(
  frame_paths,
  out_gif,
  fps = 2,
  cleanup_frames = FALSE
) {
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop("The magick package is required to write GIF animations.", call. = FALSE)
  }
  if (length(frame_paths) == 0L) {
    stop("frame_paths is empty.", call. = FALSE)
  }
  missing_frames <- frame_paths[!file.exists(frame_paths)]
  if (length(missing_frames) > 0L) {
    stop("Missing PNG frame(s): ", paste(missing_frames, collapse = ", "), call. = FALSE)
  }

  dir.create(dirname(out_gif), recursive = TRUE, showWarnings = FALSE)
  imgs <- magick::image_read(frame_paths)
  anim <- magick::image_animate(imgs, fps = fps)
  magick::image_write(anim, out_gif)

  if (isTRUE(cleanup_frames)) {
    unlink(frame_paths)
  }

  invisible(list(
    output_gif = out_gif,
    frame_paths = frame_paths,
    n_frames = length(frame_paths),
    fps = fps,
    cleanup_frames = cleanup_frames
  ))
}
