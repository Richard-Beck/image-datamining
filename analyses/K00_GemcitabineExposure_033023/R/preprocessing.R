suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(parallel)
  library(readr)
  library(stringr)
  library(tibble)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x)) || !nzchar(x)) {
    y
  } else {
    x
  }
}

safe_file_stem <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  gsub("[^A-Za-z0-9._-]+", "_", stem)
}

parse_tracking_filename <- function(path) {
  file <- basename(path)
  well <- str_extract(file, "^[A-H][0-9]{1,2}")
  col <- suppressWarnings(as.integer(str_match(file, "^[A-H][0-9]{1,2}_([0-9]+)_")[, 2]))
  row <- if (is.na(well)) NA_character_ else substr(well, 1, 1)
  plate_col <- if (is.na(well)) NA_integer_ else as.integer(str_extract(well, "[0-9]+"))
  position <- col
  plate_id <- basename(dirname(dirname(dirname(normalizePath(path, winslash = "/", mustWork = FALSE)))))

  tibble(
    file = file,
    path = normalizePath(path, winslash = "/", mustWork = FALSE),
    plate_id = plate_id,
    well = well,
    row = row,
    col = plate_col,
    position = position
  )
}

default_platemap_layout <- function() {
  rows <- LETTERS[1:8]
  cols <- 1:12
  doses <- c(0, 3.125, 6.25, 12.5, 25, 50, 100, 200, 400, 800)
  expand.grid(row = rows, col = cols, stringsAsFactors = FALSE) |>
    tibble::as_tibble() |>
    mutate(
      well = paste0(.data$row, .data$col),
      ploidy = case_when(
        .data$col %in% c(1L, 12L) ~ NA_character_,
        .data$row %in% LETTERS[1:4] ~ "2N",
        .data$row %in% LETTERS[5:8] ~ "4N",
        TRUE ~ NA_character_
      ),
      Gemcitabine = if_else(.data$col %in% 2:11, doses[.data$col - 1L], NA_real_),
      condition = if_else(!is.na(.data$ploidy) & !is.na(.data$Gemcitabine),
        paste0(.data$ploidy, "_gem_", format(.data$Gemcitabine, trim = TRUE, scientific = FALSE), "nM"),
        NA_character_
      )
    ) |>
    arrange(.data$row, .data$col)
}

parse_platemap_conditions <- function(platemap_xlsx) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    return(default_platemap_layout())
  }

  raw <- readxl::read_excel(platemap_xlsx, col_names = FALSE, .name_repair = "minimal")
  raw_mat <- as.matrix(raw)
  header <- as.character(raw_mat[1, ])
  plate_cols <- suppressWarnings(as.integer(header))
  row_labels <- as.character(raw[[1]])

  rows <- which(row_labels %in% LETTERS[1:8])
  cols <- which(!is.na(plate_cols) & plate_cols %in% 1:12)

  expand.grid(row_index = rows, col_index = cols, stringsAsFactors = FALSE) |>
    tibble::as_tibble() |>
    mutate(
      row = row_labels[.data$row_index],
      col = plate_cols[.data$col_index],
      well = paste0(.data$row, .data$col),
      treatment = as.character(raw_mat[cbind(.data$row_index, .data$col_index)]),
      treatment = na_if(trimws(.data$treatment), ""),
      ploidy = str_match(.data$treatment, "^(2N|4N):")[, 2],
      Gemcitabine = suppressWarnings(as.numeric(str_match(.data$treatment, ":\\s*([0-9.]+)\\s*nM")[, 2])),
      condition = if_else(!is.na(.data$ploidy) & !is.na(.data$Gemcitabine),
        paste0(.data$ploidy, "_gem_", format(.data$Gemcitabine, trim = TRUE, scientific = FALSE), "nM"),
        NA_character_
      )
    ) |>
    select(well, row, col, treatment, ploidy, Gemcitabine, condition) |>
    arrange(.data$row, .data$col)
}

load_tracking_file <- function(path, keep = NULL, ...) {
  if (!file.exists(path)) {
    stop("Tracking file does not exist: ", path, call. = FALSE)
  }
  if (is.null(keep)) {
    dt <- data.table::fread(path, ...)
  } else {
    dt <- data.table::fread(path, select = keep, ...)
  }
  as_tibble(dt)
}

load_tracking_directory <- function(dir_path, keep = NULL, cores = max(1L, parallel::detectCores(logical = TRUE) - 1L), bind = TRUE, ...) {
  tracking_files <- sort(list.files(dir_path, pattern = "\\.csv$", full.names = TRUE))
  if (length(tracking_files) == 0) {
    stop("No tracking CSV files found in: ", dir_path, call. = FALSE)
  }
  if (cores < 1L) {
    stop("cores must be >= 1", call. = FALSE)
  }
  message("Found ", length(tracking_files), " tracking CSV files")
  message("Loading on ", cores, " worker(s)")
  if (cores == 1L) {
    tracks <- vector("list", length(tracking_files))
    for (i in seq_along(tracking_files)) {
      message(sprintf("[%d/%d] loading %s", i, length(tracking_files), basename(tracking_files[[i]])))
      tracks[[i]] <- load_tracking_file(tracking_files[[i]], keep = keep, ...)
    }
  } else {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(readr)
        library(stringr)
        library(tibble)
      })
      `%||%` <- function(x, y) {
        if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x)) || !nzchar(x)) {
          y
        } else {
          x
        }
      }
    })
    parallel::clusterExport(cl, c("load_tracking_file"), envir = environment())
    tracks <- parallel::parLapply(cl, seq_along(tracking_files), function(i, files, keep_args, dots) {
      path <- files[[i]]
      message(sprintf("[worker %s] loading %s (%d/%d)", Sys.getpid(), basename(path), i, length(files)))
      do.call(load_tracking_file, c(list(path = path, keep = keep_args), dots))
    }, files = tracking_files, keep_args = keep, dots = list(...))
  }
  names(tracks) <- basename(tracking_files)
  if (bind) {
    bind_rows(tracks)
  } else {
    tracks
  }
}

add_nearest_object_distance <- function(
  tracks,
  frame_col = "frame",
  x_col = "Center_of_the_object_1",
  y_col = "Center_of_the_object_0",
  out_col = "nearest_object_distance"
) {
  stopifnot(all(c(frame_col, x_col, y_col) %in% names(tracks)))

  dt <- as.data.table(tracks)
  dt[, row_id_for_nearest := .I]
  dt[, (out_col) := NA_real_]

  group_cols <- frame_col
  split_keys <- unique(dt[, ..group_cols])
  message("Computing nearest-object distances for ", nrow(split_keys), " frame group(s)")

  for (i in seq_len(nrow(split_keys))) {
    if (i %% 100 == 0 || i == nrow(split_keys)) {
      message(sprintf("nearest-object distances: %d/%d groups", i, nrow(split_keys)))
    }

    idx <- rep(TRUE, nrow(dt))
    for (col in group_cols) {
      idx <- idx & dt[[col]] == split_keys[[col]][[i]]
    }
    idx <- which(idx & is.finite(dt[[x_col]]) & is.finite(dt[[y_col]]))
    if (length(idx) < 2) {
      next
    }

    coords <- as.matrix(dt[idx, c(x_col, y_col), with = FALSE])
    dist_mat <- as.matrix(stats::dist(coords))
    diag(dist_mat) <- Inf
    dt[idx, (out_col) := apply(dist_mat, 1, min)]
  }

  dt[, row_id_for_nearest := NULL]
  as_tibble(dt)
}

preprocess_tracking_file <- function(path, keep = NULL, platemap = NULL, ...) {
  tracks <- load_tracking_file(path, keep = keep, ...)
  meta <- parse_tracking_filename(path)
  if (is.null(platemap)) {
    platemap <- default_platemap_layout()
  }

  annotated <- tracks |>
    mutate(
      well = meta$well[[1]],
      row = meta$row[[1]],
      col = meta$col[[1]],
      position = meta$position[[1]]
    ) |>
    left_join(platemap |> select(well, ploidy, Gemcitabine, condition), by = "well")

  add_nearest_object_distance(annotated)
}

load_preprocessed_tracking_directory <- function(dir_path, keep = NULL, platemap = NULL, cores = max(1L, parallel::detectCores(logical = TRUE) - 1L), bind = TRUE, ...) {
  tracking_files <- sort(list.files(dir_path, pattern = "\\.csv$", full.names = TRUE))
  if (length(tracking_files) == 0) {
    stop("No tracking CSV files found in: ", dir_path, call. = FALSE)
  }
  if (cores < 1L) {
    stop("cores must be >= 1", call. = FALSE)
  }
  message("Found ", length(tracking_files), " tracking CSV files")
  message("Preprocessing on ", cores, " worker(s)")

  if (cores == 1L) {
    tracks <- vector("list", length(tracking_files))
    for (i in seq_along(tracking_files)) {
      message(sprintf("[%d/%d] preprocessing %s", i, length(tracking_files), basename(tracking_files[[i]])))
      tracks[[i]] <- preprocess_tracking_file(tracking_files[[i]], keep = keep, platemap = platemap, ...)
    }
  } else {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(data.table)
        library(dplyr)
        library(readr)
        library(stringr)
        library(tibble)
      })
    })
    parallel::clusterExport(
      cl,
      c(
        "load_tracking_file",
        "parse_tracking_filename",
        "default_platemap_layout",
        "add_nearest_object_distance",
        "preprocess_tracking_file"
      ),
      envir = environment()
    )
    tracks <- parallel::parLapply(cl, seq_along(tracking_files), function(i, files, keep_args, platemap_arg, dots) {
      path <- files[[i]]
      message(sprintf("[worker %s] preprocessing %s (%d/%d)", Sys.getpid(), basename(path), i, length(files)))
      do.call(preprocess_tracking_file, c(list(path = path, keep = keep_args, platemap = platemap_arg), dots))
    }, files = tracking_files, keep_args = keep, platemap_arg = platemap, dots = list(...))
  }

  names(tracks) <- basename(tracking_files)
  if (bind) {
    bind_rows(tracks)
  } else {
    tracks
  }
}

annotate_tracking_data <- function(tracks, platemap = NULL) {
  if (!"source_file" %in% names(tracks)) {
    stop("tracks must include source_file from load_tracking_file()", call. = FALSE)
  }
  meta <- tracks |>
    distinct(source_file) |>
    mutate(
      file_meta = lapply(source_file, parse_tracking_filename)
    ) |>
    tidyr::unnest(file_meta)

  if (is.null(platemap)) {
    platemap <- default_platemap_layout()
  }

  tracks |>
    left_join(meta |> select(source_file, well, row, col, position), by = "source_file") |>
    left_join(platemap |> select(well, ploidy, Gemcitabine, condition), by = "well")
}

summarize_tracking_frames <- function(tracks, image_pixels = 1, area_cols = c("Object_Area_0", "Size_in_pixels_0")) {
  area_col <- intersect(area_cols, names(tracks))
  if (length(area_col) == 0) {
    area_col <- NA_character_
  } else {
    area_col <- area_col[[1]]
  }

  group_cols <- intersect(c("well", "row", "col", "position", "ploidy", "Gemcitabine", "condition"), names(tracks))

  track_lengths <- tracks |>
    filter(!is.na(.data$trackId), .data$trackId >= 0) |>
    group_by(across(all_of(c(group_cols, "trackId")))) |>
    summarize(track_length_frames = n_distinct(.data$frame), .groups = "drop")

  tracks |>
    left_join(track_lengths, by = c(group_cols, "trackId")) |>
    group_by(across(all_of(c(group_cols, "frame")))) |>
    summarize(
      n_cells = n(),
      occupied_pixels = if (!is.na(area_col)) sum(.data[[area_col]], na.rm = TRUE) else NA_real_,
      average_length_current_tracks = mean(.data$track_length_frames[.data$trackId >= 0], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      confluency = if (is.finite(image_pixels) && image_pixels > 0) occupied_pixels / image_pixels else NA_real_
    ) |>
    arrange(across(all_of(c(group_cols, "frame"))))
}

summarize_isolation_frames <- function(
  tracks,
  isolation_multipliers = c(1, 2, 3),
  distance_col = "nearest_object_distance",
  diameter_col = "Diameter_0"
) {
  stopifnot(all(c(distance_col, diameter_col, "frame") %in% names(tracks)))

  group_cols <- intersect(c("well", "row", "col", "position", "ploidy", "Gemcitabine", "condition", "frame"), names(tracks))

  base_summary <- tracks |>
    group_by(across(all_of(group_cols))) |>
    summarize(n_cells = n(), .groups = "drop")

  isolation_counts <- lapply(isolation_multipliers, function(multiplier) {
    tracks |>
      group_by(across(all_of(group_cols))) |>
      summarize(
        threshold_multiplier = multiplier,
        threshold = paste0("distance > ", multiplier, "x diameter"),
        n_isolated = sum(.data[[distance_col]] > multiplier * .data[[diameter_col]], na.rm = TRUE),
        .groups = "drop"
      )
  }) |>
    bind_rows()

  isolation_counts |>
    left_join(base_summary, by = group_cols) |>
    mutate(fraction_isolated = .data$n_isolated / pmax(.data$n_cells, 1L)) |>
    arrange(across(all_of(c(group_cols, "threshold_multiplier"))))
}
