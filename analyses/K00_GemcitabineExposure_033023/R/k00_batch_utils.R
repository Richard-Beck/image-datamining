`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) {
    y
  } else {
    x
  }
}

parse_cli_args <- function(args, usage = NULL) {
  if (any(args %in% c("--help", "-h", "--help=TRUE"))) {
    if (!is.null(usage)) {
      cat(usage)
    }
    quit(save = "no", status = 0)
  }

  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) {
      stop("Arguments must be passed as --name=value: ", arg, call. = FALSE)
    }
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    out[[parts[1]]] <- paste(parts[-1], collapse = "=")
  }
  out
}

as_flag <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

repo_root_from_script <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
    return(normalizePath(file.path(dirname(script_path), "../../.."), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

safe_file_stem <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  gsub("[^A-Za-z0-9._-]+", "_", stem)
}

well_from_path <- function(path) {
  stringr::str_extract(basename(path), "^[A-H][0-9]{1,2}")
}

site_from_path <- function(path) {
  as.integer(stringr::str_match(basename(path), "^[A-H][0-9]{1,2}_([0-9]+)_")[, 2])
}

condition_id_from_parts <- function(ploidy, gemcitabine_nm) {
  dose <- format(round(gemcitabine_nm, 6), scientific = FALSE, trim = TRUE)
  dose <- sub("0+$", "", dose)
  dose <- sub("\\.$", "", dose)
  dose <- gsub("\\.", "p", dose)
  paste0(ploidy, "_gem_", dose, "nM")
}

default_platemap_layout <- function() {
  rows <- LETTERS[1:8]
  cols <- 1:12
  doses <- c(0, 3.125, 6.25, 12.5, 25, 50, 100, 200, 400, 800)
  expand.grid(plate_row = rows, plate_col = cols, stringsAsFactors = FALSE) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      well = paste0(.data$plate_row, .data$plate_col),
      ploidy = dplyr::case_when(
        .data$plate_col %in% c(1L, 12L) ~ NA_character_,
        .data$plate_row %in% LETTERS[1:4] ~ "2N",
        .data$plate_row %in% LETTERS[5:8] ~ "4N",
        TRUE ~ NA_character_
      ),
      gemcitabine_nm = dplyr::if_else(
        .data$plate_col %in% 2:11,
        doses[.data$plate_col - 1L],
        NA_real_
      ),
      condition_id = dplyr::if_else(
        !is.na(.data$ploidy) & !is.na(.data$gemcitabine_nm),
        condition_id_from_parts(.data$ploidy, .data$gemcitabine_nm),
        NA_character_
      )
    ) |>
    dplyr::arrange(.data$plate_row, .data$plate_col)
}

parse_platemap_conditions <- function(platemap_xlsx) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    warning("readxl is unavailable; using the expected K00 plate layout fallback.", call. = FALSE)
    return(default_platemap_layout())
  }

  raw <- readxl::read_excel(platemap_xlsx, col_names = FALSE, .name_repair = "minimal")
  raw_mat <- as.matrix(raw)
  header <- as.character(raw_mat[1, ])
  plate_cols <- suppressWarnings(as.integer(header))
  row_labels <- as.character(raw[[1]])

  rows <- which(row_labels %in% LETTERS[1:8])
  cols <- which(!is.na(plate_cols) & plate_cols %in% 1:12)

  cells <- expand.grid(row_index = rows, col_index = cols, stringsAsFactors = FALSE) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      plate_row = row_labels[.data$row_index],
      plate_col = plate_cols[.data$col_index],
      well = paste0(.data$plate_row, .data$plate_col),
      treatment = as.character(raw_mat[cbind(.data$row_index, .data$col_index)])
    )

  cells |>
    dplyr::mutate(
      treatment = dplyr::na_if(trimws(.data$treatment), ""),
      ploidy = stringr::str_match(.data$treatment, "^(2N|4N):")[, 2],
      gemcitabine_nm = suppressWarnings(as.numeric(stringr::str_match(
        .data$treatment,
        ":\\s*([0-9.]+)\\s*nM"
      )[, 2])),
      condition_id = dplyr::if_else(
        !is.na(.data$ploidy) & !is.na(.data$gemcitabine_nm),
        condition_id_from_parts(.data$ploidy, .data$gemcitabine_nm),
        NA_character_
      )
    ) |>
    dplyr::select("well", "plate_row", "plate_col", "ploidy", "gemcitabine_nm", "condition_id") |>
    dplyr::arrange(.data$plate_row, .data$plate_col)
}

append_delim_locked <- function(row, out_path, delim = ",", quiet = FALSE) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  lock_dir <- paste0(out_path, ".lock")
  deadline <- Sys.time() + 600
  while (!dir.create(lock_dir, showWarnings = FALSE)) {
    if (Sys.time() > deadline) {
      stop("Timed out waiting for lock: ", lock_dir, call. = FALSE)
    }
    Sys.sleep(runif(1, 0.05, 0.25))
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)

  write_header <- !file.exists(out_path) || file.info(out_path)$size == 0
  if (!write_header) {
    existing_header <- strsplit(readLines(out_path, n = 1), delim, fixed = TRUE)[[1]]
    if (!identical(existing_header, names(row))) {
      existing <- if (identical(delim, ",")) {
        readr::read_csv(out_path, show_col_types = FALSE)
      } else {
        readr::read_delim(out_path, delim = delim, show_col_types = FALSE)
      }
      all_cols <- union(names(existing), names(row))
      for (missing_col in setdiff(all_cols, names(existing))) {
        existing[[missing_col]] <- NA
      }
      for (missing_col in setdiff(all_cols, names(row))) {
        row[[missing_col]] <- NA
      }
      combined <- dplyr::bind_rows(
        existing |> dplyr::select(dplyr::all_of(all_cols)),
        row |> dplyr::select(dplyr::all_of(all_cols))
      )
      if (identical(delim, ",")) {
        readr::write_csv(combined, out_path)
      } else {
        readr::write_delim(combined, out_path, delim = delim)
      }
      if (!quiet) {
        message("Upgraded and appended ", out_path)
      }
      return(invisible(NULL))
    }
  }

  if (identical(delim, ",")) {
    readr::write_csv(row, out_path, append = !write_header, col_names = write_header)
  } else {
    readr::write_delim(row, out_path, delim = delim, append = !write_header, col_names = write_header)
  }
  if (!quiet) {
    message("Appended ", out_path)
  }
}
