#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

REPO_ROOT <- normalizePath(getwd(), mustWork = TRUE)
ANALYSIS_DIR <- file.path(REPO_ROOT, "analyses/K00_GemcitabineExposure_033023")
DATA_DIR <- file.path(ANALYSIS_DIR, "data")
FITS_CSV <- file.path(DATA_DIR, "ou_tracking_fits.csv")
TRACKS_RDS <- file.path(DATA_DIR, "tracking_data_isolated_2x_min3.rds")
OUT_DIR <- file.path(ANALYSIS_DIR, "overlay_checks/ou_fit_sanity")
IMAGES_DIR <- Sys.getenv(
  "K00_IMAGES_DIR",
  "/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Images_40Frames"
)

N_TRACKS_PER_GROUP <- as.integer(Sys.getenv("N_TRACKS_PER_GROUP", "5"))
MAX_FRAMES_PER_TRACK <- as.integer(Sys.getenv("MAX_FRAMES_PER_TRACK", "10"))
PRIMARY_PARAMETER <- Sys.getenv("PRIMARY_PARAMETER", "effective_diffusivity")

source(file.path(ANALYSIS_DIR, "R/overlays.R"))

if (!file.exists(FITS_CSV)) {
  stop("Missing OU fits CSV: ", FITS_CSV, call. = FALSE)
}
if (!file.exists(TRACKS_RDS)) {
  stop("Missing isolated tracks RDS: ", TRACKS_RDS, call. = FALSE)
}
if (!dir.exists(IMAGES_DIR)) {
  stop("Missing registered images directory: ", IMAGES_DIR, call. = FALSE)
}

safe_log2_ratio <- function(numerator, denominator) {
  ifelse(is.finite(numerator) & is.finite(denominator) & numerator > 0 & denominator > 0,
         log2(numerator / denominator), NA_real_)
}

best_condition_fits <- function(fits) {
  if (!"start_id" %in% names(fits)) {
    fits$start_id <- NA_integer_
  }
  fits |>
    group_by(ploidy, Gemcitabine) |>
    slice_max(order_by = log_likelihood, n = 1, with_ties = FALSE) |>
    ungroup()
}

choose_doses <- function(best_fits, parameter = PRIMARY_PARAMETER) {
  if (!parameter %in% names(best_fits)) {
    stop("Parameter not found in fits: ", parameter, call. = FALSE)
  }

  contrast <- best_fits |>
    select(ploidy, Gemcitabine, value = all_of(parameter)) |>
    pivot_wider(names_from = ploidy, values_from = value) |>
    mutate(
      parameter = parameter,
      log2_4N_over_2N = safe_log2_ratio(`4N`, `2N`),
      abs_log2_4N_over_2N = abs(log2_4N_over_2N)
    ) |>
    arrange(desc(abs_log2_4N_over_2N))

  max_dose <- contrast |>
    filter(is.finite(abs_log2_4N_over_2N), Gemcitabine != 0) |>
    slice_head(n = 1) |>
    pull(Gemcitabine)
  if (length(max_dose) == 0L) {
    max_dose <- contrast |>
      filter(is.finite(abs_log2_4N_over_2N)) |>
      slice_head(n = 1) |>
      pull(Gemcitabine)
  }

  selected <- unique(c(0, max_dose))
  list(contrast = contrast, doses = selected)
}

track_summary <- function(tracks) {
  tracks |>
    arrange(migration_track_id, frame) |>
    group_by(migration_track_id, ploidy, Gemcitabine, well, position) |>
    summarize(
      n_frames = n_distinct(frame),
      first_frame = min(frame, na.rm = TRUE),
      last_frame = max(frame, na.rm = TRUE),
      net_displacement_px = sqrt(
        (last(Center_of_the_object_1) - first(Center_of_the_object_1))^2 +
          (last(Center_of_the_object_0) - first(Center_of_the_object_0))^2
      ),
      path_length_px = sum(sqrt(
        diff(Center_of_the_object_1)^2 + diff(Center_of_the_object_0)^2
      ), na.rm = TRUE),
      .groups = "drop"
    ) |>
    filter(is.finite(path_length_px), n_frames >= 3)
}

select_tracks_for_group <- function(summary, dose, ploidy_value, n_tracks = N_TRACKS_PER_GROUP) {
  key_cols <- c("migration_track_id", "ploidy", "Gemcitabine", "well", "position")
  candidates <- summary |>
    filter(.data$Gemcitabine == dose, .data$ploidy == ploidy_value) |>
    mutate(path_rank = percent_rank(path_length_px))

  if (nrow(candidates) == 0L) {
    return(candidates)
  }

  high <- candidates |>
    arrange(desc(path_length_px), desc(n_frames)) |>
    slice_head(n = 2) |>
    mutate(category = "high_path")
  median_tracks <- candidates |>
    anti_join(high |> select(all_of(key_cols)), by = key_cols) |>
    arrange(abs(path_rank - 0.5), desc(n_frames)) |>
    slice_head(n = 2) |>
    mutate(category = "median_path")
  low <- candidates |>
    anti_join(bind_rows(high, median_tracks) |> select(all_of(key_cols)), by = key_cols) |>
    arrange(path_length_px, desc(n_frames)) |>
    slice_head(n = 1) |>
    mutate(category = "low_path")
  fallback <- candidates |>
    anti_join(bind_rows(high, median_tracks, low) |> select(all_of(key_cols)), by = key_cols) |>
    arrange(desc(n_frames), desc(path_length_px)) |>
    mutate(category = "fallback")

  bind_rows(high, median_tracks, low, fallback) |>
    distinct(across(all_of(key_cols)), .keep_all = TRUE) |>
    slice_head(n = min(n_tracks, nrow(candidates)))
}

render_track <- function(track_rows, selected_row, dose_label) {
  frames <- sort(unique(as.integer(track_rows$frame)))
  if (length(frames) > MAX_FRAMES_PER_TRACK) {
    center <- ceiling(length(frames) / 2)
    half_window <- floor(MAX_FRAMES_PER_TRACK / 2)
    keep_idx <- seq(
      max(1, center - half_window),
      min(length(frames), center + half_window)
    )
    frames <- frames[keep_idx][seq_len(min(length(keep_idx), MAX_FRAMES_PER_TRACK))]
  }

  track_id_safe <- gsub("[^A-Za-z0-9_=-]+", "_", selected_row$migration_track_id)
  prefix <- paste0(
    "dose_", dose_label,
    "_", selected_row$ploidy,
    "_", selected_row$category,
    "_", track_id_safe
  )
  out_dir <- file.path(OUT_DIR, paste0("dose_", dose_label))

  rendered <- tracking_overlay_png_sequence(
    tracks = track_rows,
    images_dir = IMAGES_DIR,
    out_dir = out_dir,
    output_prefix = prefix,
    frames = frames,
    crop = "tracks",
    crop_buffer_px = 40,
    arrow_mode = "current_to_next"
  )

  out_gif <- file.path(out_dir, paste0(prefix, ".gif"))
  gif_status <- tryCatch(
    {
      png_sequence_gif(rendered$frame_paths, out_gif = out_gif, fps = 2, cleanup_frames = TRUE)
      out_gif
    },
    error = function(e) {
      message("GIF creation failed for ", selected_row$migration_track_id, ": ", conditionMessage(e))
      NA_character_
    }
  )

  tibble(
    migration_track_id = selected_row$migration_track_id,
    ploidy = selected_row$ploidy,
    Gemcitabine = selected_row$Gemcitabine,
    category = selected_row$category,
    well = selected_row$well,
    position = selected_row$position,
    n_frames = selected_row$n_frames,
    path_length_px = selected_row$path_length_px,
    net_displacement_px = selected_row$net_displacement_px,
    rendered_frames = paste(frames, collapse = ";"),
    gif_path = gif_status
  )
}

fits <- read_csv(FITS_CSV, show_col_types = FALSE)
best_fits <- best_condition_fits(fits)
dose_choice <- choose_doses(best_fits)
selected_doses <- dose_choice$doses

tracks <- readRDS(TRACKS_RDS)
summary <- track_summary(tracks)

selected_tracks <- bind_rows(lapply(selected_doses, function(dose) {
  bind_rows(
    select_tracks_for_group(summary, dose = dose, ploidy_value = "2N"),
    select_tracks_for_group(summary, dose = dose, ploidy_value = "4N")
  )
})) |>
  distinct(ploidy, Gemcitabine, well, position, migration_track_id, .keep_all = TRUE)

if (nrow(selected_tracks) == 0L) {
  stop("No tracks selected for overlay rendering.", call. = FALSE)
}

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(dose_choice$contrast, file.path(OUT_DIR, "ou_fit_2N_4N_contrasts.csv"))
write_csv(selected_tracks, file.path(OUT_DIR, "selected_tracks.csv"))

message("Selected doses: ", paste(selected_doses, collapse = ", "))
message("Rendering ", nrow(selected_tracks), " track overlays into ", OUT_DIR)

rendered <- bind_rows(lapply(seq_len(nrow(selected_tracks)), function(i) {
  selected_row <- selected_tracks[i, , drop = FALSE]
  track_rows <- tracks |>
    filter(
      .data$migration_track_id == selected_row$migration_track_id[[1]],
      .data$ploidy == selected_row$ploidy[[1]],
      .data$Gemcitabine == selected_row$Gemcitabine[[1]],
      .data$well == selected_row$well[[1]],
      .data$position == selected_row$position[[1]]
    )
  dose_label <- gsub("\\.", "p", format(selected_row$Gemcitabine[[1]], trim = TRUE, scientific = FALSE))
  render_track(track_rows, selected_row, dose_label)
}))

write_csv(rendered, file.path(OUT_DIR, "rendered_overlays.csv"))
message("Wrote ", file.path(OUT_DIR, "selected_tracks.csv"))
message("Wrote ", file.path(OUT_DIR, "rendered_overlays.csv"))
