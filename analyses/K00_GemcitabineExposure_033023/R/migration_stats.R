suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
  library(tibble)
})

default_track_group_cols <- function(tracks) {
  intersect(c("well", "row", "col", "position", "ploidy", "Gemcitabine", "condition", "trackId"), names(tracks))
}

default_summary_group_cols <- function(tracks) {
  intersect(c("well", "row", "col", "position", "ploidy", "Gemcitabine", "condition"), names(tracks))
}

filter_isolated_contiguous_tracks <- function(
  tracks,
  isolation_multiplier = 2,
  min_contiguous_frames = 3,
  distance_col = "nearest_object_distance",
  diameter_col = "Diameter_0",
  frame_col = "frame",
  track_col = "trackId",
  x_col = "Center_of_the_object_1",
  y_col = "Center_of_the_object_0"
) {
  required <- c(distance_col, diameter_col, frame_col, track_col, x_col, y_col)
  stopifnot(all(required %in% names(tracks)))

  meta_cols <- default_summary_group_cols(tracks)
  track_group_cols <- c(meta_cols, track_col)

  isolated <- tracks |>
    filter(
      !is.na(.data[[track_col]]),
      .data[[track_col]] >= 0,
      is.finite(.data[[distance_col]]),
      is.finite(.data[[diameter_col]]),
      .data[[distance_col]] > isolation_multiplier * .data[[diameter_col]]
    ) |>
    arrange(across(all_of(track_group_cols)), .data[[frame_col]]) |>
    group_by(across(all_of(track_group_cols))) |>
    mutate(
      starts_segment = row_number() == 1L | .data[[frame_col]] != lag(.data[[frame_col]]) + 1L,
      contiguous_segment_index = cumsum(starts_segment),
      migration_track_id = paste(
        well,
        position,
        .data[[track_col]],
        contiguous_segment_index,
        sep = "_"
      )
    ) |>
    ungroup()

  segment_lengths <- isolated |>
    group_by(.data$migration_track_id) |>
    summarize(contiguous_track_length_frames = n_distinct(.data[[frame_col]]), .groups = "drop")

  isolated |>
    left_join(segment_lengths, by = "migration_track_id") |>
    filter(.data$contiguous_track_length_frames >= min_contiguous_frames) |>
    arrange(.data$migration_track_id, .data[[frame_col]])
}

summarize_retention_by_frame <- function(raw_tracks, filtered_tracks, frame_col = "frame", track_col = "trackId") {
  meta_cols <- default_summary_group_cols(raw_tracks)
  group_cols <- c(meta_cols, frame_col)

  raw_summary <- raw_tracks |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      n_objects_input = n(),
      n_linked_input = sum(!is.na(.data[[track_col]]) & .data[[track_col]] >= 0),
      .groups = "drop"
    )

  filtered_summary <- filtered_tracks |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      n_objects_retained = n(),
      n_tracks_retained = n_distinct(.data$migration_track_id),
      average_retained_track_length_frames = mean(.data$contiguous_track_length_frames, na.rm = TRUE),
      .groups = "drop"
    )

  raw_summary |>
    left_join(filtered_summary, by = group_cols) |>
    mutate(
      n_objects_retained = coalesce(.data$n_objects_retained, 0L),
      n_tracks_retained = coalesce(.data$n_tracks_retained, 0L),
      fraction_objects_retained = .data$n_objects_retained / pmax(.data$n_objects_input, 1L)
    )
}

compute_steps <- function(
  tracks,
  frame_col = "frame",
  x_col = "Center_of_the_object_1",
  y_col = "Center_of_the_object_0",
  segment_col = "migration_track_id"
) {
  meta_cols <- default_summary_group_cols(tracks)

  tracks |>
    arrange(.data[[segment_col]], .data[[frame_col]]) |>
    group_by(.data[[segment_col]]) |>
    mutate(
      next_frame = lead(.data[[frame_col]]),
      dx = lead(.data[[x_col]]) - .data[[x_col]],
      dy = lead(.data[[y_col]]) - .data[[y_col]],
      step_frames = next_frame - .data[[frame_col]],
      step_distance_px = sqrt(dx^2 + dy^2),
      speed_px_frame = step_distance_px / step_frames
    ) |>
    ungroup() |>
    filter(step_frames == 1L, is.finite(step_distance_px)) |>
    select(all_of(c(meta_cols, segment_col, frame_col, "next_frame", "dx", "dy", "step_distance_px", "speed_px_frame")))
}

summarize_steps_by_frame <- function(steps) {
  group_cols <- c(default_summary_group_cols(steps), "frame")

  steps |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      n_steps = n(),
      mean_step_distance_px = mean(.data$step_distance_px, na.rm = TRUE),
      median_step_distance_px = median(.data$step_distance_px, na.rm = TRUE),
      mean_speed_px_frame = mean(.data$speed_px_frame, na.rm = TRUE),
      median_speed_px_frame = median(.data$speed_px_frame, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_turning_angles <- function(steps, segment_col = "migration_track_id") {
  meta_cols <- default_summary_group_cols(steps)

  steps |>
    arrange(.data[[segment_col]], .data$frame) |>
    group_by(.data[[segment_col]]) |>
    mutate(
      next_step_frame = lead(.data$frame),
      dx_next = lead(.data$dx),
      dy_next = lead(.data$dy),
      next_step_distance_px = lead(.data$step_distance_px),
      cos_turning_angle = (.data$dx * dx_next + .data$dy * dy_next) /
        pmax(.data$step_distance_px * next_step_distance_px, 1e-9),
      cos_turning_angle = pmax(pmin(cos_turning_angle, 1), -1),
      turning_angle_degrees = acos(cos_turning_angle) * 180 / pi
    ) |>
    ungroup() |>
    filter(next_step_frame == frame + 1L, is.finite(cos_turning_angle)) |>
    select(all_of(c(meta_cols, segment_col, "frame", "cos_turning_angle", "turning_angle_degrees")))
}

summarize_turning_by_frame <- function(turning) {
  group_cols <- c(default_summary_group_cols(turning), "frame")

  turning |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      n_turns = n(),
      mean_cos_turning_angle = mean(.data$cos_turning_angle, na.rm = TRUE),
      mean_turning_angle_degrees = mean(.data$turning_angle_degrees, na.rm = TRUE),
      median_turning_angle_degrees = median(.data$turning_angle_degrees, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_msd_by_lag <- function(
  tracks,
  max_lag_frames = 10,
  frame_col = "frame",
  x_col = "Center_of_the_object_1",
  y_col = "Center_of_the_object_0",
  segment_col = "migration_track_id"
) {
  meta_cols <- default_summary_group_cols(tracks)

  tracks |>
    arrange(.data[[segment_col]], .data[[frame_col]]) |>
    group_by(.data[[segment_col]]) |>
    group_modify(~ {
      df <- .x |> arrange(.data[[frame_col]])
      frames <- df[[frame_col]]
      x <- df[[x_col]]
      y <- df[[y_col]]

      bind_rows(lapply(seq_len(max_lag_frames), function(lag_frame) {
        j <- match(frames + lag_frame, frames)
        ok <- !is.na(j)
        if (!any(ok)) {
          return(NULL)
        }
        tibble(
          frame = frames[ok],
          lag_frames = lag_frame,
          squared_displacement_px2 = (x[j[ok]] - x[ok])^2 + (y[j[ok]] - y[ok])^2,
          displacement_px = sqrt(squared_displacement_px2)
        )
      }))
    }) |>
    ungroup() |>
    left_join(tracks |> select(all_of(c(meta_cols, segment_col))) |> distinct(), by = segment_col)
}

summarize_msd_by_frame_lag <- function(msd) {
  group_cols <- c(default_summary_group_cols(msd), "frame", "lag_frames")

  msd |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      n_displacements = n(),
      mean_msd_px2 = mean(.data$squared_displacement_px2, na.rm = TRUE),
      median_msd_px2 = median(.data$squared_displacement_px2, na.rm = TRUE),
      mean_displacement_px = mean(.data$displacement_px, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_velocity_autocorrelation <- function(steps, max_lag_frames = 10, segment_col = "migration_track_id") {
  meta_cols <- default_summary_group_cols(steps)

  steps |>
    arrange(.data[[segment_col]], .data$frame) |>
    group_by(.data[[segment_col]]) |>
    group_modify(~ {
      df <- .x |> arrange(.data$frame)
      frames <- df$frame

      bind_rows(lapply(seq_len(max_lag_frames), function(lag_frame) {
        j <- match(frames + lag_frame, frames)
        ok <- !is.na(j)
        if (!any(ok)) {
          return(NULL)
        }
        dot_product <- df$dx[ok] * df$dx[j[ok]] + df$dy[ok] * df$dy[j[ok]]
        norm_product <- pmax(df$step_distance_px[ok] * df$step_distance_px[j[ok]], 1e-9)
        tibble(
          frame = frames[ok],
          lag_frames = lag_frame,
          velocity_dot_product = dot_product,
          velocity_cos_autocorrelation = dot_product / norm_product
        )
      }))
    }) |>
    ungroup() |>
    left_join(steps |> select(all_of(c(meta_cols, segment_col))) |> distinct(), by = segment_col)
}

summarize_autocorrelation_by_frame_lag <- function(autocorrelation) {
  group_cols <- c(default_summary_group_cols(autocorrelation), "frame", "lag_frames")

  autocorrelation |>
    group_by(across(all_of(group_cols))) |>
    summarize(
      n_pairs = n(),
      mean_velocity_dot_product = mean(.data$velocity_dot_product, na.rm = TRUE),
      mean_velocity_cos_autocorrelation = mean(.data$velocity_cos_autocorrelation, na.rm = TRUE),
      median_velocity_cos_autocorrelation = median(.data$velocity_cos_autocorrelation, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_migration_summaries_for_group <- function(group_tracks, max_lag_frames = 10) {
  steps <- compute_steps(group_tracks)
  turning <- compute_turning_angles(steps)
  msd <- compute_msd_by_lag(group_tracks, max_lag_frames = max_lag_frames)
  autocorrelation <- compute_velocity_autocorrelation(steps, max_lag_frames = max_lag_frames)

  list(
    step_summary = summarize_steps_by_frame(steps),
    turning_summary = summarize_turning_by_frame(turning),
    msd_summary = summarize_msd_by_frame_lag(msd),
    autocorrelation_summary = summarize_autocorrelation_by_frame_lag(autocorrelation)
  )
}

compute_migration_summaries_parallel <- function(tracks, max_lag_frames = 10, cores = 1L) {
  split_cols <- default_summary_group_cols(tracks)
  grouped_tracks <- tracks |>
    group_by(across(all_of(split_cols))) |>
    group_split()

  message("Computing migration summaries for ", length(grouped_tracks), " group(s) on ", cores, " worker(s)")

  if (cores <= 1L) {
    results <- lapply(grouped_tracks, compute_migration_summaries_for_group, max_lag_frames = max_lag_frames)
  } else {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(dplyr)
        library(tibble)
      })
    })
    parallel::clusterExport(
      cl,
      c(
        "default_summary_group_cols",
        "compute_steps",
        "summarize_steps_by_frame",
        "compute_turning_angles",
        "summarize_turning_by_frame",
        "compute_msd_by_lag",
        "summarize_msd_by_frame_lag",
        "compute_velocity_autocorrelation",
        "summarize_autocorrelation_by_frame_lag",
        "compute_migration_summaries_for_group"
      ),
      envir = environment()
    )
    results <- parallel::parLapply(cl, grouped_tracks, compute_migration_summaries_for_group, max_lag_frames = max_lag_frames)
  }

  list(
    step_summary = bind_rows(lapply(results, `[[`, "step_summary")),
    turning_summary = bind_rows(lapply(results, `[[`, "turning_summary")),
    msd_summary = bind_rows(lapply(results, `[[`, "msd_summary")),
    autocorrelation_summary = bind_rows(lapply(results, `[[`, "autocorrelation_summary"))
  )
}
