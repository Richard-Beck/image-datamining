suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

if (!requireNamespace("Rcpp", quietly = TRUE)) {
  stop("Rcpp is required for the OU velocity likelihood.", call. = FALSE)
}

Rcpp::cppFunction(
  includes = "#include <cmath>",
  code = '
double ou_loglik_cpp(Rcpp::NumericVector log_params, Rcpp::List segments) {
  auto kalman_dim_loglik_scalar = [](
    const Rcpp::NumericVector& y,
    const Rcpp::NumericVector& dt,
    const double tau,
    const double velocity_scale,
    const double obs_noise
  ) {
    const int n = y.size();
    if (n < 2) {
      return 0.0;
    }

    double m_x = y[0];
    double m_v = 0.0;
    double P_xx = obs_noise * obs_noise;
    double P_xv = 0.0;
    double P_vv = velocity_scale * velocity_scale;
    const double R = obs_noise * obs_noise;
    const double velocity_var = velocity_scale * velocity_scale;
    double loglik = 0.0;

    for (int i = 1; i < n; ++i) {
      const double dt_i = dt[i - 1];
      const double a = dt_i / tau;
      const double phi = std::exp(-a);
      const double one_minus_phi = -std::expm1(-a);
      const double one_minus_phi2 = -std::expm1(-2.0 * a);
      const double f12 = tau * one_minus_phi;

      double q22 = velocity_var * one_minus_phi2;
      const double q12 = velocity_var * tau * one_minus_phi * one_minus_phi;
      double q11_unit;
      if (a < 1e-4) {
        q11_unit = a * a * a / 3.0 - a * a * a * a / 4.0 + 7.0 * a * a * a * a * a / 60.0;
      } else {
        q11_unit = a - 2.0 * one_minus_phi + 0.5 * one_minus_phi2;
      }
      double q11 = 2.0 * velocity_var * tau * tau * q11_unit;

      if (q11 < 0.0) {
        q11 = 0.0;
      }
      if (q22 < 0.0) {
        q22 = 0.0;
      }
      if (q22 > 0.0 && q11 * q22 < q12 * q12) {
        q11 = q12 * q12 / q22 + std::numeric_limits<double>::epsilon();
      }

      const double m_x_pred = m_x + f12 * m_v;
      const double m_v_pred = phi * m_v;

      const double P_xx_pred = P_xx + 2.0 * f12 * P_xv + f12 * f12 * P_vv + q11;
      const double P_xv_pred = phi * (P_xv + f12 * P_vv) + q12;
      const double P_vv_pred = phi * phi * P_vv + q22;

      const double innovation = y[i] - m_x_pred;
      const double S = P_xx_pred + R;
      if (!R_finite(S) || S <= 0.0) {
        return R_NegInf;
      }
      loglik += -0.5 * (std::log(2.0 * M_PI) + std::log(S) + innovation * innovation / S);

      const double K_x = P_xx_pred / S;
      const double K_v = P_xv_pred / S;
      m_x = m_x_pred + K_x * innovation;
      m_v = m_v_pred + K_v * innovation;

      const double one_minus_Kx = 1.0 - K_x;
      P_xx = one_minus_Kx * one_minus_Kx * P_xx_pred + K_x * K_x * R;
      P_xv = one_minus_Kx * (P_xv_pred - K_v * P_xx_pred) + K_x * K_v * R;
      P_vv = K_v * K_v * P_xx_pred - 2.0 * K_v * P_xv_pred + P_vv_pred + K_v * K_v * R;
    }

    return loglik;
  };

  const double tau = std::exp(log_params[0]);
  const double velocity_scale = std::exp(log_params[1]);
  const double obs_noise = std::exp(log_params[2]);
  if (!R_finite(tau) || !R_finite(velocity_scale) || !R_finite(obs_noise)) {
    return R_NegInf;
  }

  double total = 0.0;
  const int n_segments = segments.size();
  for (int i = 0; i < n_segments; ++i) {
    Rcpp::List seg = segments[i];
    Rcpp::NumericMatrix y = seg["y"];
    Rcpp::NumericVector dt = seg["dt"];
    total += kalman_dim_loglik_scalar(y(Rcpp::_, 0), dt, tau, velocity_scale, obs_noise);
    total += kalman_dim_loglik_scalar(y(Rcpp::_, 1), dt, tau, velocity_scale, obs_noise);
  }

  return total;
}
')

resolve_track_paths <- function(track_paths, manifest_path, repo_root) {
  vapply(track_paths, function(path) {
    if (grepl("^/", path)) {
      return(normalizePath(path, mustWork = TRUE))
    }

    manifest_relative <- file.path(dirname(manifest_path), path)
    if (file.exists(manifest_relative)) {
      return(normalizePath(manifest_relative, mustWork = TRUE))
    }

    repo_relative <- file.path(repo_root, path)
    normalizePath(repo_relative, mustWork = TRUE)
  }, character(1))
}

load_isolated_track_rds <- function(track_paths, target_wells = NULL) {
  tracks <- bind_rows(lapply(track_paths, function(path) {
    readRDS(path) |>
      mutate(
        source_rds = path,
        file = as.character(file),
        well_from_file = str_extract(file, "^[A-H][0-9]{1,2}"),
        well = dplyr::coalesce(as.character(well), well_from_file),
        site = as.integer(site),
        split_track_id = paste(file, split_track_id, sep = "::")
      )
  })) |>
    filter(!is.na(frame), !is.na(x), !is.na(y), !is.na(split_track_id)) |>
    mutate(
      frame = as.numeric(frame),
      x = as.numeric(x),
      y = as.numeric(y)
    )

  if (!is.null(target_wells)) {
    tracks <- tracks |> filter(.data$well_from_file %in% target_wells | .data$well %in% target_wells)
  }
  tracks
}

split_ou_segments <- function(tracks, frame_interval) {
  raw_segments <- tracks |>
    arrange(split_track_id, frame) |>
    group_by(split_track_id) |>
    summarize(
      frame = list(frame),
      x = list(x),
      y = list(y),
      .groups = "drop"
    )

  segments <- list()
  total_track_time <- 0
  for (i in seq_len(nrow(raw_segments))) {
    frames <- raw_segments$frame[[i]]
    coords <- cbind(raw_segments$x[[i]], raw_segments$y[[i]])
    keep <- is.finite(frames) & is.finite(coords[, 1]) & is.finite(coords[, 2])
    frames <- frames[keep]
    coords <- coords[keep, , drop = FALSE]
    if (length(frames) < 2) {
      next
    }
    ord <- order(frames)
    frames <- frames[ord]
    coords <- coords[ord, , drop = FALSE]
    unique_keep <- !duplicated(frames)
    frames <- frames[unique_keep]
    coords <- coords[unique_keep, , drop = FALSE]
    if (length(frames) < 2) {
      next
    }
    dt <- diff(frames) * frame_interval
    valid <- is.finite(dt) & dt > 0
    if (!all(valid)) {
      next
    }
    total_track_time <- total_track_time + sum(dt)
    segments[[length(segments) + 1]] <- list(y = coords, dt = dt)
  }

  attr(segments, "total_track_time") <- total_track_time
  segments
}

ou_transition_r <- function(dt, tau, velocity_scale) {
  a <- dt / tau
  phi <- exp(-a)
  one_minus_phi <- -expm1(-a)
  one_minus_phi2 <- -expm1(-2 * a)
  f12 <- tau * one_minus_phi
  q22 <- velocity_scale^2 * one_minus_phi2
  q12 <- velocity_scale^2 * tau * one_minus_phi^2
  q11_unit <- if (a < 1e-4) {
    a^3 / 3 - a^4 / 4 + 7 * a^5 / 60
  } else {
    a - 2 * one_minus_phi + 0.5 * one_minus_phi2
  }
  q11 <- 2 * velocity_scale^2 * tau^2 * q11_unit

  q11 <- max(q11, 0)
  q22 <- max(q22, 0)
  if (q22 > 0 && q11 * q22 < q12^2) {
    q11 <- q12^2 / q22 + .Machine$double.eps
  }

  list(
    F = matrix(c(1, 0, f12, phi), nrow = 2),
    Q = matrix(c(q11, q12, q12, q22), nrow = 2)
  )
}

kalman_dim_loglik_r <- function(y, dt, tau, velocity_scale, obs_noise) {
  n <- length(y)
  if (n < 2) {
    return(0)
  }

  m <- c(y[1], 0)
  P <- matrix(c(obs_noise^2, 0, 0, velocity_scale^2), nrow = 2)
  loglik <- 0
  R <- obs_noise^2
  H <- matrix(c(1, 0), nrow = 1)
  I <- diag(2)

  for (i in 2:n) {
    tr <- ou_transition_r(dt[i - 1], tau, velocity_scale)
    m <- drop(tr$F %*% m)
    P <- tr$F %*% P %*% t(tr$F) + tr$Q
    P <- (P + t(P)) / 2

    innovation <- y[i] - m[1]
    S <- P[1, 1] + R
    if (!is.finite(S) || S <= 0) {
      return(-Inf)
    }
    loglik <- loglik - 0.5 * (log(2 * pi) + log(S) + innovation^2 / S)

    K <- P[, 1] / S
    m <- m + K * innovation
    KH <- K %*% H
    P <- (I - KH) %*% P %*% t(I - KH) + tcrossprod(K) * R
    P <- (P + t(P)) / 2
  }

  loglik
}

ou_loglik_r <- function(log_params, segments) {
  tau <- exp(log_params[1])
  velocity_scale <- exp(log_params[2])
  obs_noise <- exp(log_params[3])
  if (!all(is.finite(c(tau, velocity_scale, obs_noise)))) {
    return(-Inf)
  }

  sum(vapply(segments, function(seg) {
    kalman_dim_loglik_r(seg$y[, 1], seg$dt, tau, velocity_scale, obs_noise) +
      kalman_dim_loglik_r(seg$y[, 2], seg$dt, tau, velocity_scale, obs_noise)
  }, numeric(1)))
}

ou_loglik <- function(log_params, segments) {
  ou_loglik_cpp(log_params, segments)
}

initial_ou_params <- function(segments) {
  all_dt <- unlist(lapply(segments, `[[`, "dt"), use.names = FALSE)
  displacements <- do.call(rbind, lapply(segments, function(seg) {
    diff(seg$y)
  }))
  step_dt <- rep(all_dt, each = 1)
  speeds <- sqrt(rowSums(displacements^2)) / step_dt
  speeds <- speeds[is.finite(speeds) & speeds > 0]
  velocity_scale <- median(speeds, na.rm = TRUE) / sqrt(2)
  if (!is.finite(velocity_scale) || velocity_scale <= 0) {
    velocity_scale <- 1
  }
  obs_noise <- median(sqrt(rowSums(displacements^2)), na.rm = TRUE) / 10
  if (!is.finite(obs_noise) || obs_noise <= 0) {
    obs_noise <- 0.1
  }
  tau <- max(median(all_dt, na.rm = TRUE) * 2, min(all_dt, na.rm = TRUE))
  c(tau = tau, velocity_scale = velocity_scale, obs_noise = obs_noise)
}

positive_range <- function(x, fallback, probs = c(0.05, 0.95)) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) {
    return(sort(fallback))
  }
  out <- unname(stats::quantile(x, probs, na.rm = TRUE, names = FALSE))
  out <- sort(out)
  if (!all(is.finite(out)) || any(out <= 0)) {
    out <- sort(fallback)
  }
  if (out[1] == out[2]) {
    out <- out * c(0.5, 2)
  }
  out
}

with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

random_ou_starts <- function(segments, n_starts, seed = 17L) {
  empirical <- initial_ou_params(segments)
  n_starts <- as.integer(n_starts)
  if (is.na(n_starts) || n_starts < 1) {
    stop("n_starts must be a positive integer", call. = FALSE)
  }
  if (n_starts == 1L) {
    return(as.data.frame(as.list(empirical)))
  }

  all_dt <- unlist(lapply(segments, `[[`, "dt"), use.names = FALSE)
  all_dt <- all_dt[is.finite(all_dt) & all_dt > 0]
  displacements <- do.call(rbind, lapply(segments, function(seg) {
    diff(seg$y)
  }))
  step_displacements <- sqrt(rowSums(displacements^2))
  step_dt <- unlist(lapply(segments, `[[`, "dt"), use.names = FALSE)
  speeds <- step_displacements / step_dt

  dt_min <- min(all_dt, na.rm = TRUE)
  dt_med <- median(all_dt, na.rm = TRUE)
  if (!is.finite(dt_min) || dt_min <= 0) {
    dt_min <- empirical[["tau"]] / 2
  }
  if (!is.finite(dt_med) || dt_med <= 0) {
    dt_med <- empirical[["tau"]] / 2
  }
  tau_range <- sort(c(dt_min, max(dt_min * 2, dt_med * 50)))
  speed_range <- positive_range(speeds, fallback = empirical[["velocity_scale"]] * sqrt(2) * c(0.25, 4))
  disp_range <- positive_range(step_displacements, fallback = empirical[["obs_noise"]] * 10 * c(0.25, 4))

  random <- with_seed(seed, {
    data.frame(
      tau = exp(runif(n_starts - 1L, log(tau_range[1]), log(tau_range[2]))),
      velocity_scale = exp(runif(n_starts - 1L, log(speed_range[1]), log(speed_range[2]))) / sqrt(2),
      obs_noise = exp(runif(n_starts - 1L, log(disp_range[1] * 0.02), log(disp_range[2] * 0.5)))
    )
  })

  dplyr::bind_rows(as.data.frame(as.list(empirical)), random)
}

fit_one_ou_start <- function(segments, start, start_id) {
  lower <- log(c(1e-4, 1e-8, 1e-8))
  upper <- log(c(1e6, 1e6, 1e6))
  start <- pmin(pmax(log(unname(unlist(start[c("tau", "velocity_scale", "obs_noise")]))), lower), upper)
  fit <- tryCatch(
    optim(
      par = start,
      fn = function(par) -ou_loglik(par, segments),
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 1000)
    ),
    error = function(e) {
      list(
        par = rep(NA_real_, 3),
        value = Inf,
        convergence = NA_integer_,
        error_message = conditionMessage(e)
      )
    }
  )

  params <- exp(fit$par)
  tibble(
    start_id = start_id,
    start_tau = exp(start[1]),
    start_velocity_scale = exp(start[2]),
    start_obs_noise = exp(start[3]),
    tau = unname(params[1]),
    velocity_scale = unname(params[2]),
    effective_diffusivity = unname(params[2]^2 * params[1]),
    obs_noise = unname(params[3]),
    log_likelihood = -fit$value,
    fit_status = if (!is.null(fit$error_message)) {
      paste0("error: ", fit$error_message)
    } else if (fit$convergence == 0) {
      "converged"
    } else {
      paste0("optim_code_", fit$convergence)
    }
  )
}

fit_ou_velocity <- function(segments, n_starts = 25L, seed = 17L, start = NULL, start_id = NULL) {
  starts <- if (!is.null(start)) {
    as.data.frame(as.list(start))
  } else {
    random_ou_starts(segments, n_starts = n_starts, seed = seed)
  }
  starts$start_id <- seq_len(nrow(starts))

  if (!is.null(start_id)) {
    start_id <- as.integer(start_id)
    if (is.na(start_id) || start_id < 1 || start_id > nrow(starts)) {
      stop("start_id must be between 1 and n_starts", call. = FALSE)
    }
    starts <- starts[starts$start_id == start_id, , drop = FALSE]
  }

  fits <- dplyr::bind_rows(lapply(seq_len(nrow(starts)), function(i) {
    fit_one_ou_start(segments, starts[i, , drop = FALSE], start_id = starts$start_id[[i]])
  }))

  fits |>
    arrange(.data$start_id) |>
    mutate(
      n_starts = as.integer(n_starts),
      seed = seed
    ) |>
    select(
      "tau", "velocity_scale", "effective_diffusivity", "obs_noise",
      "log_likelihood", "fit_status", "n_starts", "seed", "start_id",
      "start_tau", "start_velocity_scale", "start_obs_noise"
    )
}
