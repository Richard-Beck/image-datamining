#!/usr/bin/env Rscript

usage <- paste0(
  "Usage: render_ou_tracking_fit_report.R [options]\n\n",
  "Options:\n",
  "  --repo_root=/path/to/repo\n",
  "  --input=/path/to/ou_tracking_fit_report.Rmd\n",
  "  --output_file=ou_tracking_fit_report.html\n"
)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE))
analysis_dir_guess <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(analysis_dir_guess, "R/k00_batch_utils.R"))

args <- parse_cli_args(commandArgs(trailingOnly = TRUE), usage)
repo_root <- normalizePath(args$repo_root %||% normalizePath(file.path(analysis_dir_guess, "../.."), mustWork = TRUE), mustWork = TRUE)
analysis_dir <- file.path(repo_root, "analyses/K00_GemcitabineExposure_033023")
input <- normalizePath(args$input %||% file.path(analysis_dir, "ou_tracking_fit_report.Rmd"), mustWork = TRUE)
output_file <- args$output_file %||% "ou_tracking_fit_report.html"

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The rmarkdown package is required to render the report.", call. = FALSE)
}

rmarkdown::render(
  input,
  output_file = output_file,
  output_dir = dirname(input),
  quiet = FALSE
)
