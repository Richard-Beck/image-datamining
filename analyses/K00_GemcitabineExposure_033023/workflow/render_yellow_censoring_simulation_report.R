#!/usr/bin/env Rscript

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
analysis_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

rmarkdown::render(
  file.path(analysis_dir, "yellow_censoring_simulation_report.Rmd"),
  output_dir = analysis_dir,
  clean = TRUE
)
