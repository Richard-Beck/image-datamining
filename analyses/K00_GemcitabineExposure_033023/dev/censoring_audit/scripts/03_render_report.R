#!/usr/bin/env Rscript

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)
audit_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

rmarkdown::render(
  file.path(audit_dir, "ou_censoring_simulation_report.Rmd"),
  output_dir = audit_dir,
  clean = TRUE
)
