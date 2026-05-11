#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rmarkdown)
})

script_file <- normalizePath(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]),
  mustWork = TRUE
)
report_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
input <- file.path(report_dir, "alternative_motility_model_report.Rmd")

rmarkdown::render(
  input = input,
  output_dir = report_dir,
  clean = TRUE,
  quiet = FALSE
)
