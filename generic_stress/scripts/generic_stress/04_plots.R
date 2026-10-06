#!/usr/bin/env Rscript
# 04 — figures (descriptive names in plots/mammalian_generic_stress/).
.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
here <- "scripts/generic_stress"
run <- function(s) {
  st <- system2("Rscript", s)
  if (!identical(st, 0L)) stop(paste(s, "failed"))
}
run(file.path(here, "04_locked_plots.R"))
run(file.path(here, "04_volcano.R"))
