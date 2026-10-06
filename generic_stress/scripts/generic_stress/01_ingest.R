#!/usr/bin/env Rscript
# 01 — discovery ingest. ASTRA Zenodo + extra cell GEO. Drive held out.
.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
here <- "scripts/generic_stress"
run <- function(cmd, args) {
  st <- system2(cmd, args)
  if (!identical(st, 0L)) stop(paste(cmd, paste(args, collapse = " "), "failed"))
}
run("Rscript", file.path(here, "01_astra.R"))
run("Rscript", file.path(here, "01_geo.R"))
run("python3", file.path(here, "01_metadata.py"))
