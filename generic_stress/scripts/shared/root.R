# Analysis root and Drive DGE path. Override with FLEX_ROOT / DRIVE_DGE.
# Scripts live in scripts/generic_stress/; this file lives in scripts/shared/.

flex_this_file <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit)) {
    return(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = FALSE))
  }
  NA_character_
}

flex_root <- function() {
  env <- Sys.getenv("FLEX_ROOT", unset = "")
  if (nzchar(env)) return(normalizePath(env, mustWork = TRUE))
  f <- flex_this_file()
  if (!is.na(f) && nzchar(f) && file.exists(f)) {
    return(normalizePath(file.path(dirname(f), "..", "..")))
  }
  d <- normalizePath(getwd())
  for (i in seq_len(8L)) {
    if (dir.exists(file.path(d, "scripts", "generic_stress"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Set FLEX_ROOT to the analysis root (the folder that contains scripts/generic_stress/).")
}

ROOT <- flex_root()
setwd(ROOT)

DRIVE_DGE <- Sys.getenv("DRIVE_DGE", unset = "")
if (!nzchar(DRIVE_DGE)) DRIVE_DGE <- file.path(ROOT, "drive_dge")
DRIVE_DGE <- normalizePath(DRIVE_DGE, mustWork = FALSE)

list_drive_dge <- function() {
  if (!dir.exists(DRIVE_DGE)) {
    stop(
      "Drive DGE folder not found: ", DRIVE_DGE,
      "\nPut pairwise DESeq2 CSVs there, or set DRIVE_DGE.",
      "\nPull: bash scripts/pull_drive_dge.sh"
    )
  }
  list.files(DRIVE_DGE, pattern = "\\.csv$", full.names = TRUE)
}
