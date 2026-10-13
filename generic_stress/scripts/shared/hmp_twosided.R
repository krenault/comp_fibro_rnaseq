# Direction-aware two-sided HMP. The function body is share/hmp_combine.R,
# which is the file to send on its own. This wrapper keeps the old source path.

if (!exists("ROOT", inherits = TRUE)) {
  stop("source scripts/shared/root.R before scripts/shared/hmp_twosided.R")
}
source(file.path(get("ROOT", inherits = TRUE), "share", "hmp_combine.R"), local = FALSE)
