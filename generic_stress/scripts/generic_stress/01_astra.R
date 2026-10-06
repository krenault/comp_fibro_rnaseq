#!/usr/bin/env Rscript
# Download / cache ASTRA Zenodo tables and write a comparison catalog + per-gene z.
# Drive is not touched here.

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
})

.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
source("scripts/shared/geo_helpers.R")

RAW <- file.path(ROOT, "raw/astra")
OUT <- file.path(ROOT, "results/astra")
dir.create(RAW, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ZENODO <- "https://zenodo.org/records/15885686/files/%s?download=1"
need <- c("datasets.tsv", "differential_expression.tsv")

copy_if_present <- function(name) {
  dest <- file.path(RAW, name)
  if (file.exists(dest) && file.info(dest)$size > 100) return(dest)
  tmp <- file.path("/tmp", paste0("astra_", sub("\\.tsv$", "", name), if (name == "differential_expression.tsv") "" else ""), name)
  # known /tmp names from earlier download
  alts <- c(
    file.path("/tmp", name),
    file.path("/tmp", paste0("astra_", name)),
    "/tmp/astra_de.tsv",
    "/tmp/astra_datasets.tsv"
  )
  if (name == "differential_expression.tsv") alts <- c("/tmp/astra_de.tsv", alts)
  if (name == "datasets.tsv") alts <- c("/tmp/astra_datasets.tsv", alts)
  hit <- alts[file.exists(alts)]
  if (length(hit)) {
    ok <- file.copy(hit[[1]], dest, overwrite = TRUE)
    if (ok) return(dest)
  }
  url <- sprintf(ZENODO, name)
  message("Downloading ", name)
  download.file(url, dest, mode = "wb", quiet = TRUE)
  dest
}

ds_path <- copy_if_present("datasets.tsv")
de_path <- copy_if_present("differential_expression.tsv")

cat("=== ASTRA ingest ===\n")
ds <- fread(ds_path)
setnames(ds, names(ds), gsub(" ", "_", names(ds)))
ds[, record_id := as.character(id)]
ds[, Cell_group := gsub("\\s+", " ", trimws(Cell_group))]
ds[, stress_class := fcase(
  stress_factor == "heat", "heat",
  stress_factor == "hypoxia", "hypoxia",
  stress_factor == "h2o2", "oxidative",
  stress_factor == "uv", "dna_damage",
  default = NA_character_
)]
ds[, is_nondiased_wt := Cell_group == "Reference Group (Non-diseased)" &
     Genotype == "WT" & Genomic_condition == "Wild-Type"]
ds[, experiment_id := paste(GEO_id, Cell_type, stress_class, Genotype, sep = "|")]

de <- fread(de_path)
de[, record_id := as.character(record_id)]
de[, gene := strip(gene_id)]
de <- de[grepl("^ENSG", gene) & is.finite(logfc)]
de[, z := as.numeric(scale(logfc)), by = record_id]

cat <- ds[record_id %in% unique(de$record_id)]
cat("Comparisons in DE: ", nrow(cat), "  genes x rows: ", nrow(de), "\n", sep = "")
print(cat[, .N, by = .(stress_class, Cell_group)])

fwrite(as.data.frame(cat), file.path(OUT, "comparison_catalog.csv"))
saveRDS(de[, .(record_id, gse_id, gene, logfc, pvalue, z)], file.path(OUT, "astra_z_long.rds"))
saveRDS(as.data.frame(cat), file.path(OUT, "comparison_catalog.rds"))
cat("Wrote results/astra/\n")
