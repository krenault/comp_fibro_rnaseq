#!/usr/bin/env Rscript
# DESeq2 per species × treatment axis.
#
# Usage:
#   Rscript run_deseq2.R --counts FILE --outdir DIR [--mode pairwise|multi]
#
# pairwise (default, production): one fit per contrast; filter within contrast samples.
# multi: one fit per axis with all levels together.
suppressPackageStartupMessages({ library(tidyverse); library(DESeq2) })

args <- commandArgs(trailingOnly = TRUE)
opt <- list(counts = NULL, outdir = "DESeq_results", mode = "pairwise")
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[[i]])
  if (!key %in% names(opt) || i == length(args)) stop("bad arg: ", args[[i]])
  opt[[key]] <- args[[i + 1L]]
  i <- i + 2L
}
if (is.null(opt$counts)) stop("required: --counts FILE")
opt$mode <- tolower(opt$mode)
if (opt$mode %in% c("multilevel", "multi-level", "pooled")) opt$mode <- "multi"
if (opt$mode == "singleside") opt$mode <- "pairwise"
if (!opt$mode %in% c("multi", "pairwise")) stop("--mode must be pairwise or multi")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

parse_sample_name <- function(sample_name) {
  parts <- strsplit(sample_name, "_")[[1]]
  n <- length(parts)
  list(
    species    = paste(parts[1:(n - 4)], collapse = "_"),
    individual = parts[n - 3],
    temp       = parts[n - 2],
    glucose    = parts[n - 1],
    hypoxia    = parts[n]
  )
}

treatment_axes   <- c("glucose", "hypoxia", "temperature")
treatment_levels <- list(
  glucose = c("8mM", "2.5mM", "30mM"),
  hypoxia = c("0", "6H", "24H"),
  temperature = c("37C", "32C", "41C")
)
treatment_column <- list(glucose = "glucose", hypoxia = "hypoxia", temperature = "temp")
contrasts_vs_control <- list(
  glucose = list(c("2.5mM", "8mM"), c("30mM", "8mM")),
  hypoxia = list(c("6H", "0"), c("24H", "0")),
  temperature = list(c("32C", "37C"), c("41C", "37C"))
)

filter_genes <- function(counts, min_count = 5, min_samples = 3) {
  counts[rowSums(counts >= min_count) >= min_samples, , drop = FALSE]
}

save_contrast <- function(res_df, species_i, axis_j, comparison, n_samples, n_ind, summary_rows) {
  write_csv(res_df, file.path(opt$outdir, paste0(species_i, "_", axis_j, "_", comparison, ".csv")))
  n_degs <- sum(res_df$padj < 0.05, na.rm = TRUE)
  message("  ", comparison, ": ", n_degs, " DEGs (padj < 0.05)")
  summary_rows[[comparison]] <- tibble(
    species = species_i, treatment_type = axis_j, comparison = comparison,
    n_samples = n_samples, n_individuals = n_ind, n_genes_tested = nrow(res_df),
    n_DEGs_padj05 = n_degs, design_formula = "~ individual + treatment",
    mode = opt$mode
  )
  summary_rows
}

fit_deseq <- function(counts, meta) {
  counts[is.na(counts)] <- 0
  counts <- filter_genes(counts)
  if (nrow(counts) < 10) return(NULL)
  DESeq(DESeqDataSetFromMatrix(round(counts), meta, ~ individual + treatment), quiet = TRUE)
}

extract_res <- function(dds, trt, ctl) {
  as.data.frame(results(dds, contrast = c("treatment", trt, ctl))) %>%
    rownames_to_column("gene") %>%
    select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    arrange(padj)
}

axis_meta <- function(species_meta, axis_j) {
  if (axis_j == "glucose") species_meta %>% filter(temp == "37C", hypoxia == "0")
  else if (axis_j == "hypoxia") species_meta %>% filter(temp == "37C", glucose == "8mM")
  else species_meta %>% filter(glucose == "8mM", hypoxia == "0")
}

cat("mode:", opt$mode, " ->", opt$outdir, "\n")
all_counts <- read.csv(opt$counts, row.names = 1, check.names = FALSE)
sample_meta <- map_dfr(colnames(all_counts), function(s) {
  p <- parse_sample_name(s)
  tibble(sample = s, species = p$species, individual = p$individual,
         temp = p$temp, glucose = p$glucose, hypoxia = p$hypoxia)
})
species_list <- unique(sample_meta$species)
cat("loaded", ncol(all_counts), "samples,", nrow(all_counts), "genes;",
    length(species_list), "species\n")

all_summaries <- list()
for (species_i in species_list) {
  cat("\n===", toupper(species_i), "===\n")
  for (axis_j in treatment_axes) {
    cat("processing", axis_j, "...\n")
    tryCatch({
      experiment_meta <- axis_meta(filter(sample_meta, species == species_i), axis_j)
      treat_col <- treatment_column[[axis_j]]
      if (nrow(experiment_meta) < 4 || length(unique(experiment_meta[[treat_col]])) < 2) {
        message("  skip ", axis_j)
        next
      }
      experiment_meta$treatment <- factor(
        experiment_meta[[treat_col]],
        levels = intersect(treatment_levels[[axis_j]], unique(experiment_meta[[treat_col]]))
      )
      experiment_meta$individual <- factor(experiment_meta$individual)
      summary_rows <- list()

      if (opt$mode == "multi") {
        dds <- fit_deseq(all_counts[, experiment_meta$sample, drop = FALSE], experiment_meta)
        if (is.null(dds)) next
        for (ct in contrasts_vs_control[[axis_j]]) {
          if (!all(ct %in% levels(experiment_meta$treatment))) next
          summary_rows <- save_contrast(
            extract_res(dds, ct[[1]], ct[[2]]), species_i, axis_j,
            paste0(ct[[1]], "_vs_", ct[[2]]),
            nrow(experiment_meta), n_distinct(experiment_meta$individual), summary_rows
          )
        }
      } else {
        for (ct in contrasts_vs_control[[axis_j]]) {
          trt <- ct[[1]]; ctl <- ct[[2]]
          if (!all(c(trt, ctl) %in% levels(experiment_meta$treatment))) next
          cm <- experiment_meta %>%
            filter(treatment %in% c(trt, ctl)) %>%
            mutate(treatment = factor(as.character(treatment), levels = c(ctl, trt)))
          if (nrow(cm) < 4) next
          dds <- fit_deseq(all_counts[, cm$sample, drop = FALSE], cm)
          if (is.null(dds)) next
          summary_rows <- save_contrast(
            extract_res(dds, trt, ctl), species_i, axis_j, paste0(trt, "_vs_", ctl),
            nrow(cm), n_distinct(cm$individual), summary_rows
          )
        }
      }
      if (length(summary_rows)) {
        all_summaries[[paste0(species_i, "_", axis_j)]] <- bind_rows(summary_rows)
      }
    }, error = function(e) cat("  error:", conditionMessage(e), "\n"))
  }
}

if (length(all_summaries)) {
  final_summary <- bind_rows(all_summaries)
  write_csv(final_summary, file.path(opt$outdir, "DEG_summary_all_species_treatments.csv"))
  print(final_summary, n = Inf)
} else {
  cat("no analyses completed.\n")
}
