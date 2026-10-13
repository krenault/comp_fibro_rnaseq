#!/usr/bin/env Rscript
# Sensitivity: primary discovery vs the same discovery without GSE328685.
# GSE328685 is wild-type tendon fibroblast hypoxia for 4 days, not a 6–24 h pulse.
# The primary signature keeps it. This script asks whether membership, sign, or
# FH experiment hold-out GSEA calls move when it is removed.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(scales)
})

.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
source("scripts/shared/plot_aesthetics.R")
source("scripts/generic_stress/plot_helpers.R")

PRI <- file.path(ROOT, "results/mammalian_generic_stress")
SEN <- Sys.getenv("GENERIC_RES", unset = "")
if (!nzchar(SEN)) SEN <- file.path(PRI, "sens_no_GSE328685")
if (!dir.exists(SEN)) stop("Sensitivity folder missing: ", SEN)

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (!length(a) && !length(b)) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

call_of <- function(nes, padj) {
  fcase(
    is.finite(nes) & is.finite(padj) & nes > 0 & padj < 0.05, "Induced",
    is.finite(nes) & is.finite(padj) & nes < 0 & padj < 0.05, "Opposite",
    default = "Unclear"
  )
}

read_gen <- function(dir, file) {
  p <- file.path(dir, file)
  if (!file.exists(p)) return(data.table())
  x <- fread(p)
  if (!"in_generic" %in% names(x) && "in_set" %in% names(x)) x[, in_generic := in_set]
  if (!"direction" %in% names(x)) x[, direction := NA_character_]
  x[in_generic %in% TRUE, .(gene, direction)]
}

cmp_lists <- function(a, b, label) {
  if (!nrow(a) && !nrow(b)) {
    return(data.table(
      comparison = label, n_a = 0L, n_b = 0L, n_shared = 0L,
      jaccard = 1, n_same_dir = 0L, n_sign_flips = 0L
    ))
  }
  shared <- merge(a, b, by = "gene", suffixes = c("_a", "_b"))
  data.table(
    comparison = label,
    n_a = uniqueN(a$gene), n_b = uniqueN(b$gene),
    n_shared = nrow(shared),
    jaccard = jaccard(a$gene, b$gene),
    n_same_dir = sum(shared$direction_a == shared$direction_b, na.rm = TRUE),
    n_sign_flips = sum(shared$direction_a != shared$direction_b, na.rm = TRUE)
  )
}

pri_pan <- read_gen(PRI, "gene_human_mouse_generic.csv")
sen_pan <- read_gen(SEN, "gene_human_mouse_generic.csv")
pri_fib <- read_gen(PRI, "gene_fibroblast_generic.csv")
sen_fib <- read_gen(SEN, "gene_fibroblast_generic.csv")
old_1344 <- read_gen(PRI, "mammalian_generic_genes_pre_fibro_hypoxia.csv")
if (!nrow(old_1344) && file.exists(file.path(PRI, "mammalian_generic_genes_pre_fibro_hypoxia.csv"))) {
  old_1344 <- fread(file.path(PRI, "mammalian_generic_genes_pre_fibro_hypoxia.csv"))
  if (!"direction" %in% names(old_1344)) old_1344[, direction := NA_character_]
  old_1344 <- old_1344[, .(gene, direction)]
}

jac <- rbindlist(list(
  cmp_lists(pri_pan, sen_pan, "primary_vs_sens_no_GSE328685"),
  cmp_lists(pri_fib, sen_fib, "fibroblast_primary_vs_sens_no_GSE328685"),
  cmp_lists(pri_pan, pri_fib, "pancell_vs_fibroblast"),
  cmp_lists(pri_pan, old_1344, "primary_vs_pre_fibro_1344")
), fill = TRUE)
fwrite(jac, file.path(PRI, "generic_qc_drop_jaccard.csv"))
fwrite(jac, file.path(SEN, "generic_qc_drop_jaccard.csv"))
cat("Gene-list Jaccard:\n")
print(jac)

read_honest <- function(dir, file) {
  p <- file.path(dir, file)
  if (!file.exists(p)) return(data.table())
  x <- fread(p)
  x <- x[!exposure %in% c("glucose_low", "glucose_high") & species != "rousette"]
  x[, call := call_of(NES_up, padj_up)]
  x[]
}

gsea_cmp <- function(a, b, tag) {
  if (!nrow(a) || !nrow(b)) return(data.table())
  m <- merge(
    a[, .(species, exposure, NES_pri = NES_up, padj_pri = padj_up, call_pri = call)],
    b[, .(species, exposure, NES_sens = NES_up, padj_sens = padj_up, call_sens = call)],
    by = c("species", "exposure")
  )
  m[, tag := tag]
  m[, call_same := call_pri == call_sens]
  m[]
}

gsea <- rbindlist(list(
  gsea_cmp(
    read_honest(PRI, "drive_gsea_honest_holdout.csv"),
    read_honest(SEN, "drive_gsea_honest_holdout.csv"),
    "pan-cell"
  ),
  gsea_cmp(
    read_honest(PRI, "drive_gsea_honest_holdout_fibroblast.csv"),
    read_honest(SEN, "drive_gsea_honest_holdout_fibroblast.csv"),
    "fibroblast"
  )
), fill = TRUE)
if (nrow(gsea)) {
  fwrite(gsea, file.path(PRI, "sens_no_GSE328685_drive_gsea_calls.csv"))
  fwrite(gsea, file.path(SEN, "sens_no_GSE328685_drive_gsea_calls.csv"))
  cat("\nFH experiments honest GSEA call concordance:\n")
  print(gsea[, .(
    n = .N,
    n_same_call = sum(call_same),
    n_flip = sum(!call_same),
    n_induced_pri = sum(call_pri == "Induced"),
    n_induced_sens = sum(call_sens == "Induced")
  ), by = tag])
}

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
jac_plot <- jac[comparison %in% c(
  "primary_vs_sens_no_GSE328685",
  "fibroblast_primary_vs_sens_no_GSE328685"
)]
jac_plot[, lab := fcase(
  comparison == "primary_vs_sens_no_GSE328685", "Pan-cell generic\n(drop 4-day tendon hypoxia)",
  comparison == "fibroblast_primary_vs_sens_no_GSE328685", "Fibroblast-only generic\n(drop 4-day tendon hypoxia)",
  default = comparison
)]
if (nrow(jac_plot)) {
  p_j <- ggplot(jac_plot, aes(jaccard, lab)) +
    geom_col(width = 0.62, fill = unname(treatment_colors["stress"]),
             color = "grey25", linewidth = 0.35) +
    geom_text(
      aes(label = sprintf("Jaccard %.2f   %d of %d shared   %d sign flips",
                          jaccard, n_shared, n_a, n_sign_flips)),
      hjust = -0.02, fontface = "bold", size = 3.2
    ) +
    coord_cartesian(xlim = c(0, 1.55)) +
    labs(
      title = "Dropping GSE328685 does not rewrite the generic list",
      subtitle = "Leave-one-series: 4-day tendon fibroblast hypoxia vs the rest of discovery. Sign flips counted on shared genes.",
      x = "Jaccard overlap with the primary list", y = NULL
    ) +
    theme_pub(12)
  save_png(p_j, "sens_drop_GSE328685_jaccard.png", 9.2, 4.6)
}

if (nrow(gsea)) {
  species_levels <- c(
    "human", "gelada", "rat", "squirrel", "little_brown_bat",
    "honey_badger", "seal", "dolphin", "whale", "rhino",
    "dromedary_camel", "bactrian_camel"
  )
  call_cols <- c(Induced = "#7CB577", Unclear = "#D0D4D8", Opposite = "#2166AC")
  gd <- gsea[tag == "pan-cell"]
  gd[, species_lab := factor(pretty_species(species), levels = pretty_species(species_levels))]
  gd[, exp_lab := factor(pretty_exp(exposure),
                         levels = c("Hypoxia 6H", "Hypoxia 24H", "Cold 32C", "Heat 41C"))]
  n_flip <- sum(!gd$call_same)
  p_c <- ggplot(gd, aes(exp_lab, species_lab, fill = call_pri)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = ifelse(call_same, "", "x")),
              fontface = "bold", size = 3, color = "grey15") +
    scale_fill_manual(values = call_cols, name = "Primary call", drop = FALSE) +
    labs(
      title = "FH experiment hold-out GSEA calls are unchanged without GSE328685",
      subtitle = sprintf(
        "Fill = primary pan-cell honest GSEA. x = call differs after dropping GSE328685 (%d of %d contrasts).",
        n_flip, nrow(gd)
      ),
      x = NULL, y = NULL
    ) +
    theme_pub(12) +
    theme(axis.text.x = element_text(face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid = element_blank(),
          legend.position = "bottom")
  save_png(p_c, "sens_drop_GSE328685_drive_calls.png", 8.6, 6.8)
}

cat("Sensitivity comparison written to ", PRI, " and ", SEN, "\n", sep = "")
