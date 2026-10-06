#!/usr/bin/env Rscript
# Leave-one-class-out generic CSR, then score held-out Drive.
#
# Drive libraries are not in discovery. Heat and hypoxia *classes* are:
# ASTRA+GEO discovery includes public heat and hypoxia, so Drive heat/hypoxia
# recovering the full generic list is a matched-class positive control, not
# evidence of a class-independent CSR. Cold has no discovery counterpart.
#
# Honest tests, from the already-fit class mixed models (no remeta):
#   no_heat     = HMP of hypoxia + oxidative + DNA-damage
#   no_hypoxia  = HMP of heat + oxidative + DNA-damage
#   unmatched   = HMP of oxidative + DNA-damage only
# Score Drive heat on no_heat, Drive hypoxia on no_hypoxia, Drive cold on
# the full generic list (and on unmatched as a shared cross-class core).

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(fgsea)
})

.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
source("scripts/shared/geo_helpers.R")

RES <- file.path(ROOT, "results/mammalian_generic_stress")
CLASSES <- c("heat", "hypoxia", "oxidative", "dna_damage")
SIGN_AGREE <- 0.80
HMP_PADJ <- 0.05

hmp_twosided <- function(p, sign_ec, L = length(p)) {
  ok <- is.finite(p) & is.finite(sign_ec) & p > 0
  p <- p[ok]; s <- sign(sign_ec[ok])
  if (length(p) < 2L) return(list(p = NA_real_, dir = 0, n = length(p), frac_agree = NA_real_))
  s[s == 0] <- 1
  one_pos <- ifelse(s > 0, pmin(pmax(p / 2, 1e-300), 1 - 1e-16),
                    pmin(pmax(1 - p / 2, 1e-300), 1 - 1e-16))
  one_neg <- ifelse(s < 0, pmin(pmax(p / 2, 1e-300), 1 - 1e-16),
                    pmin(pmax(1 - p / 2, 1e-300), 1 - 1e-16))
  hmp_one <- function(ps) {
    if (requireNamespace("harmonicmeanp", quietly = TRUE)) {
      as.numeric(harmonicmeanp::p.hmp(ps, w = rep(1 / length(ps), length(ps)),
                                      L = L, multilevel = FALSE))
    } else {
      h <- length(ps) / sum(1 / ps)
      mu <- log(L) + 1 - digamma(1)
      sig <- pi / 2
      z <- (1 / h - mu) / sig
      min(1, max(h, stats::pnorm(z, lower.tail = FALSE) * 1.5))
    }
  }
  pp <- hmp_one(one_pos)
  pn <- hmp_one(one_neg)
  if (!is.finite(pp)) pp <- 1
  if (!is.finite(pn)) pn <- 1
  if (pp <= pn) list(p = min(1, 2 * pp), dir = 1, n = length(p),
                     frac_agree = mean(s > 0))
  else list(p = min(1, 2 * pn), dir = -1, n = length(p),
            frac_agree = mean(s < 0))
}

hmp_from_classes <- function(wide, class_keep, min_n, label) {
  p_cols <- paste0("p_class_", class_keep)
  ec_cols <- paste0("ec_class_", class_keep)
  L <- length(class_keep)
  rows <- lapply(seq_len(nrow(wide)), function(i) {
    pv <- unlist(wide[i, ..p_cols], use.names = FALSE)
    ev <- unlist(wide[i, ..ec_cols], use.names = FALSE)
    h <- hmp_twosided(pv, ev, L = L)
    data.table(
      gene = wide$gene[i],
      hmp_p = h$p, hmp_dir = h$dir, n_class = h$n,
      frac_agree = h$frac_agree,
      ec_mean = mean(ev[is.finite(ev)])
    )
  })
  out <- rbindlist(rows)
  out[, hmp_padj := p.adjust(hmp_p, method = "BH")]
  out[, in_set := is.finite(hmp_padj) & hmp_padj < HMP_PADJ &
        n_class >= min_n & frac_agree >= SIGN_AGREE]
  out[, direction := fcase(
    in_set & hmp_dir > 0, "up",
    in_set & hmp_dir < 0, "down",
    default = "none"
  )]
  out[, set := label]
  out[]
}

sigs <- fread(file.path(RES, "cell_meta_signatures.csv"))
mam <- fread(file.path(RES, "gene_human_mouse_generic.csv"))
class_wide <- dcast(
  sigs[signature %in% paste0("class_", CLASSES), .(gene, signature, p, ec)],
  gene ~ signature, value.var = c("p", "ec")
)

hold <- rbindlist(list(
  hmp_from_classes(class_wide, CLASSES, 3L, "full"),
  hmp_from_classes(class_wide, setdiff(CLASSES, "heat"), 2L, "no_heat"),
  hmp_from_classes(class_wide, setdiff(CLASSES, "hypoxia"), 2L, "no_hypoxia"),
  hmp_from_classes(class_wide, c("oxidative", "dna_damage"), 2L, "unmatched")
))
fwrite(hold, file.path(RES, "cell_meta_loso_class_genes.csv"))

set_n <- hold[, .(
  n = .N,
  n_sig = sum(in_set),
  n_up = sum(direction == "up"),
  n_down = sum(direction == "down")
), by = set]
fwrite(set_n, file.path(RES, "cell_meta_loso_class_counts.csv"))
cat("Leave-one-class-out gene-set sizes:\n")
print(set_n)

full_up <- hold[set == "full" & direction == "up", gene]
full_dn <- hold[set == "full" & direction == "down", gene]
# sanity vs script-25 generic
cat("Full HMP overlap with script-25 generic:",
    length(intersect(full_up, mam$gene[mam$in_generic & mam$direction == "up"])),
    "/", length(full_up), " up\n")

pathways <- list(
  full_up = hold[set == "full" & direction == "up", gene],
  full_down = hold[set == "full" & direction == "down", gene],
  no_heat_up = hold[set == "no_heat" & direction == "up", gene],
  no_heat_down = hold[set == "no_heat" & direction == "down", gene],
  no_hypoxia_up = hold[set == "no_hypoxia" & direction == "up", gene],
  no_hypoxia_down = hold[set == "no_hypoxia" & direction == "down", gene],
  unmatched_up = hold[set == "unmatched" & direction == "up", gene],
  unmatched_down = hold[set == "unmatched" & direction == "down", gene]
)
pathways <- pathways[vapply(pathways, length, integer(1)) >= 8L]
cat("GSEA sets:", paste(sprintf("%s=%d", names(pathways), lengths(pathways)), collapse = ", "), "\n")

# Drive DGE (glucose skipped)
drive_files <- list_drive_dge()
drive_long <- rbindlist(lapply(drive_files, function(f) {
  bn <- tools::file_path_sans_ext(basename(f))
  tt <- str_extract(bn, "glucose|hypoxia|temperature")
  cmp <- str_extract(bn, "(2\\.5mM_vs_8mM|30mM_vs_8mM|6H_vs_0|24H_vs_0|32C_vs_37C|41C_vs_37C)$")
  if (is.na(tt) || is.na(cmp) || tt == "glucose") return(NULL)
  sp <- str_remove(bn, paste0("_", tt, "_", cmp, "$"))
  exposure <- fcase(
    tt == "hypoxia" & cmp == "6H_vs_0", "hypoxia_6H",
    tt == "hypoxia" & cmp == "24H_vs_0", "hypoxia_24H",
    tt == "temperature" & grepl("32C", cmp), "cold",
    tt == "temperature" & grepl("41C", cmp), "heat",
    default = NA_character_
  )
  d <- as.data.table(read_csv(f, show_col_types = FALSE))
  pcol <- if ("pvalue" %in% names(d)) "pvalue" else if ("padj" %in% names(d)) "padj" else NA
  data.table(
    gene = strip(d$gene),
    lfc = as.numeric(d$log2FoldChange),
    p = if (!is.na(pcol)) as.numeric(d[[pcol]]) else NA_real_,
    species = sp,
    exposure = exposure,
    driver = paste0("drive_", sp, "_", exposure)
  )[is.finite(lfc) & grepl("^ENSG", gene) & !is.na(exposure)]
}))
drive_long[, rankstat := -log10(pmin(pmax(fifelse(is.finite(p) & p > 0, p, 1), 1e-300), 1)) * sign(lfc)]
drive_long[p == 0, rankstat := 300 * sign(lfc)]
cat("Drive contrasts:", uniqueN(drive_long$driver), " species:", uniqueN(drive_long$species), "\n")

gsea_one <- function(sub) {
  st <- setNames(sub$rankstat, sub$gene)
  st <- st[is.finite(st)]
  st <- st[!duplicated(names(st))]
  if (length(st) < 200) return(NULL)
  set.seed(1)
  fg <- tryCatch(
    fgsea::fgseaMultilevel(pathways = pathways, stats = st, minSize = 8, maxSize = 5000),
    error = function(e) NULL
  )
  if (is.null(fg) || !nrow(fg)) return(NULL)
  as.data.table(fg)[, `:=`(species = sub$species[1], exposure = sub$exposure[1],
                           driver = sub$driver[1])]
}

gsea_drv <- rbindlist(lapply(split(drive_long, by = "driver"), gsea_one), fill = TRUE)
gsea_drv[, padj := p.adjust(pval, method = "BH")]
fwrite(gsea_drv[, .(driver, species, exposure, pathway, pval, padj, NES, size)],
       file.path(RES, "drive_gsea_class_holdout.csv"))

# Honest query: each Drive exposure vs the gene set that does not contain its class
honest_map <- data.table(
  exposure = c("heat", "hypoxia_6H", "hypoxia_24H", "cold"),
  honest_up = c("no_heat_up", "no_hypoxia_up", "no_hypoxia_up", "full_up"),
  honest_down = c("no_heat_down", "no_hypoxia_down", "no_hypoxia_down", "full_down"),
  note = c(
    "Drive heat vs generic built without public heat",
    "Drive hypoxia vs generic built without public hypoxia",
    "Drive hypoxia vs generic built without public hypoxia",
    "Drive cold vs full generic (cold is not a discovery class)"
  )
)
honest <- rbindlist(lapply(seq_len(nrow(honest_map)), function(i) {
  m <- honest_map[i]
  u <- gsea_drv[exposure == m$exposure & pathway == m$honest_up]
  d <- gsea_drv[exposure == m$exposure & pathway == m$honest_down]
  merge(
    u[, .(species, exposure, NES_up = NES, padj_up = padj, set_up = pathway, size_up = size)],
    d[, .(species, exposure, NES_down = NES, padj_down = padj, set_down = pathway, size_down = size)],
    by = c("species", "exposure"), all = TRUE
  )[, note := m$note]
}))
fwrite(honest, file.path(RES, "drive_gsea_honest_holdout.csv"))

# Per-species Spearman: Drive LFC vs public generic EC (full HMP ranking)
ec <- mam[, .(gene, ec_generic = ec_mean, in_generic)]
rho_sp <- drive_long[, {
  x <- merge(.SD, ec, by = "gene")
  x <- x[is.finite(lfc) & is.finite(ec_generic)]
  all <- if (nrow(x) >= 30) {
    ct <- suppressWarnings(cor.test(x$lfc, x$ec_generic, method = "spearman", exact = FALSE))
    list(n_all = nrow(x), rho_all = unname(ct$estimate), p_all = ct$p.value)
  } else list(n_all = nrow(x), rho_all = NA_real_, p_all = NA_real_)
  g <- x[in_generic == TRUE]
  gen <- if (nrow(g) >= 20) {
    ct <- suppressWarnings(cor.test(g$lfc, g$ec_generic, method = "spearman", exact = FALSE))
    list(n_generic = nrow(g), rho_generic = unname(ct$estimate), p_generic = ct$p.value)
  } else list(n_generic = nrow(g), rho_generic = NA_real_, p_generic = NA_real_)
  c(all, gen)
}, by = .(species, exposure)]
fwrite(rho_sp, file.path(RES, "drive_vs_generic_spearman_by_species.csv"))

# Human vs other mammals
cmp <- merge(
  gsea_drv[pathway == "full_up", .(species, exposure, NES_full_up = NES, padj_full_up = padj)],
  honest[, .(species, exposure, NES_honest_up = NES_up, padj_honest_up = padj_up, note)],
  by = c("species", "exposure"), all = TRUE
)
cmp[, group := fifelse(species == "human", "Human Drive", "Other Drive mammals")]
fwrite(cmp, file.path(RES, "drive_human_vs_other_mammals.csv"))

cat("\nMedian full-generic-up NES by exposure × group:\n")
print(cmp[, .(med = median(NES_full_up, na.rm = TRUE), n = .N),
          by = .(exposure, group)][order(exposure, group)])
cat("\nMedian honest-up NES by exposure × group:\n")
print(cmp[, .(med = median(NES_honest_up, na.rm = TRUE), n = .N),
          by = .(exposure, group)][order(exposure, group)])
cat("\nMedian Spearman (generic genes) by exposure × human vs other:\n")
print(rho_sp[, .(
  med = median(rho_generic, na.rm = TRUE),
  n = .N
), by = .(exposure, group = fifelse(species == "human", "human", "other"))][order(exposure, group)])

writeLines(c(
  "Leave-one-class-out generic CSR vs Drive",
  "",
  "Discovery and Drive share heat and hypoxia as classes, not as samples.",
  "Drive recovering the full generic list after heat/hypoxia is a matched-class",
  "positive control. The claim of a class-independent CSR uses:",
  "  - Drive heat vs no_heat generic (public heat dropped from HMP)",
  "  - Drive hypoxia vs no_hypoxia generic",
  "  - Drive cold vs full generic (cold is absent from discovery)",
  "  - optional unmatched core = oxidative + DNA-damage only",
  "",
  "Cell 2023 does not leave-one-dataset-out of the 92-dataset aging meta.",
  "They fit tissue- and species-stratified signatures, and use species LOO only",
  "for Elastic Net lifespan prediction. Class hold-out is the analog here."
), file.path(RES, "README_class_holdout.txt"))

cat("Done.\n")
