#!/usr/bin/env Rscript
# REML stage-2 class models + cell-type sensitivities.
# Uses already Deming-scaled ECs from 02_public_hmp.R (no Deming rerun).
# Primary here: drop C2C12, REML tau^2. Sensitivities: no cancer; fibroblast-only.
# Does not overwrite the 02 generic list; writes parallel tables.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(metafor)
})

.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
source("scripts/shared/plot_aesthetics.R")

RES <- file.path(ROOT, "results/mammalian_generic_stress")
source("scripts/generic_stress/plot_helpers.R")

CLASSES <- c("heat", "hypoxia", "oxidative", "dna_damage")
SIGN_AGREE <- 0.80
HMP_PADJ <- 0.05
NCORES <- as.integer(Sys.getenv("NCORES", unset = max(1L, parallel::detectCores() - 1L)))

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
      z <- (1 / h - mu) / (pi / 2)
      min(1, max(h, stats::pnorm(z, lower.tail = FALSE) * 1.5))
    }
  }
  pp <- hmp_one(one_pos); pn <- hmp_one(one_neg)
  if (!is.finite(pp)) pp <- 1
  if (!is.finite(pn)) pn <- 1
  if (pp <= pn) list(p = min(1, 2 * pp), dir = 1, n = length(p), frac_agree = mean(s > 0))
  else list(p = min(1, 2 * pn), dir = -1, n = length(p), frac_agree = mean(s < 0))
}

fit_rma <- function(d) {
  d <- d[is.finite(yi) & is.finite(vi) & vi > 0]
  if (nrow(d) < 3L) return(NULL)
  d[, vi := pmin(pmax(vi, 1e-12), 1e4)]
  g <- d[, {
    w <- 1 / vi
    sw <- sum(w)
    .(yi = sum(w * yi) / sw, vi = 1 / sw)
  }, by = gse_id]
  if (nrow(g) < 3L) g <- d[, .(yi, vi)]
  m <- tryCatch(
    metafor::rma(yi, vi, data = g, method = "REML"),
    error = function(e) tryCatch(
      metafor::rma(yi, vi, data = g, method = "DL"),
      error = function(e2) tryCatch(
        metafor::rma(yi, vi, data = g, method = "EE"),
        error = function(e3) NULL
      )
    )
  )
  if (is.null(m)) return(NULL)
  c(ec = as.numeric(m$beta), se = as.numeric(m$se),
    p = as.numeric(m$pval), n = nrow(d))
}

fit_rma_mv <- function(d) {
  d <- d[is.finite(yi) & is.finite(vi) & vi > 0]
  if (nrow(d) < 3L) return(NULL)
  d[, vi := pmin(pmax(vi, 1e-12), 1e4)]
  d[, `:=`(
    gse_f = factor(gse_id),
    tissue_f = factor(tissue),
    species_f = factor(species)
  )]
  n_g <- uniqueN(d$gse_f)
  n_t <- uniqueN(d$tissue_f)
  n_s <- uniqueN(d$species_f)
  rand <- list(~ 1 | gse_f)
  if (n_t >= 3L) rand <- c(rand, list(~ 1 | tissue_f))
  if (n_s >= 2L) rand <- c(rand, list(~ 1 | species_f))
  m <- tryCatch(
    metafor::rma.mv(yi, vi, random = rand, data = d, method = "REML",
                    verbose = FALSE),
    error = function(e) tryCatch(
      metafor::rma.mv(yi, vi, random = ~ 1 | gse_f, data = d, method = "REML",
                      verbose = FALSE),
      error = function(e2) tryCatch(
        metafor::rma(yi, vi, data = d, method = "REML"),
        error = function(e3) tryCatch(
          metafor::rma(yi, vi, data = d, method = "DL"),
          error = function(e4) NULL
        )
      )
    )
  )
  if (is.null(m)) return(NULL)
  c(ec = as.numeric(m$beta)[1], se = as.numeric(m$se)[1],
    p = as.numeric(m$pval)[1], n = nrow(d))
}

meta_by_gene <- function(dt, min_n = 3L, label = "", fit = fit_rma) {
  cat("  mixed-effect", label, " contrasts=", uniqueN(dt$dataset_id),
      " GSE=", uniqueN(dt$gse_id), " genes=", uniqueN(dt$gene), "\n")
  spl <- split(dt, by = "gene", drop = TRUE)
  keep <- vapply(spl, nrow, integer(1)) >= min_n
  spl <- spl[keep]
  cat("    genes with n>=", min_n, ": ", length(spl), "  cores=", NCORES, "\n", sep = "")
  res <- parallel::mclapply(spl, fit, mc.cores = NCORES)
  ok <- !vapply(res, is.null, logical(1))
  if (!any(ok)) return(data.table())
  mat <- do.call(rbind, res[ok])
  out <- data.table(
    gene = names(spl)[ok],
    ec = as.numeric(mat[, "ec"]),
    se = as.numeric(mat[, "se"]),
    p = as.numeric(mat[, "p"]),
    n = as.integer(mat[, "n"])
  )
  out <- out[is.finite(p) & is.finite(ec)]
  out[, padj := p.adjust(p, method = "BH")]
  out[, significant := padj < 0.05]
  out[, signature := label]
  out[]
}

hmp_generic <- function(sigs, label) {
  wide <- dcast(
    sigs[signature %in% paste0("class_", CLASSES), .(gene, signature, p, ec)],
    gene ~ signature, value.var = c("p", "ec")
  )
  p_cols <- paste0("p_class_", CLASSES)
  ec_cols <- paste0("ec_class_", CLASSES)
  rows <- lapply(seq_len(nrow(wide)), function(i) {
    pv <- unlist(wide[i, ..p_cols], use.names = FALSE)
    ev <- unlist(wide[i, ..ec_cols], use.names = FALSE)
    h <- hmp_twosided(pv, ev, L = 4L)
    data.table(
      gene = wide$gene[i], hmp_p = h$p, hmp_dir = h$dir,
      n_class = h$n, frac_agree = h$frac_agree,
      ec_mean = mean(ev[is.finite(ev)])
    )
  })
  out <- rbindlist(rows)
  out[, hmp_padj := p.adjust(hmp_p, method = "BH")]
  out[, in_generic := is.finite(hmp_padj) & hmp_padj < HMP_PADJ &
        n_class >= 3L & frac_agree >= SIGN_AGREE]
  out[, direction := fcase(
    in_generic & hmp_dir > 0, "up",
    in_generic & hmp_dir < 0, "down",
    default = "none"
  )]
  out[, panel := label]
  out[]
}

run_panel <- function(dt, label, fit = fit_rma) {
  cat("\n=== Panel:", label, "===\n")
  sigs <- rbindlist(lapply(CLASSES, function(cl) {
    meta_by_gene(dt[stress_class == cl], min_n = 3L,
                 label = paste0("class_", cl), fit = fit)
  }), fill = TRUE)
  sigs[, panel := label]
  gen <- hmp_generic(sigs, label)
  cat("  generic:", sum(gen$in_generic),
      " up", sum(gen$direction == "up"),
      " down", sum(gen$direction == "down"), "\n")
  list(sigs = sigs, gen = gen)
}

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (!length(a) && !length(b)) return(1)
  length(intersect(a, b)) / length(union(a, b))
}

# -----------------------------------------------------------------------------
cat("Loading scaled EC matrix\n")
long <- as.data.table(readRDS(file.path(RES, "cell_meta_ec_long.rds")))
contr <- fread(file.path(RES, "cell_meta_contrasts.csv"))
old <- fread(file.path(RES, "gene_human_mouse_generic.csv"))
old_gen <- old[in_generic == TRUE, gene]
old_up <- old[in_generic == TRUE & direction == "up", gene]
old_dn <- old[in_generic == TRUE & direction == "down", gene]

fib_pat <- "fibroblast|WI-38|MRC5|IMR90|HLF|HPF|HDF|NIH 3T3|NIH3T3|MEF"
muscle_pat <- "C2C12|skeletal muscle|myoblast"
cancer_pat <- "HeLa|MCF7|HCT116|MDA-MB|LNCaP|U2OS|TK6|RAMOS|RCC4|Bewo|Pterygium"
contr[, cell_group := fcase(
  grepl("fibroblast-derived neuron", cell_type, ignore.case = TRUE), "other",
  grepl(muscle_pat, paste(cell_type, tissue), ignore.case = TRUE), "myoblast",
  grepl(fib_pat, paste(cell_type, tissue), ignore.case = TRUE), "fibroblast",
  grepl(cancer_pat, cell_type, ignore.case = TRUE), "cancer",
  default = "other"
)]
long <- merge(long, contr[, .(dataset_id, cell_group)], by = "dataset_id", all.x = TRUE)
long[is.na(cell_group), cell_group := "other"]
cat("Contrasts by group:\n")
print(contr[, .N, by = cell_group])

# Fast concordance: REML on the existing 1344 generic genes only (all contrasts)
cat("\n=== Fast REML concordance on the DL generic 1344 ===\n")
old_long <- long[gene %in% old_gen]
conc_sigs <- rbindlist(lapply(CLASSES, function(cl) {
  meta_by_gene(old_long[stress_class == cl], min_n = 3L, label = paste0("class_", cl))
}), fill = TRUE)
conc <- merge(
  old[in_generic == TRUE, .(gene, direction, ec_mean, hmp_padj)],
  conc_sigs[, .(gene, signature, ec_reml = ec, p_reml = p, padj_reml = padj)],
  by = "gene", all.x = TRUE, allow.cartesian = TRUE
)
sign_tab <- conc_sigs[, {
  w <- dcast(.SD, gene ~ signature, value.var = "ec")
  ev <- as.matrix(w[, paste0("class_", CLASSES), with = FALSE])
  data.table(
    gene = w$gene,
    reml_mean = rowMeans(ev, na.rm = TRUE),
    reml_n = rowSums(is.finite(ev))
  )
}]
sign_tab <- merge(sign_tab, old[in_generic == TRUE, .(gene, direction, ec_mean)], by = "gene")
sign_tab[, same_sign := sign(reml_mean) == sign(ec_mean) & reml_mean != 0 & ec_mean != 0]
cat("  genes remeta'd:", nrow(sign_tab),
    " same sign vs DL EC:", sum(sign_tab$same_sign, na.rm = TRUE),
    sprintf(" (%.1f%%)\n", 100 * mean(sign_tab$same_sign, na.rm = TRUE)))
fwrite(sign_tab, file.path(RES, "reml_concordance_generic1344.csv"))

# Full panels
panels <- list()
panels$reml_noc2c12 <- run_panel(long[cell_group != "myoblast"], "reml_noc2c12")
panels$reml_nocancer <- run_panel(long[cell_group != "cancer" & cell_group != "myoblast"], "reml_nocancer")
panels$reml_fibroblast <- run_panel(long[cell_group == "fibroblast"], "reml_fibroblast")
panels$rmamv_noc2c12 <- run_panel(
  long[cell_group != "myoblast"], "rmamv_noc2c12", fit = fit_rma_mv
)

all_sigs <- rbindlist(lapply(panels, `[[`, "sigs"), fill = TRUE)
all_gen <- rbindlist(lapply(panels, `[[`, "gen"), fill = TRUE)
fwrite(all_sigs, file.path(RES, "cell_meta_signatures_reml_celltype.csv"))
fwrite(all_gen, file.path(RES, "cell_meta_loso_celltype_genes.csv"))

sets <- list(
  dl_full = old_gen,
  reml_noc2c12 = all_gen[panel == "reml_noc2c12" & in_generic == TRUE, gene],
  reml_nocancer = all_gen[panel == "reml_nocancer" & in_generic == TRUE, gene],
  reml_fibroblast = all_gen[panel == "reml_fibroblast" & in_generic == TRUE, gene],
  rmamv_noc2c12 = all_gen[panel == "rmamv_noc2c12" & in_generic == TRUE, gene]
)
nm <- names(sets)
jac <- rbindlist(lapply(seq_along(nm), function(i) {
  rbindlist(lapply(seq_along(nm), function(j) {
    if (j < i) return(NULL)
    data.table(
      a = nm[i], b = nm[j],
      n_a = length(sets[[i]]), n_b = length(sets[[j]]),
      n_shared = length(intersect(sets[[i]], sets[[j]])),
      jaccard = jaccard(sets[[i]], sets[[j]])
    )
  }))
}))
fwrite(jac, file.path(RES, "cell_meta_celltype_jaccard.csv"))
cat("\nJaccard vs DL generic:\n")
print(jac[a == "dl_full" | b == "dl_full"])

counts <- all_gen[, .(
  n_sig = sum(in_generic),
  n_up = sum(direction == "up"),
  n_down = sum(direction == "down")
), by = panel]
counts <- rbind(
  data.table(panel = "dl_full", n_sig = length(old_gen),
             n_up = length(old_up), n_down = length(old_dn)),
  counts
)
fwrite(counts, file.path(RES, "cell_meta_celltype_counts.csv"))
print(counts)

# Direction flips among shared genes (DL vs REML no-C2C12)
shared <- intersect(old_gen, sets$reml_noc2c12)
flip <- merge(
  old[gene %in% shared, .(gene, dir_dl = direction)],
  all_gen[panel == "reml_noc2c12" & gene %in% shared, .(gene, dir_reml = direction)],
  by = "gene"
)
cat("Shared with REML no-C2C12:", length(shared),
    " direction flips:", sum(flip$dir_dl != flip$dir_reml), "\n")

# Plot
jac_plot <- jac[a == "dl_full" & b != "dl_full"] %>%
  mutate(
    lab = recode(
      b,
      reml_noc2c12 = "REML, drop C2C12",
      reml_nocancer = "REML, drop cancer + C2C12",
      reml_fibroblast = "REML, fibroblast-only",
      rmamv_noc2c12 = "rma.mv REML, drop C2C12"
    ),
    lab = fct_reorder(lab, jaccard)
  )
p <- ggplot(jac_plot, aes(jaccard, lab)) +
  geom_col(width = 0.62, fill = unname(treatment_colors["stress"]),
           color = "grey25", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f  (%d / %d)", jaccard, n_shared, n_b)),
            hjust = -0.05, fontface = "bold", size = 3.6) +
  coord_cartesian(xlim = c(0, 1.15)) +
  labs(
    title = "Generic-gene overlap after REML and cell-type filters",
    subtitle = "Jaccard vs the original DL 1,344-gene HMP list. Primary discovery stays the full multi-cell-type set.",
    x = "Jaccard overlap", y = NULL
  ) +
  theme_pub(13)
save_png(p, "reml_filter_jaccard.png", 8.4, 4.6)

n_plot <- counts %>%
  mutate(
    lab = recode(
      panel,
      dl_full = "DL, all cell types",
      reml_noc2c12 = "REML, drop C2C12",
      reml_nocancer = "REML, drop cancer + C2C12",
      reml_fibroblast = "REML, fibroblast-only",
      rmamv_noc2c12 = "rma.mv REML, drop C2C12"
    )
  ) %>%
  pivot_longer(c(n_up, n_down), names_to = "dir", values_to = "n") %>%
  mutate(
    dir = recode(dir, n_up = "Up", n_down = "Down"),
    lab = factor(lab, levels = c(
      "DL, all cell types", "REML, drop C2C12",
      "REML, drop cancer + C2C12", "REML, fibroblast-only",
      "rma.mv REML, drop C2C12"
    ))
  )
p2 <- ggplot(n_plot, aes(lab, n, fill = dir)) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.3) +
  scale_fill_manual(values = c(Up = unname(treatment_colors["heat"]),
                               Down = unname(treatment_colors["hypoxia"])),
                    name = NULL) +
  labs(
    title = "Generic CSR size under REML and cell-type filters",
    subtitle = "Fibroblast-only is the Drive-matched check; hypoxia is thin (3 GSE).",
    x = NULL, y = "Generic genes"
  ) +
  theme_pub(13) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "bottom")
save_png(p2, "reml_filter_sizes.png", 8.2, 5.2)

writeLines(c(
  "REML stage-2 + cell-type sensitivities (06_reml_celltype.R)",
  "",
  "Primary discovery list remains the 02_public_hmp DL HMP (1,344 genes).",
  "This script recomputes the four class mixed models with REML tau^2",
  "on the existing Deming-scaled ECs, dropping C2C12 from the primary REML panel.",
  "Sensitivities: no cancer lines; fibroblast-like cells only.",
  "",
  "Full rma.mv with crossed tissue/species intercepts is still not run."
), file.path(RES, "README_reml_celltype.txt"))

cat("Done.\n")
