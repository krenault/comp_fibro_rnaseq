#!/usr/bin/env Rscript
# Cell 2023 (Tyshkovskiy / Gladyshev, PMC11192172) meta-analysis, adapted to
# cellular stress. Drive is held out of discovery and used only as the GSEA query.
#
# STAR Methods mapping:
#   1. Per-contrast EC + SE + p  (limma in the paper; here ASTRA edgeR p, GEO
#      Wald/padj, or conservative MAD-SE when no p exists)
#   2. Denoised Spearman on the union of top-250 genes per pair
#   3. Multiple Deming scale-normalization (optim L-BFGS-B, 10 random starts)
#   4. metafor mixed-effects on mean EC + SE; random GSE, tissue, species
#      (and stress class for generic models)
#   5. Signature = BH p.adj < 0.05
#   6. HMP across the four class signatures, two-sided, ≥80% sign agreement
#   7. Fisher combination for shared vs class-distinct
#   8. Honest Drive GSEA is 03_drive_holdout.R (class hold-out), not this script.
#
# Discovery units: ASTRA record_id (not collapsed) + extra human/mouse cell GEO.
# In vivo GEO and Drive libraries are not in discovery.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(metafor)
  library(fgsea)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
source("scripts/shared/plot_aesthetics.R")
source("scripts/shared/geo_helpers.R")

OUT <- file.path(ROOT, "results/mammalian_generic_stress")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

CLASSES <- c("heat", "hypoxia", "oxidative", "dna_damage")
NCORES <- max(1L, parallel::detectCores() - 1L)
TOP_K <- 250L
SIGN_AGREE <- 0.80
HMP_PADJ <- 0.05
META_PADJ <- 0.05

# -----------------------------------------------------------------------------
# Helpers (Cell STAR Methods)
# -----------------------------------------------------------------------------
se_from_p <- function(ec, p) {
  p <- pmin(pmax(as.numeric(p), 1e-300), 1 - 1e-16)
  z <- stats::qnorm(p / 2, lower.tail = FALSE)
  se <- abs(as.numeric(ec)) / pmax(z, 1e-8)
  se[!is.finite(se) | se <= 0] <- NA_real_
  se
}

# Direction-aware harmonic mean p (Yoon / Wilson, as in the Cell paper).
# Landau tail via pgamma approximation of the HMP null for small L is
# conservative enough after BH; if harmonicmeanp is present, use it.
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
      as.numeric(harmonicmeanp::p.hmp(ps, w = rep(1 / length(ps), length(ps)), L = L, multilevel = FALSE))
    } else {
      h <- length(ps) / sum(1 / ps)
      # Wilson 2019 Landau location log(L)+γ, scale π/2; p ≈ 1−F_Landau(1/h)
      mu <- log(L) + 1 - digamma(1)
      sig <- pi / 2
      z <- (1 / h - mu) / sig
      # Landau CDF via the stable-distribution relation used by FMStable/pEstable
      # Fallback: 2*h (anti-conservative) then BH across genes.
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

fisher_twosided <- function(p, sign_ec) {
  ok <- is.finite(p) & is.finite(sign_ec) & p > 0
  p <- p[ok]; s <- sign(sign_ec[ok]); s[s == 0] <- 1
  if (length(p) < 2L) return(list(p = NA_real_, dir = 0))
  comb <- function(ps) stats::pchisq(-2 * sum(log(pmax(ps, 1e-300))),
                                     df = 2 * length(ps), lower.tail = FALSE)
  one_pos <- ifelse(s > 0, pmin(p / 2, 1 - 1e-16), pmin(1 - p / 2, 1 - 1e-16))
  one_neg <- ifelse(s < 0, pmin(p / 2, 1 - 1e-16), pmin(1 - p / 2, 1 - 1e-16))
  pp <- comb(one_pos); pn <- comb(one_neg)
  if (pp <= pn) list(p = min(1, 2 * pp), dir = 1) else list(p = min(1, 2 * pn), dir = -1)
}

# Two-stage random-effects (Cell: EC+SE into metafor; random study term).
# Stage 1: inverse-variance collapse of records within GSE (dose/time/cell
#   from the same study). Stage 2: REML rma() across GSE (DL/EE fallback).
# Full rma.mv REML with crossed tissue/species intercepts is the paper's
# model but is too slow for ~15k genes × 15 signatures; univariate REML
# on GSE-level estimates still uses SE. Cell-type sensitivities and the
# DL-vs-REML concordance live in scripts/generic_stress/06_reml_celltype.R.
fit_rma <- function(d) {
  d <- d[is.finite(yi) & is.finite(vi) & vi > 0]
  if (nrow(d) < 3L) return(NULL)
  d[, vi := pmin(pmax(vi, 1e-12), 1e4)]
  g <- d[, {
    w <- 1 / vi
    sw <- sum(w)
    yi_g <- sum(w * yi) / sw
    vi_g <- 1 / sw
    .(yi = yi_g, vi = vi_g)
  }, by = gse_id]
  if (nrow(g) < 3L) {
    g <- d[, .(yi, vi)]
  }
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

meta_by_gene <- function(dt, min_n = 3L, label = "") {
  cat("  mixed-effect", label, " contrasts=", uniqueN(dt$dataset_id),
      " genes=", uniqueN(dt$gene), "\n")
  spl <- split(dt, by = "gene", drop = TRUE)
  keep <- vapply(spl, nrow, integer(1)) >= min_n
  spl <- spl[keep]
  cat("    genes with n>=", min_n, ": ", length(spl), "  cores=", NCORES, "\n", sep = "")
  res <- parallel::mclapply(spl, fit_rma, mc.cores = NCORES)
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
  out[, significant := padj < META_PADJ]
  out[, signature := label]
  out[]
}

# -----------------------------------------------------------------------------
# 1) Per-contrast EC + SE
# -----------------------------------------------------------------------------
cat("=== Build EC + SE long table ===\n")
astra_cat <- as.data.table(read_csv("results/astra/comparison_catalog.csv", show_col_types = FALSE))
astra_cat[, record_id := as.character(record_id)]
astra <- as.data.table(readRDS("results/astra/astra_z_long.rds"))
astra[, record_id := as.character(record_id)]
astra <- merge(
  astra,
  astra_cat[, .(record_id, stress_class, experiment_id, GEO_id, Cell_type,
                Tissue, Cell_group, Genotype)],
  by = "record_id"
)
astra <- astra[!is.na(stress_class) & is.finite(logfc)]
astra[, `:=`(
  dataset_id = paste0("astra_", record_id),
  gse_id = GEO_id,
  species = "human",
  tissue = Tissue,
  cell_type = Cell_type,
  source = "ASTRA",
  yi = as.numeric(logfc),
  p_raw = as.numeric(pvalue)
)]
astra[, se := se_from_p(yi, p_raw)]
astra <- astra[is.finite(se)]

geo_disc <- tribble(
  ~file, ~dataset_id, ~stress_class, ~species,
  "GSE98906_mouse_HS_vs_UN.csv", "mouse_GSE98906_heat", "heat", "mouse",
  "GSE65636_mouse_NIH3T3_HS_45_vs_control.csv", "mouse_GSE65636_HS45", "heat", "mouse",
  "GSE65636_mouse_NIH3T3_HS_135_vs_control.csv", "mouse_GSE65636_HS135", "heat", "mouse",
  "GSE98906_mouse_H2O2_vs_UN.csv", "mouse_GSE98906_H2O2", "oxidative", "mouse",
  "GSE262772_mouse_20Gy_vs_Con_24h.csv", "mouse_GSE262772_20Gy", "dna_damage", "mouse",
  "GSE215293_MEF_WT_HS_vs_ctrl.csv", "mouse_GSE215293_heat", "heat", "mouse",
  "GSE243284_NIH3T3_siCtrl_HS_vs_ctrl.csv", "mouse_GSE243284_heat", "heat", "mouse",
  "GSE197536_NIH3T3_HS_2S_vs_ctrl.csv", "mouse_GSE197536_HS_2S", "heat", "mouse",
  "GSE197536_NIH3T3_HS_8M_vs_ctrl.csv", "mouse_GSE197536_HS_8M", "heat", "mouse",
  "GSE197536_NIH3T3_H2O2_2h_vs_ctrl.csv", "mouse_GSE197536_H2O2_2h", "oxidative", "mouse",
  "GSE197536_NIH3T3_H2O2_7h_vs_ctrl.csv", "mouse_GSE197536_H2O2_7h", "oxidative", "mouse",
  "GSE281390_NIH3T3_UVA_vs_ctrl.csv", "mouse_GSE281390_UVA", "dna_damage", "mouse",
  "GSE66286_NIH3T3_UVC_6h_vs_ctrl.csv", "mouse_GSE66286_UVC", "dna_damage", "mouse",
  "GSE114086_C2C12_GM_hypoxia_vs_normoxia.csv", "mouse_GSE114086_GM_hyp", "hypoxia", "mouse",
  "GSE114086_C2C12_DM_hypoxia_vs_normoxia.csv", "mouse_GSE114086_DM_hyp", "hypoxia", "mouse",
  "GSE171871_mESC_hypoxia_vs_normoxia.csv", "mouse_GSE171871_mESC_hyp", "hypoxia", "mouse",
  "GSE205381_MEF_hypoxia_vs_ctrl.csv", "mouse_GSE205381_hyp", "hypoxia", "mouse",
  "GSE205381_MEF_oxidative_vs_ctrl.csv", "mouse_GSE205381_ox", "oxidative", "mouse",
  "GSE281533_human_hypoxia_vs_normoxia.csv", "human_GSE281533_hyp", "hypoxia", "human",
  "GSE139963_human_hypoxia_vs_PBS.csv", "human_GSE139963_hyp", "hypoxia", "human",
  "GSE225095_IMR90_H2O2_vs_Prolif.csv", "human_GSE225095_H2O2", "oxidative", "human",
  "GSE262772_human_HDF_20Gy_vs_Con_24h.csv", "human_GSE262772_20Gy", "dna_damage", "human"
)

meta <- as.data.table(read_tsv(file.path(OUT, "discovery_metadata.tsv"), show_col_types = FALSE))

geo_long <- rbindlist(lapply(seq_len(nrow(geo_disc)), function(i) {
  r <- geo_disc[i, ]
  path <- file.path("results/geo_deseq", r$file)
  if (!file.exists(path)) {
    message("missing ", r$file)
    return(NULL)
  }
  d <- as.data.table(read_csv(path, show_col_types = FALSE))
  lfc_col <- intersect(c("log2FoldChange", "lfc", "logFC"), names(d))[1]
  p_col <- intersect(c("pvalue", "P.Value", "p_value"), names(d))[1]
  padj_col <- intersect(c("padj", "adj.P.Val", "FDR"), names(d))[1]
  se_col <- intersect(c("lfcSE", "SE", "se"), names(d))[1]
  method <- if ("method" %in% names(d)) as.character(d$method[1]) else NA_character_
  out <- data.table(
    gene = strip(d$gene),
    yi = as.numeric(d[[lfc_col]])
  )
  if (!is.na(se_col)) out$se_given <- as.numeric(d[[se_col]]) else out$se_given <- NA_real_
  if (!is.na(p_col)) out$p_raw <- as.numeric(d[[p_col]]) else out$p_raw <- NA_real_
  if (!is.na(padj_col) && all(is.na(out$p_raw))) out$p_raw <- as.numeric(d[[padj_col]])
  out <- out[is.finite(yi) & grepl("^ENSG", gene)]
  out[, se := fifelse(is.finite(se_given) & se_given > 0, se_given, se_from_p(yi, p_raw))]
  if (mean(is.finite(out$se)) < 0.5) {
    mad_se <- max(stats::mad(out$yi, na.rm = TRUE), 0.35)
    out[!is.finite(se), se := mad_se]
  }
  # author-sig-only / n=1: inflate SE so they cannot dominate ASTRA
  meth <- ifelse(is.na(method) || !nzchar(method), "", method)
  if (grepl("author_DESeq|log2abund_n1|log2FPKM|log2RPKM|log2Cuff", meth)) {
    out[, se := se * 1.5]
  }
  if (r$dataset_id %in% c("mouse_GSE205381_hyp", "mouse_GSE205381_ox", "mouse_GSE66286_UVC")) {
    out[, se := se * 1.5]
  }
  out[, `:=`(
    dataset_id = r$dataset_id, stress_class = r$stress_class, species = r$species,
    gse_id = sub(".*(GSE[0-9]+).*", "\\1", r$file),
    source = "GEO", se_given = NULL
  )]
  out[is.finite(se) & se > 0]
}), fill = TRUE)

# tissue from metadata (GEO drivers) or ASTRA Tissue
geo_meta <- meta[source == "GEO", .(dataset_id = experiment_id, tissue = tissue_of_origin,
                                    cell_type)]
geo_long <- merge(geo_long, geo_meta, by = "dataset_id", all.x = TRUE)
geo_long[is.na(tissue), tissue := "unknown"]

long <- rbindlist(list(
  astra[, .(gene, dataset_id, gse_id, species, tissue, cell_type, stress_class,
            source, yi, se, p_raw)],
  geo_long[, .(gene, dataset_id, gse_id, species, tissue, cell_type, stress_class,
               source, yi, se, p_raw)]
), fill = TRUE)
long[, vi := se^2]
long <- long[is.finite(yi) & is.finite(vi) & vi > 0]
cat("Contrasts:", uniqueN(long$dataset_id),
    "  genes×rows:", nrow(long),
    "  human:", uniqueN(long[species == "human", dataset_id]),
    "  mouse:", uniqueN(long[species == "mouse", dataset_id]), "\n")
fwrite(unique(long[, .(dataset_id, gse_id, species, tissue, cell_type, stress_class, source)]),
       file.path(OUT, "cell_meta_contrasts.csv"))

# -----------------------------------------------------------------------------
# 2) Denoised pairwise Spearman (top 250)
# -----------------------------------------------------------------------------
cat("=== Denoised pairwise Spearman ===\n")
ds <- sort(unique(long$dataset_id))
wide_p <- dcast(long, gene ~ dataset_id, value.var = "p_raw")
wide_y <- dcast(long, gene ~ dataset_id, value.var = "yi")
genes <- wide_y$gene
Y <- as.matrix(wide_y[, -1])
P <- as.matrix(wide_p[, -1])
storage.mode(Y) <- "double"
# rank genes within dataset by p (missing p → |yi| rank)
rank_mat <- vapply(seq_along(ds), function(j) {
  pj <- P[, j]; yj <- Y[, j]
  sc <- ifelse(is.finite(pj), pj, 1 / (1e-6 + abs(yj)))
  sc[!is.finite(yj)] <- NA_real_
  rank(sc, na.last = "keep", ties.method = "average")
}, numeric(nrow(Y)))
colnames(rank_mat) <- ds

pair_n <- length(ds) * (length(ds) - 1L) / 2L
pair_i <- integer(pair_n); pair_j <- integer(pair_n)
rho <- numeric(pair_n); pv <- numeric(pair_n)
k <- 0L
for (a in seq_len(length(ds) - 1L)) {
  for (b in (a + 1L):length(ds)) {
    k <- k + 1L
    pair_i[k] <- a; pair_j[k] <- b
    top <- which(rank_mat[, a] <= TOP_K | rank_mat[, b] <= TOP_K)
    xa <- Y[top, a]; xb <- Y[top, b]
    ok <- is.finite(xa) & is.finite(xb)
    if (sum(ok) < 40L) {
      rho[k] <- NA_real_; pv[k] <- NA_real_; next
    }
    ct <- suppressWarnings(cor.test(xa[ok], xb[ok], method = "spearman", exact = FALSE))
    rho[k] <- unname(ct$estimate); pv[k] <- ct$p.value
  }
}
pairs <- data.table(
  dataset_i = ds[pair_i], dataset_j = ds[pair_j],
  rho = rho, p = pv
)
pairs[, padj := p.adjust(p, method = "BH")]
pairs[, keep_deming := is.finite(rho) & padj < 0.05 & abs(rho) > 0.1]
fwrite(pairs, file.path(OUT, "cell_meta_pairwise_spearman.csv"))
cat("Pairs:", nrow(pairs), "  Deming-eligible:", sum(pairs$keep_deming, na.rm = TRUE),
    "  median rho:", median(pairs$rho, na.rm = TRUE), "\n")

# -----------------------------------------------------------------------------
# 3) Multiple Deming scale (geomean s = 1, s1 fixed via log-sum)
# -----------------------------------------------------------------------------
cat("=== Multiple Deming ===\n")
keepp <- pairs[keep_deming == TRUE]
if (nrow(keepp) >= 10L) {
  ds_index <- setNames(seq_along(ds), ds)
  pair_chunks <- vector("list", nrow(keepp))
  for (r in seq_len(nrow(keepp))) {
    a <- ds_index[[keepp$dataset_i[r]]]
    b <- ds_index[[keepp$dataset_j[r]]]
    top <- which(rank_mat[, a] <= TOP_K | rank_mat[, b] <= TOP_K)
    xa <- Y[top, a]; xb <- Y[top, b]
    ok <- is.finite(xa) & is.finite(xb)
    if (!any(ok)) next
    pair_chunks[[r]] <- data.table(xi = a, xj = b, yi = xa[ok], yj = xb[ok])
  }
  pair_dt <- rbindlist(pair_chunks, fill = TRUE)
  xi <- pair_dt$xi; xj <- pair_dt$xj; yi <- pair_dt$yi; yj <- pair_dt$yj
  n_ds <- length(ds)
  loss <- function(par) {
    s <- numeric(n_ds)
    s[-1] <- par
    s[1] <- 1
    sum((s[xi] * yi - s[xj] * yj)^2)
  }
  lower <- rep(0.01, n_ds - 1L)
  upper <- rep(100, n_ds - 1L)
  best <- NULL
  best_val <- Inf
  set.seed(1)
  inits <- vector("list", 10L)
  inits[[1]] <- rep(1, n_ds - 1L)
  mad_s <- apply(Y, 2, function(z) {
    m <- stats::mad(z, na.rm = TRUE)
    if (!is.finite(m) || m < 1e-6) 1 else 1 / m
  })
  mad_s <- mad_s / mad_s[1]
  inits[[2]] <- pmin(pmax(mad_s[-1], 0.01), 100)
  for (t in 3:10) inits[[t]] <- 10^runif(n_ds - 1L, -1.5, 1.5)
  for (t in seq_along(inits)) {
    fit <- tryCatch(
      stats::optim(inits[[t]], loss, method = "L-BFGS-B",
                   lower = lower, upper = upper,
                   control = list(maxit = 200, factr = 1e9)),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    if (fit$value < best_val) {
      best_val <- fit$value
      best <- fit$par
    }
    cat("  start", t, " loss=", signif(fit$value, 5), "\n", sep = "")
  }
  s <- rep(1, n_ds); names(s) <- ds
  if (!is.null(best)) {
    s[-1] <- best
    s <- s / exp(mean(log(s)))
  }
  scale_tbl <- data.table(dataset_id = ds, deming_s = as.numeric(s[ds]))
} else {
  message("Too few correlated pairs; MAD scaling only")
  scale_tbl <- data.table(
    dataset_id = ds,
    deming_s = apply(Y, 2, function(z) {
      m <- stats::mad(z, na.rm = TRUE)
      if (!is.finite(m) || m < 1e-6) 1 else 1 / m
    })
  )
  scale_tbl[, deming_s := deming_s / exp(mean(log(deming_s)))]
}
fwrite(scale_tbl, file.path(OUT, "cell_meta_deming_scales.csv"))
long <- merge(long, scale_tbl, by = "dataset_id", all.x = TRUE)
long[is.na(deming_s), deming_s := 1]
long[, `:=`(yi = yi * deming_s, se = se * deming_s, vi = (se)^2)]
saveRDS(long, file.path(OUT, "cell_meta_ec_long.rds"))

# -----------------------------------------------------------------------------
# 4) Mixed-effect signatures
# -----------------------------------------------------------------------------
cat("=== Mixed-effect signatures ===\n")
sig_list <- list()

for (cl in CLASSES) {
  sig_list[[paste0("class_", cl)]] <- meta_by_gene(
    long[stress_class == cl], min_n = 3L, label = paste0("class_", cl)
  )
}
sig_list[["human"]] <- meta_by_gene(
  long[species == "human"], min_n = 5L, label = "human"
)
sig_list[["mouse"]] <- meta_by_gene(
  long[species == "mouse"], min_n = 3L, label = "mouse"
)
sig_list[["global"]] <- meta_by_gene(long, min_n = 5L, label = "global")

for (sp in c("human", "mouse")) {
  for (cl in CLASSES) {
    lab <- paste(sp, cl, sep = "_")
    sig_list[[lab]] <- meta_by_gene(
      long[species == sp & stress_class == cl], min_n = 3L, label = lab
    )
  }
}

sigs <- rbindlist(sig_list, fill = TRUE)
fwrite(sigs, file.path(OUT, "cell_meta_signatures.csv"))
cat("Signature sizes (padj<0.05):\n")
print(sigs[, .(n_sig = sum(significant), n_up = sum(significant & ec > 0),
               n_dn = sum(significant & ec < 0)), by = signature])

# -----------------------------------------------------------------------------
# 5) HMP generic (four class signatures) + 80% sign agreement
# -----------------------------------------------------------------------------
cat("=== HMP generic stress ===\n")
class_wide <- dcast(
  sigs[signature %in% paste0("class_", CLASSES), .(gene, signature, p, ec, padj)],
  gene ~ signature, value.var = c("p", "ec", "padj")
)
p_cols <- paste0("p_class_", CLASSES)
ec_cols <- paste0("ec_class_", CLASSES)
hmp_rows <- lapply(seq_len(nrow(class_wide)), function(i) {
  pv <- unlist(class_wide[i, ..p_cols], use.names = FALSE)
  ev <- unlist(class_wide[i, ..ec_cols], use.names = FALSE)
  h <- hmp_twosided(pv, ev, L = 4L)
  data.table(
    gene = class_wide$gene[i],
    hmp_p = h$p, hmp_dir = h$dir, n_class = h$n, frac_agree = h$frac_agree,
    ec_mean = mean(ev[is.finite(ev)])
  )
})
generic <- rbindlist(hmp_rows)
generic[, hmp_padj := p.adjust(hmp_p, method = "BH")]
generic[, in_generic := is.finite(hmp_padj) & hmp_padj < HMP_PADJ &
          n_class >= 3L & frac_agree >= SIGN_AGREE]
generic[, direction := fcase(
  in_generic & hmp_dir > 0, "up",
  in_generic & hmp_dir < 0, "down",
  default = "none"
)]

# species-generic: HMP of human_* and mouse_* class signatures
species_hmp <- function(sp) {
  labs <- paste(sp, CLASSES, sep = "_")
  w <- dcast(sigs[signature %in% labs, .(gene, signature, p, ec)],
             gene ~ signature, value.var = c("p", "ec"))
  pc <- paste0("p_", labs); ec <- paste0("ec_", labs)
  rbindlist(lapply(seq_len(nrow(w)), function(i) {
    pv <- unlist(w[i, ..pc], use.names = FALSE)
    ev <- unlist(w[i, ..ec], use.names = FALSE)
    h <- hmp_twosided(pv, ev, L = 4L)
    data.table(gene = w$gene[i], p = h$p, dir = h$dir, n = h$n,
               frac_agree = h$frac_agree, ec_mean = mean(ev, na.rm = TRUE))
  }))
}
hum_g <- species_hmp("human")
mou_g <- species_hmp("mouse")
hum_g[, padj := p.adjust(p, method = "BH")]
mou_g[, padj := p.adjust(p, method = "BH")]
hum_g[, in_sp := padj < HMP_PADJ & n >= 3L & frac_agree >= SIGN_AGREE]
mou_g[, in_sp := padj < HMP_PADJ & n >= 3L & frac_agree >= SIGN_AGREE]

mam <- merge(
  hum_g[, .(gene, p_h = p, padj_h = padj, dir_h = dir, agree_h = frac_agree,
            n_h = n, ec_h = ec_mean, human_generic = in_sp)],
  mou_g[, .(gene, p_m = p, padj_m = padj, dir_m = dir, agree_m = frac_agree,
            n_m = n, ec_m = ec_mean, mouse_generic = in_sp)],
  by = "gene", all = TRUE
)
mam <- merge(mam, generic[, .(gene, hmp_p, hmp_padj, hmp_dir, frac_agree,
                              n_class, in_generic, direction, ec_mean)],
             by = "gene", all = TRUE)
mam[, same_dir := is.finite(dir_h) & is.finite(dir_m) & dir_h == dir_m & dir_h != 0]
mam[, in_mammalian := human_generic == TRUE & mouse_generic == TRUE & same_dir == TRUE]
# primary Cell-style generic = HMP of the four class mixed models (both species)
# plus the stricter both-species call
mam[is.na(in_generic), in_generic := FALSE]
mam[is.na(in_mammalian), in_mammalian := FALSE]

syms <- as.data.table(ensg_to_symbol(unique(mam$gene)))
mam <- merge(mam, syms, by = "gene", all.x = TRUE)

# class-specific: class mixed-model significant, not generic, same sign in both species
class_sp <- rbindlist(lapply(CLASSES, function(cl) {
  h <- sigs[signature == paste0("human_", cl), .(gene, ec_h = ec, p_h = p, padj_h = padj, sig_h = significant)]
  m <- sigs[signature == paste0("mouse_", cl), .(gene, ec_m = ec, p_m = p, padj_m = padj, sig_m = significant)]
  both <- merge(h, m, by = "gene")
  both[, stress_class := cl]
  both[, in_class_mammalian := sig_h == TRUE & sig_m == TRUE & sign(ec_h) == sign(ec_m) & ec_h != 0]
  both[]
}))
class_sp <- merge(class_sp, mam[, .(gene, in_generic, in_mammalian)], by = "gene", all.x = TRUE)
class_sp[is.na(in_generic), in_generic := FALSE]
class_sp[, in_class_specific := in_class_mammalian == TRUE & in_generic == FALSE]
class_sp <- merge(class_sp, syms, by = "gene", all.x = TRUE)

# Fisher: generic vs each class (shared / opposite)
fisher_tbl <- rbindlist(lapply(CLASSES, function(cl) {
  s <- sigs[signature == paste0("class_", cl), .(gene, p_cl = p, ec_cl = ec)]
  g <- generic[, .(gene, p_g = hmp_p, ec_g = ec_mean)]
  x <- merge(s, g, by = "gene")
  sh <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
    f <- fisher_twosided(c(x$p_g[i], x$p_cl[i]), c(x$ec_g[i], x$ec_cl[i]))
    opp <- fisher_twosided(c(x$p_g[i], x$p_cl[i]), c(x$ec_g[i], -x$ec_cl[i]))
    data.table(gene = x$gene[i], p_shared = f$p, dir_shared = f$dir,
               p_opposite = opp$p, dir_opposite = opp$dir)
  }))
  sh[, stress_class := cl]
  sh[, padj_shared := p.adjust(p_shared, method = "BH")]
  sh[, padj_opposite := p.adjust(p_opposite, method = "BH")]
  sh[]
}))

fwrite(sigs, file.path(OUT, "cell_meta_signatures.csv"))
fwrite(mam, file.path(OUT, "gene_human_mouse_generic.csv"))
fwrite(mam[in_generic == TRUE], file.path(OUT, "mammalian_generic_genes.csv"))
fwrite(mam[in_mammalian == TRUE], file.path(OUT, "mammalian_generic_both_species.csv"))
fwrite(class_sp, file.path(OUT, "gene_class_human_mouse.csv"))
fwrite(class_sp[in_class_specific == TRUE], file.path(OUT, "class_specific_mammalian_genes.csv"))
fwrite(fisher_tbl, file.path(OUT, "cell_meta_fisher_shared_distinct.csv"))

set_counts <- rbindlist(list(
  mam[, .(n = .N), by = .(in_generic, direction)][, set := "generic_hmp"],
  mam[, .(n = sum(in_mammalian == TRUE), set = "generic_both_species")],
  class_sp[in_class_specific == TRUE, .(n = .N), by = stress_class][, set := "class_specific"]
), fill = TRUE)
fwrite(set_counts, file.path(OUT, "set_counts.csv"))
cat("Generic HMP:", sum(mam$in_generic),
    " up", sum(mam$direction == "up"), " down", sum(mam$direction == "down"), "\n")
cat("Both-species generic:", sum(mam$in_mammalian, na.rm = TRUE), "\n")
cat("Class-specific:\n")
print(class_sp[in_class_specific == TRUE, .N, by = stress_class])

# -----------------------------------------------------------------------------
# 6) Ranked generic list (Hallmark + Reactome + GO BP GSEA is written below)
# -----------------------------------------------------------------------------
rank_generic <- mam[!is.na(hmp_p) & is.finite(ec_mean)]
rank_generic[, stat := -log10(pmax(hmp_p, 1e-300)) * sign(ec_mean)]
rvec <- setNames(rank_generic$stat, rank_generic$gene)
rvec <- rvec[is.finite(rvec)]
rvec <- rvec[!duplicated(names(rvec))]

# -----------------------------------------------------------------------------
# 7) Drive conserved assignment (held out). Honest GSEA is 03_drive_holdout.R.
# -----------------------------------------------------------------------------
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
drive_long <- drive_long[species != "rousette"]
drive_med <- drive_long[, .(
  z_drive = median(lfc / stats::mad(lfc, na.rm = TRUE), na.rm = TRUE),
  n_sp = uniqueN(species),
  frac_agree = frac_agree_fun(lfc)
), by = .(gene, exposure)]
# per-exposure, scale within species already mixed; use median LFC
drive_med2 <- drive_long[, .(
  lfc_med = median(lfc), n_sp = uniqueN(species),
  frac_agree = frac_agree_fun(lfc)
), by = .(gene, exposure)]
cmp <- merge(drive_med2, mam[, .(gene, ec_mean, in_generic, direction, symbol, hmp_p)],
             by = "gene")
rho_tbl <- cmp[, {
  a <- safe_spearman(lfc_med, ec_mean)
  g <- safe_spearman(lfc_med[in_generic == TRUE], ec_mean[in_generic == TRUE], min_n = 20L)
  data.table(n_all = a$n, rho_all = a$rho, p_all = a$p,
             n_generic = g$n, rho_generic = g$rho, p_generic = g$p)
}, by = exposure]
fwrite(rho_tbl, file.path(OUT, "drive_vs_generic_spearman.csv"))

drive_cons <- drive_med2[n_sp >= 5 & frac_agree >= 0.8 & abs(lfc_med) >= 0.25]
drive_cons[, in_generic := gene %in% mam$gene[mam$in_generic == TRUE]]
drive_cons[, matched_class := fcase(
  exposure == "heat", gene %in% class_sp$gene[class_sp$stress_class == "heat" & class_sp$in_class_specific == TRUE],
  exposure %in% c("hypoxia_6H", "hypoxia_24H"),
  gene %in% class_sp$gene[class_sp$stress_class == "hypoxia" & class_sp$in_class_specific == TRUE],
  default = FALSE
)]
drive_cons[, assignment := fcase(
  in_generic == TRUE, "generic_mammalian_stress",
  matched_class == TRUE, "matched_class_specific",
  default = "drive_residual"
)]
drive_cons <- merge(drive_cons, syms, by = "gene", all.x = TRUE)
fwrite(drive_cons, file.path(OUT, "drive_conserved_assignment.csv"))
asg <- drive_cons[, .N, by = .(exposure, assignment)]
setnames(asg, "N", "n")
fwrite(asg, file.path(OUT, "drive_conserved_assignment_counts.csv"))
cat("Drive conserved assignment:\n")
print(asg)

focal <- c("FASN", "ACACA", "SREBF1", "SREBF2", "SCD", "HSPA1A", "HSPA1B",
           "HIF1A", "VEGFA", "DDIT3", "ATF3", "ATF4", "CDKN1A", "GADD45A",
           "HMOX1", "SOD2", "TXN")
focal_tab <- mam[symbol %in% focal, .(
  gene, symbol, ec_h, ec_m, ec_mean, hmp_p, hmp_padj, in_generic, in_mammalian, direction
)]
fwrite(focal_tab, file.path(OUT, "focal_genes.csv"))

# Plots are drawn by scripts/generic_stress/04_plots.R.
# Combined Hallmark + Reactome + GO BP GSEA (figure public_pathways).
gsea_combo <- file.path(OUT, "cell_meta_gsea_hallmark_reactome_gobp.csv")
if (!file.exists(gsea_combo) && length(rvec) > 200) {
  msig_sets <- function(collection, subcollection = NULL) {
    x <- tryCatch({
      if (is.null(subcollection)) msigdbr::msigdbr(species = "Homo sapiens", collection = collection)
      else msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = subcollection)
    }, error = function(e) {
      if (is.null(subcollection)) msigdbr::msigdbr(species = "Homo sapiens", category = collection)
      else msigdbr::msigdbr(species = "Homo sapiens", category = collection, subcategory = subcollection)
    })
    gcol <- intersect(c("ensembl_gene", "human_ensembl_gene"), names(x))[1]
    split(strip(x[[gcol]]), x$gs_name)
  }
  pw <- c(msig_sets("H"), msig_sets("C2", "CP:REACTOME"), msig_sets("C5", "GO:BP"))
  set.seed(1)
  fg_all <- as.data.table(fgsea::fgseaMultilevel(
    pathways = pw, stats = rvec, minSize = 15, maxSize = 400
  ))
  fg_all[, collection := fcase(
    startsWith(pathway, "HALLMARK"), "Hallmark",
    startsWith(pathway, "REACTOME"), "Reactome",
    startsWith(pathway, "GOBP"), "GO BP",
    default = "Other"
  )]
  fwrite(fg_all[, .(collection, pathway, pval, padj, NES, size)], gsea_combo)
}

n_generic <- sum(mam$in_generic, na.rm = TRUE)
n_both <- sum(mam$in_mammalian, na.rm = TRUE)
cat("Done. Generic HMP n=", n_generic, " both-species n=", n_both, "\n", sep = "")
cat("Outputs in results/mammalian_generic_stress\n")
