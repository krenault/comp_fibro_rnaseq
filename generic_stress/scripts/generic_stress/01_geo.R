#!/usr/bin/env Rscript
# Download and DE mouse *cell* stress RNA-seq with usable count/FPKM tables.
# Writes results/geo_deseq/ in the same schema as existing GEO DGE files.
# Drive is not used.

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(babelgene)
})

.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
source("scripts/shared/geo_helpers.R")

RAW <- file.path(ROOT, "raw/mouse_cell_stress")
OUT <- file.path(ROOT, "results/geo_deseq")
dir.create(RAW, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

write_hum <- function(res, map, dataset, panel, comparison, method) {
  hum <- collapse_to_human(res, map)
  if (!nrow(hum)) {
    message("EMPTY after ortholog map: ", dataset)
    return(invisible(NULL))
  }
  write_geo_deseq(hum, file.path(OUT, paste0(dataset, ".csv")),
                  dataset, panel, comparison, method)
}

run_one <- function(label, expr) {
  cat("\n======== ", label, " ========\n", sep = "")
  tryCatch(expr, error = function(e) {
    message("FAILED ", label, ": ", conditionMessage(e))
  })
}

# -----------------------------------------------------------------------------
run_one("GSE215293 WT MEF heat 43C (FPKM)", {
  f <- download_geo_suppl("GSE215293", "GSE215293_gene_fpkm.txt.gz", RAW)
  d <- read_tsv(f, show_col_types = FALSE)
  gn <- if ("gene_id" %in% names(d)) "gene_id" else names(d)[1]
  keep <- grep("^MEF[0-9]+(_HS)?$", names(d), value = TRUE)
  mat <- as.matrix(d[, keep, drop = FALSE])
  rownames(mat) <- strip(d[[gn]])
  ctrl <- grep("^MEF[0-9]+$", colnames(mat), value = TRUE)
  hs <- grep("^MEF[0-9]+_HS$", colnames(mat), value = TRUE)
  if (!length(ctrl) || !length(hs)) {
    message("cols: ", paste(colnames(mat), collapse = ", "))
    stop("could not find WT MEF / MEF-HS columns")
  }
  res <- lfc_from_abundance(mat, hs, ctrl)
  map <- map_mouse_ids_to_ensg(res$gene_raw)
  write_hum(res, map, "GSE215293_MEF_WT_HS_vs_ctrl", "heat", "HS43C_1h_vs_37C", "log2FPKM_LFC")
})

# -----------------------------------------------------------------------------
run_one("GSE243284 NIH3T3 siCtrl heat (featureCounts)", {
  f <- download_geo_suppl("GSE243284", "GSE243284_BB3mRNA-all.featCnt.mm10.txt.gz", RAW)
  d <- read_tsv(f, comment = "#", show_col_types = FALSE)
  gene_col <- intersect(c("Geneid", "GeneID", "gene_id", "gene"), names(d))[1]
  sample_cols <- setdiff(names(d), c(gene_col, "Chr", "Start", "End", "Strand", "Length"))
  mat <- as.matrix(d[, sample_cols, drop = FALSE])
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat <- round(mat)
  rownames(mat) <- strip(d[[gene_col]])
  cn <- colnames(mat)
  # Sample descriptions BB3-1..8: odd siCtrl pairs nohs/hs by title mapping
  # Prefer column names containing nohs/hs or siCtrl
  message("count columns: ", paste(cn, collapse = " | "))
  pick <- function(pat) which(grepl(pat, cn, ignore.case = TRUE))
  ctrl_i <- pick("nohs.*siCtrl|siCtrl.*nohs|no hs.*siCtrl")
  hs_i <- pick("(^|[^a-z])hs, siCtrl|siCtrl.*[^n]hs|hs_siCtrl|siCtrl_hs")
  if (!length(ctrl_i) || !length(hs_i)) {
    # BB3-1,3,5,7 from series matrix: 1 and 5 nohs siCtrl; 3 and 7 hs siCtrl
    # FeatureCounts often uses BAM basenames. Fall back to 1,5 vs 3,7 if 8 cols.
    if (length(cn) == 8) {
      ctrl_i <- c(1L, 5L)
      hs_i <- c(3L, 7L)
    } else stop("cannot identify siCtrl HS vs noHS columns")
  }
  group <- rep(NA_character_, ncol(mat))
  group[ctrl_i] <- "control"
  group[hs_i] <- "HS"
  keep <- !is.na(group)
  res <- run_deseq_pair(mat[, keep, drop = FALSE], group[keep], "HS")
  map <- map_mouse_ids_to_ensg(res$gene_raw)
  write_hum(res, map, "GSE243284_NIH3T3_siCtrl_HS_vs_ctrl", "heat", "HS45C_10min_vs_37C", "DESeq2")
})

# -----------------------------------------------------------------------------
run_one("GSE197536 NIH3T3 RNA-seq heat + H2O2 (RSEM)", {
  f <- download_geo_suppl("GSE197536", "GSE197536_Rnaseq.rsem.genes.expected_count.txt.gz", RAW)
  d <- read_tsv(f, show_col_types = FALSE)
  gene_col <- intersect(c("Gene_Id", "gene_id", "Geneid", "id"), names(d))[1]
  mat <- as.matrix(d[, setdiff(names(d), gene_col), drop = FALSE])
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat <- round(mat)
  rownames(mat) <- sub("_.*$", "", strip(d[[gene_col]]))
  cn <- colnames(mat)
  message("RSEM columns: ", paste(cn, collapse = ", "))
  contrasts <- list(
    list(treat = "^hs_2S_", ctrl = "^hs_ctrl_", name = "GSE197536_NIH3T3_HS_2S_vs_ctrl",
         panel = "heat", cmp = "HS_2S_vs_ctrl"),
    list(treat = "^hs_8M_", ctrl = "^hs_ctrl_", name = "GSE197536_NIH3T3_HS_8M_vs_ctrl",
         panel = "heat", cmp = "HS_8M_vs_ctrl"),
    list(treat = "^h2o2_2h_", ctrl = "^h2o2_kcl_ctrl_", name = "GSE197536_NIH3T3_H2O2_2h_vs_ctrl",
         panel = "oxidative", cmp = "H2O2_2h_vs_ctrl"),
    list(treat = "^h2o2_7h_", ctrl = "^h2o2_kcl_ctrl_", name = "GSE197536_NIH3T3_H2O2_7h_vs_ctrl",
         panel = "oxidative", cmp = "H2O2_7h_vs_ctrl")
  )
  map <- map_mouse_ids_to_ensg(rownames(mat))
  for (ct in contrasts) {
    ti <- grepl(ct$treat, cn)
    ci <- grepl(ct$ctrl, cn)
    if (sum(ti) < 2 || sum(ci) < 2) {
      message("skip ", ct$name, " treat=", sum(ti), " ctrl=", sum(ci))
      next
    }
    group <- ifelse(ti, "treat", ifelse(ci, "control", NA_character_))
    keep <- !is.na(group)
    res <- run_deseq_pair(mat[, keep, drop = FALSE], group[keep], "treat")
    write_hum(res, map, ct$name, ct$panel, ct$cmp, "DESeq2")
  }
})

# -----------------------------------------------------------------------------
run_one("GSE281390 NIH3T3 UVA (FPKM)", {
  f <- download_geo_suppl("GSE281390", "GSE281390_all.fpkm_anno.txt.gz", RAW)
  d <- read_tsv(f, show_col_types = FALSE)
  message("cols: ", paste(names(d), collapse = ", "))
  gene_col <- if ("gene_id" %in% names(d)) "gene_id" else names(d)[1]
  fpkm <- grep("_FPKM$", names(d), value = TRUE)
  mat <- as.matrix(d[, fpkm, drop = FALSE])
  rownames(mat) <- strip(d[[gene_col]])
  ctrl <- grep("^NC[0-9]+_FPKM$", colnames(mat), value = TRUE)
  uva <- grep("^MOCK[0-9]+_FPKM$", colnames(mat), value = TRUE)
  if (!length(ctrl) || !length(uva)) stop("cannot find NC vs MOCK FPKM columns")
  res <- lfc_from_abundance(mat, uva, ctrl)
  map <- map_mouse_ids_to_ensg(res$gene_raw)
  write_hum(res, map, "GSE281390_NIH3T3_UVA_vs_ctrl", "dna_damage", "UVA_10Jcm2_vs_ctrl", "log2FPKM_LFC")
})

# -----------------------------------------------------------------------------
run_one("GSE66286 NIH3T3 UVC (author DESeq)", {
  f <- download_geo_suppl("GSE66286", "GSE66286_DESeq_0h-6h_UV.csv.gz", RAW)
  d <- read_tsv(f, show_col_types = FALSE)
  message("cols: ", paste(names(d), collapse = ", "))
  gene_col <- intersect(c("id", "gene", "Gene", "gene_id"), names(d))[1]
  lfc_col <- intersect(c("log2FoldChange", "log2FC", "logFC"), names(d))[1]
  padj_col <- intersect(c("padj", "FDR", "adj.P.Val"), names(d))[1]
  res <- tibble(
    gene_raw = strip(d[[gene_col]]),
    log2FoldChange = as.numeric(d[[lfc_col]]),
    padj = if (!is.na(padj_col)) as.numeric(d[[padj_col]]) else NA_real_
  ) %>% filter(is.finite(log2FoldChange))
  map <- map_mouse_ids_to_ensg(res$gene_raw)
  write_hum(res, map, "GSE66286_NIH3T3_UVC_6h_vs_ctrl", "dna_damage", "UVC_80Jm2_6h_vs_ctrl", "author_DESeq")
})

# -----------------------------------------------------------------------------
run_one("GSE114086 C2C12 hypoxia in GM and DM (counts)", {
  f <- download_geo_suppl("GSE114086", "GSE114086_Cuff_Gene_Counts.csv.gz", RAW)
  d <- read_csv(f, show_col_types = FALSE)
  message("cols: ", paste(names(d), collapse = ", "))
  gene_col <- names(d)[1]
  mat <- as.matrix(d[, -1, drop = FALSE])
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat <- round(mat)
  rownames(mat) <- strip(d[[gene_col]])
  cn <- colnames(mat)
  map <- map_mouse_ids_to_ensg(rownames(mat))
  for (media in c("GM", "DM")) {
    nor <- grep(paste0("^", media, "_NOR_"), cn, value = TRUE)
    hyp <- grep(paste0("^", media, "_HYP_"), cn, value = TRUE)
    if (length(nor) < 2 || length(hyp) < 2) {
      message(media, " skip nor=", length(nor), " hyp=", length(hyp))
      next
    }
    res <- lfc_from_abundance(mat, hyp, nor)
    write_hum(res, map,
              paste0("GSE114086_C2C12_", media, "_hypoxia_vs_normoxia"),
              "hypoxia", paste0(media, "_2pctO2_vs_21pct"), "log2Cuff_LFC")
  }
})

# -----------------------------------------------------------------------------
run_one("GSE171871 mESC hypoxia (counts)", {
  f <- download_geo_suppl("GSE171871", "GSE171871_readcount_matrix.txt.gz", RAW)
  d <- read_tsv(f, show_col_types = FALSE)
  if ("...1" %in% names(d)) d <- d %>% dplyr::select(-`...1`)
  gene_col <- if ("gene_id" %in% names(d)) "gene_id" else names(d)[1]
  sample_cols <- setdiff(names(d), c(gene_col, "Gene.Name", "gene_name"))
  mat <- as.matrix(d[, sample_cols, drop = FALSE])
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat <- round(mat)
  rownames(mat) <- strip(d[[gene_col]])
  cn <- colnames(mat)
  message("cols: ", paste(cn, collapse = ", "))
  hyp <- grep("H4", cn, value = TRUE)
  nor <- grep("N4", cn, value = TRUE)
  if (length(hyp) < 2 || length(nor) < 2) stop("cannot split H vs N")
  group <- ifelse(cn %in% hyp, "hypoxia", ifelse(cn %in% nor, "control", NA_character_))
  keep <- !is.na(group)
  res <- run_deseq_pair(mat[, keep, drop = FALSE], group[keep], "hypoxia")
  map <- map_mouse_ids_to_ensg(res$gene_raw)
  write_hum(res, map, "GSE171871_mESC_hypoxia_vs_normoxia", "hypoxia", "H4_vs_N4", "DESeq2")
})

# -----------------------------------------------------------------------------
run_one("GSE205381 MEF panel n=1 (abundance LFC)", {
  f <- download_geo_suppl("GSE205381", "GSE205381_MEF_stress_expression.csv.gz", RAW)
  d <- read_csv(f, show_col_types = FALSE)
  gene_col <- names(d)[grepl("ENSMUS|gene", names(d), ignore.case = TRUE)][1]
  mat <- as.matrix(d[, setdiff(names(d), c(gene_col, "Gene name", "Gene.Name")), drop = FALSE])
  rownames(mat) <- strip(d[[gene_col]])
  if (!"MEF_ctr" %in% colnames(mat)) stop("no MEF_ctr")
  map <- map_mouse_ids_to_ensg(rownames(mat))
  pairs <- list(
    list(col = "MEF_hypoxia_10h", name = "GSE205381_MEF_hypoxia_vs_ctrl", panel = "hypoxia", cmp = "hypoxia_10h_vs_ctrl"),
    list(col = "MEF_Oxidative_3h", name = "GSE205381_MEF_oxidative_vs_ctrl", panel = "oxidative", cmp = "oxidative_3h_vs_ctrl")
  )
  for (p in pairs) {
    if (!p$col %in% colnames(mat)) next
    res <- lfc_from_abundance(mat, p$col, "MEF_ctr")
    write_hum(res, map, p$name, p$panel, p$cmp, "log2abund_n1")
  }
})

cat("\nDone. New mouse files:\n")
print(list.files(OUT, pattern = "GSE(215293|243284|197536|281390|66286|114086|171871|205381)"))
