#!/usr/bin/env Rscript
# Download and DE mouse *cell* stress RNA-seq with usable count/FPKM tables.
# Writes results/geo_deseq/ in the same schema as existing GEO DGE files.
# The FH experiments are not used.

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
  only <- Sys.getenv("GEO_ONLY", unset = "")
  if (nzchar(only) && !grepl(only, label)) {
    return(invisible(NULL))
  }
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

# Per-sample counts (symbol, count). The series-level DESeq file is significant
# genes only and is not used.
run_one("GSE66286 NIH3T3 UVC (counts)", {
  specs <- tibble::tribble(
    ~gsm, ~file, ~sample, ~group,
    "GSM1618725", "GSM1618725_0h_R1.txt.gz", "ctrl_1", "control",
    "GSM1618726", "GSM1618726_0h_R2.txt.gz", "ctrl_2", "control",
    "GSM1618727", "GSM1618727_0h_R3.txt.gz", "ctrl_3", "control",
    "GSM1618728", "GSM1618728_6h_UV_R1.txt.gz", "uv_1", "uv",
    "GSM1618729", "GSM1618729_6h_UV_R2.txt.gz", "uv_2", "uv",
    "GSM1618730", "GSM1618730_6h_UV_R3.txt.gz", "uv_3", "uv"
  )
  dest_dir <- file.path(RAW, "GSE66286_counts")
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  cols <- lapply(seq_len(nrow(specs)), function(i) {
    gsm <- specs$gsm[i]
    fn <- specs$file[i]
    n <- as.integer(sub("^GSM", "", gsm))
    url <- sprintf(
      "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM%dnnn/%s/suppl/%s",
      n %/% 1000L, gsm, fn
    )
    dest <- file.path(dest_dir, fn)
    if (!file.exists(dest) || file.info(dest)$size < 100) {
      message("Downloading ", url)
      utils::download.file(url, dest, mode = "wb", quiet = TRUE)
    }
    d <- read_tsv(dest, col_names = c("gene_id", "count"), col_types = "cd",
                  show_col_types = FALSE)
    # Keep symbols intact. strip() is for Ensembl versions and would merge
    # Epb4.1 with Epb4.2 (and the H2-M10 family) before the join.
    d$gene_id <- as.character(d$gene_id)
    agg <- tapply(d$count, d$gene_id, sum)
    tibble(gene_id = names(agg), !!specs$sample[i] := as.numeric(agg))
  })
  mat_df <- Reduce(function(a, b) full_join(a, b, by = "gene_id"), cols)
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat <- round(mat)
  rownames(mat) <- mat_df$gene_id
  group <- specs$group[match(colnames(mat), specs$sample)]
  if (sum(group == "uv") != 3L || sum(group == "control") != 3L) {
    stop("UVC groups: ", paste(colnames(mat), group, sep = "=", collapse = ", "))
  }
  res <- run_deseq_pair(mat, group, "uv")
  map <- map_mouse_ids_to_ensg(res$gene_raw)
  write_hum(
    res, map, "GSE66286_NIH3T3_UVC_6h_vs_ctrl",
    "dna_damage", "UVC_80Jm2_6h_vs_ctrl", "DESeq2"
  )
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

# -----------------------------------------------------------------------------
# WT tendon fibroblasts, 1% vs 20% O2. Drop HIF1α cKO and Cre-gesicle arms.
# Count columns are library IDs; titles from the GEO series matrix:
#   Ctrl_20_01/02/03 = 13107-SS-1 / SS-2 / SS-3
#   Ctrl_1_01/02/03  = 13107-SS-10 / SS-11 / SS-12
run_one("GSE328685 WT tendon fibroblast hypoxia (counts)", {
  f <- download_geo_suppl(
    "GSE328685", "GSE328685_03_inVitro_RNA-seq_GeneCounts.txt.gz", RAW
  )
  d <- read_tsv(f, show_col_types = FALSE)
  gene_col <- if ("gene_id" %in% names(d)) "gene_id" else names(d)[1]
  nor <- c("13107-SS-1", "13107-SS-2", "13107-SS-3")
  hyp <- c("13107-SS-10", "13107-SS-11", "13107-SS-12")
  miss <- setdiff(c(nor, hyp), names(d))
  if (length(miss)) stop("missing count columns: ", paste(miss, collapse = ", "))
  mat <- as.matrix(d[, c(nor, hyp), drop = FALSE])
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat <- round(mat)
  rownames(mat) <- strip(d[[gene_col]])
  group <- ifelse(colnames(mat) %in% hyp, "hypoxia", "control")
  res <- run_deseq_pair(mat, group, "hypoxia")
  map <- map_mouse_ids_to_ensg(res$gene_raw)
  write_hum(
    res, map, "GSE328685_tendon_fibroblast_hypoxia_vs_normoxia",
    "hypoxia", "WT_1pct_vs_20pct", "DESeq2"
  )
})

# -----------------------------------------------------------------------------
# Immortalized mouse fibroblasts, WT 1% vs 21% O2 24 h. Drop mir210 KO.
# Author-normalized integer table; DESeq2 on rounded values, method flagged
# so 02 inflates SE 1.5× (n=2 and not raw counts).
run_one("GSE282023 WT immortalized fibroblast hypoxia (normalized counts)", {
  f <- download_geo_suppl(
    "GSE282023", "GSE282023_NormalizedTableCount_RNA.txt.gz", RAW
  )
  d <- read_tsv(f, show_col_types = FALSE)
  gene_col <- names(d)[1]
  nor <- grep("^WT_NX_", names(d), value = TRUE)
  hyp <- grep("^WT_HX_", names(d), value = TRUE)
  if (length(nor) < 2L || length(hyp) < 2L) {
    stop("WT NX/HX columns: nor=", length(nor), " hyp=", length(hyp))
  }
  mat <- as.matrix(d[, c(nor, hyp), drop = FALSE])
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat <- round(pmax(mat, 0))
  rownames(mat) <- strip(d[[gene_col]])
  group <- ifelse(colnames(mat) %in% hyp, "hypoxia", "control")
  res <- run_deseq_pair(mat, group, "hypoxia")
  map <- map_mouse_ids_to_ensg(res$gene_raw)
  write_hum(
    res, map, "GSE282023_immortalized_fibroblast_hypoxia_vs_normoxia",
    "hypoxia", "WT_1pct_vs_21pct_24h", "DESeq2_normalized"
  )
})

# -----------------------------------------------------------------------------
# WT MEFs, 1% vs 20% O2. Drop the tafazzin-knockout arm (TAZH / TAZN).
# Per-sample geneCounts live in GSE119775_RAW.tar (no header; ENSMUSG, count).
run_one("GSE119775 WT MEF hypoxia (counts)", {
  tarf <- download_geo_suppl("GSE119775", "GSE119775_RAW.tar", RAW)
  ex <- file.path(RAW, "GSE119775_raw")
  if (!dir.exists(ex)) {
    dir.create(ex, recursive = TRUE)
    utils::untar(tarf, exdir = ex)
  }
  files <- list.files(ex, pattern = "geneCounts\\.txt\\.gz$", full.names = TRUE, recursive = TRUE)
  wt <- files[grepl("_WT[HN][123]_", basename(files))]
  if (length(wt) != 6L) {
    stop("expected 6 WT count files, found ", length(wt), ": ",
         paste(basename(files), collapse = ", "))
  }
  cols <- lapply(wt, function(f) {
    d <- read_tsv(f, col_names = c("gene_id", "count"), col_types = "cd",
                  show_col_types = FALSE)
    nm <- sub(".*_(WT[HN][123])_.*", "\\1", basename(f))
    tibble(gene_id = strip(d$gene_id), !!nm := d$count)
  })
  mat_df <- Reduce(function(a, b) full_join(a, b, by = "gene_id"), cols)
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  mat <- round(mat)
  rownames(mat) <- mat_df$gene_id
  group <- ifelse(grepl("^WTH", colnames(mat)), "hypoxia", "control")
  if (sum(group == "hypoxia") != 3L || sum(group == "control") != 3L) {
    stop("WT groups: ", paste(colnames(mat), group, sep = "=", collapse = ", "))
  }
  res <- run_deseq_pair(mat, group, "hypoxia")
  map <- map_mouse_ids_to_ensg(res$gene_raw)
  write_hum(
    res, map, "GSE119775_MEF_WT_hypoxia_vs_normoxia",
    "hypoxia", "WT_1pct_vs_20pct", "DESeq2"
  )
})

cat("\nDone. New mouse files:\n")
print(list.files(OUT, pattern = "GSE(215293|243284|197536|281390|66286|114086|171871|205381|328685|282023|119775)"))
