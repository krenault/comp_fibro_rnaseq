#!/usr/bin/env Rscript
# Gene-level Spearman heatmap.
# Fibroblast longevity differences at each condition, Ma et al. slopes,
# and Tyshkovskiy Signatures 4.7. Style matches the BBM external
# correlation heatmaps (corrplot, blue-white-red, coefficient in the cell,
# BH stars on the correlation).

suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
})

root <- "/Users/katiarenault/Documents/GitHub/comp_fibro_rnaseq"
source(file.path(root, "generic_stress/scripts/shared/plot_aesthetics.R"))

plot_dir <- file.path(root, "longevity/results/static_genes/plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

static_files <- c(
  "32C" = "longevity/data/longevity_genes_static_32C.csv",
  "37C" = "longevity/data/longevity_genes_static_37C.csv",
  "41C" = "longevity/data/longevity_genes_static_41C.csv",
  "Hypoxia 6H" = "longevity/data/longevity_genes_static_hypoxia_6h.csv",
  "Hypoxia 24H" = "longevity/data/longevity_genes_static_hypoxia_24h.csv"
)

fibro <- imap_dfr(static_files, function(path, lab) {
  read_csv(file.path(root, path), show_col_types = FALSE) %>%
    filter(passes_all_filters %in% TRUE, is.finite(longevity_log2diff), hgnc_symbol != "") %>%
    transmute(symbol = toupper(hgnc_symbol), dataset = lab, effect = longevity_log2diff) %>%
    group_by(symbol, dataset) %>%
    slice_max(abs(effect), n = 1, with_ties = FALSE) %>%
    ungroup()
})

ma <- read_csv(file.path(root, "longevity/results/static_genes/ma_allgene_slopes_vs_37C.csv"),
               show_col_types = FALSE) %>%
  transmute(
    symbol = toupper(symbol),
    `Ma ML` = ma_slope_logML,
    `Ma MLres` = ma_slope_MLres
  ) %>%
  filter(symbol != "") %>%
  group_by(symbol) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_longer(-symbol, names_to = "dataset", values_to = "effect") %>%
  filter(is.finite(effect))

sig <- readRDS("/Users/katiarenault/PhD/Databases/Signatures_4.7.rds")
sig_long <- imap_dfr(sig, function(d, nm) {
  m <- regmatches(nm, regexec("^Species: (Brain|Kidney|Liver) (ML|MLres)$", nm))[[1]]
  if (length(m) < 3) return(NULL)
  tissue <- m[2]
  trait <- if (m[3] == "ML") "ML" else "MLres"
  tibble(
    symbol = toupper(trimws(as.character(d$Mouse_Gene))),
    dataset = paste(tissue, trait),
    effect = as.numeric(d$Slope)
  ) %>%
    filter(symbol != "", is.finite(effect)) %>%
    group_by(symbol, dataset) %>%
    slice(1) %>%
    ungroup()
})

effects <- bind_rows(fibro, ma, sig_long)
order_cols <- c(
  "32C", "37C", "41C", "Hypoxia 6H", "Hypoxia 24H",
  "Ma ML", "Ma MLres",
  "Brain ML", "Brain MLres", "Kidney ML", "Kidney MLres", "Liver ML", "Liver MLres"
)
wide <- effects %>%
  filter(dataset %in% order_cols) %>%
  dplyr::select(symbol, dataset, effect) %>%
  pivot_wider(names_from = dataset, values_from = effect) %>%
  dplyr::select(symbol, all_of(order_cols))

mat <- as.matrix(wide[, -1])
rownames(mat) <- wide$symbol

n <- ncol(mat)
cm <- matrix(NA_real_, n, n, dimnames = list(order_cols, order_cols))
pv <- cm
nn <- cm
diag(cm) <- 1
diag(pv) <- 0
diag(nn) <- nrow(mat)

for (i in seq_len(n - 1)) {
  for (j in (i + 1):n) {
    ok <- is.finite(mat[, i]) & is.finite(mat[, j])
    if (sum(ok) < 20) next
    tt <- cor.test(mat[ok, i], mat[ok, j], method = "spearman", exact = FALSE)
    cm[i, j] <- cm[j, i] <- unname(tt$estimate)
    pv[i, j] <- pv[j, i] <- tt$p.value
    nn[i, j] <- nn[j, i] <- sum(ok)
  }
}
fdr <- pv
fdr[upper.tri(fdr)] <- p.adjust(pv[upper.tri(pv)], method = "BH")
fdr[lower.tri(fdr)] <- t(fdr)[lower.tri(fdr)]

write_csv(as.data.frame(cm) %>% mutate(dataset = order_cols, .before = 1),
          file.path(plot_dir, "corr_fibro_ma_sig47.csv"))
write_csv(as.data.frame(fdr) %>% mutate(dataset = order_cols, .before = 1),
          file.path(plot_dir, "corr_fibro_ma_sig47_fdr.csv"))
write_csv(as.data.frame(nn) %>% mutate(dataset = order_cols, .before = 1),
          file.path(plot_dir, "corr_fibro_ma_sig47_n.csv"))

label_col <- c(
  "32C" = unname(treatment_colors["cold"]),
  "37C" = "#3F545E",
  "41C" = unname(treatment_colors["heat"]),
  "Hypoxia 6H" = unname(treatment_colors["hypoxia"]),
  "Hypoxia 24H" = "#67A9CF",
  "Ma ML" = "#7CB577",
  "Ma MLres" = "#3D6B3A",
  "Brain ML" = "#31a6ad",
  "Brain MLres" = "#31a6ad",
  "Kidney ML" = "#D99F6A",
  "Kidney MLres" = "#D99F6A",
  "Liver ML" = "#9b383a",
  "Liver MLres" = "#9b383a"
)

draw_corrplot <- function() {
  corrplot(
    cm, method = "color", type = "full",
    col = colorRampPalette(c(gradient_low, "white", gradient_high))(200),
    tl.col = unname(label_col[colnames(cm)]),
    tl.srt = 45, tl.cex = 1.05,
    addCoef.col = "black", number.cex = 0.72,
    mar = c(0, 0, 5.8, 0),
    addgrid.col = "white", is.corr = TRUE,
    title = paste0(
      "Culture levels plus raw and residual lifespan\n",
      sprintf(
        "Ma correlations use the genes present in both columns. 37C versus Ma MLres, n = %s.",
        format(nn["37C", "Ma MLres"], big.mark = ",")
      )
    )
  )
  n_cols <- ncol(cm)
  for (i in seq_len(n_cols)) {
    for (j in seq_len(n_cols)) {
      if (i == j) next
      fv <- fdr[i, j]
      if (!is.finite(fv)) next
      star <- if (fv < 0.001) "***" else if (fv < 0.01) "**" else if (fv < 0.05) "*" else ""
      if (nzchar(star)) text(j, n_cols - i + 1 - 0.32, star, cex = 0.62, col = "black")
    }
  }
  ma_col <- grepl("^Ma", colnames(cm))
  for (j in which(ma_col)) {
    for (i in seq_len(n_cols)) {
      if (i == j || !is.finite(nn[i, j])) next
      text(j, n_cols - i + 1 + 0.28, nn[i, j], cex = 0.48, col = "grey30")
    }
  }
}

png(file.path(plot_dir, "corr_fibro_ma_sig47.png"),
    width = 13.5, height = 12.8, units = "in", res = PUB_DPI)
draw_corrplot()
dev.off()
cairo_pdf(file.path(plot_dir, "corr_fibro_ma_sig47.pdf"), width = 13.5, height = 12.8)
draw_corrplot()
dev.off()

message("Pairwise n range: ", min(nn[upper.tri(nn)]), " to ", max(nn[upper.tri(nn)]))
ext <- c("Ma ML", "Ma MLres", "Brain ML", "Brain MLres", "Kidney ML", "Kidney MLres", "Liver ML", "Liver MLres")
ours <- c("32C", "37C", "41C", "Hypoxia 6H", "Hypoxia 24H")
message("Fibro vs external, 37C row:")
print(round(cm["37C", ext], 3))
message("Wrote corr_fibro_ma_sig47.png")
