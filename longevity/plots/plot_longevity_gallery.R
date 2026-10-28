#!/usr/bin/env Rscript
# Summary figures for the fibroblast longevity analyses.
# Colors, tiles, and correlation heatmaps follow the BBM figure style.

suppressPackageStartupMessages({
  library(tidyverse)
  library(corrplot)
  library(patchwork)
})

root <- "/Users/katiarenault/Documents/GitHub/comp_fibro_rnaseq"
source(file.path(root, "generic_stress/scripts/shared/plot_aesthetics.R"))
plot_root <- file.path(root, "longevity/plots")
fig <- file.path(plot_root, "supplement")
talk <- file.path(plot_root, "talk")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
dir.create(talk, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, file, w, h) {
  ggsave(file.path(fig, paste0(file, ".png")), p, width = w, height = h, dpi = PUB_DPI, bg = "white")
  ggsave(file.path(fig, paste0(file, ".pdf")), p, width = w, height = h, bg = "white")
  message("Wrote ", file)
}

static_conds <- c("32C", "37C", "41C", "hypoxia_6h", "hypoxia_24h")
static_labs <- c("32C" = "32C", "37C" = "37C", "41C" = "41C",
                 "hypoxia_6h" = "Hypoxia 6H", "hypoxia_24h" = "Hypoxia 24H")
dyn_conds <- c("32C_vs_37C", "41C_vs_37C", "hypoxia_6h_vs_0h", "hypoxia_24h_vs_0h")
dyn_labs <- c(
  "32C_vs_37C" = "Cold 32C", "41C_vs_37C" = "Heat 41C",
  "hypoxia_6h_vs_0h" = "Hypoxia 6H", "hypoxia_24h_vs_0h" = "Hypoxia 24H"
)
cond_colors <- c(
  "32C" = unname(treatment_colors["cold"]),
  "37C" = "#3F545E",
  "41C" = unname(treatment_colors["heat"]),
  "Hypoxia 6H" = unname(treatment_colors["hypoxia"]),
  "Hypoxia 24H" = "#67A9CF",
  "Cold 32C" = unname(treatment_colors["cold"]),
  "Heat 41C" = unname(treatment_colors["heat"])
)
dir_colors <- c(
  higher = unname(treatment_colors["heat"]),
  lower = unname(treatment_colors["hypoxia"])
)

read_static <- function(cond) {
  read_csv(file.path(root, "longevity/data", paste0("longevity_genes_static_", cond, ".csv")),
           show_col_types = FALSE) %>%
    transmute(
      gene, symbol = toupper(hgnc_symbol),
      diff = longevity_log2diff,
      called = effecttype_fc1_1 %in% c("higher_in_long_lived", "lower_in_long_lived"),
      direction = case_when(
        effecttype_fc1_1 == "higher_in_long_lived" ~ "higher",
        effecttype_fc1_1 == "lower_in_long_lived" ~ "lower",
        TRUE ~ NA_character_
      ),
      condition = cond
    ) %>%
    filter(symbol != "", is.finite(diff))
}
static <- map_dfr(static_conds, read_static) %>%
  mutate(condition_lab = factor(static_labs[condition], levels = unname(static_labs)))

read_dyn <- function(cond) {
  read_csv(file.path(root, "longevity/results/dynamic_genes",
                     paste0("longevity_genes_dynamic_", cond, ".csv")),
           show_col_types = FALSE) %>%
    transmute(
      gene, symbol = toupper(hgnc_symbol),
      diff = longevity_logFC_difference,
      effect = effecttype_fc1_1,
      called = !is.na(effecttype_fc1_1) & effecttype_fc1_1 != "",
      condition = cond
    ) %>%
    filter(symbol != "", is.finite(diff))
}
dynamic <- map_dfr(dyn_conds, read_dyn) %>%
  mutate(condition_lab = factor(dyn_labs[condition], levels = unname(dyn_labs)))

theme_tile <- function(base = 13) {
  theme_pub(base) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      legend.position = "right"
    )
}

# --- 01 and 17 correlation heatmaps ------------------------------------------
pair_corr <- function(wide, cols) {
  m <- as.matrix(wide[, cols])
  n <- length(cols)
  cm <- matrix(NA_real_, n, n, dimnames = list(cols, cols))
  pv <- cm
  nn <- cm
  diag(cm) <- 1
  diag(pv) <- 0
  diag(nn) <- colSums(is.finite(m))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      ok <- is.finite(m[, i]) & is.finite(m[, j])
      nn[i, j] <- nn[j, i] <- sum(ok)
      if (sum(ok) < 20) next
      tt <- cor.test(m[ok, i], m[ok, j], method = "spearman", exact = FALSE)
      cm[i, j] <- cm[j, i] <- unname(tt$estimate)
      pv[i, j] <- pv[j, i] <- tt$p.value
    }
  }
  fdr <- pv
  fdr[upper.tri(fdr)] <- p.adjust(pv[upper.tri(pv)], method = "BH")
  fdr[lower.tri(fdr)] <- t(fdr)[lower.tri(fdr)]
  list(cm = cm, fdr = fdr, nn = nn)
}

save_corrplot <- function(cm, fdr, colors, title, file, w, h, num_cex, mar_top = 3.6, nn = NULL) {
  draw <- function() {
    corrplot(
      cm, method = "color", type = "full",
      col = colorRampPalette(c(gradient_low, "white", gradient_high))(200),
      tl.col = unname(colors[colnames(cm)]), tl.srt = 45, tl.cex = 1.45,
      addCoef.col = "black", number.cex = num_cex,
      mar = c(0, 0, mar_top, 0), addgrid.col = "white", is.corr = TRUE,
      title = title
    )
    n_cols <- ncol(cm)
    for (i in seq_len(n_cols)) {
      for (j in seq_len(n_cols)) {
        if (i == j) next
        fv <- fdr[i, j]
        if (!is.finite(fv) || fv >= 0.05) next
        star <- if (fv < 0.001) "***" else if (fv < 0.01) "**" else "*"
        text(j, n_cols - i + 1 - 0.28, star, cex = 1.15, col = "black")
      }
    }
    if (!is.null(nn)) {
      ma_col <- grepl("^Ma", colnames(cm))
      for (j in which(ma_col)) {
        for (i in seq_len(n_cols)) {
          if (i == j || !is.finite(nn[i, j])) next
          text(j, n_cols - i + 1 + 0.30, nn[i, j], cex = num_cex * 0.62, col = "grey30")
        }
      }
    }
  }
  png(file.path(fig, paste0(file, ".png")), width = w, height = h, units = "in", res = PUB_DPI)
  draw(); dev.off()
  cairo_pdf(file.path(fig, paste0(file, ".pdf")), width = w, height = h)
  draw(); dev.off()
  message("Wrote ", file)
}

static_wide <- static %>%
  dplyr::select(symbol, condition_lab, diff) %>%
  pivot_wider(names_from = condition_lab, values_from = diff) %>%
  group_by(symbol) %>% slice(1) %>% ungroup()
static_cols <- unname(static_labs)
sc <- pair_corr(static_wide, static_cols)
save_corrplot(
  sc$cm, sc$fdr, cond_colors,
  "Static level rankings compared with each other\nSpearman correlation",
  "static_spearman", 8.2, 7.4, 1.4
)

dyn_wide <- dynamic %>%
  dplyr::select(symbol, condition_lab, diff) %>%
  pivot_wider(names_from = condition_lab, values_from = diff) %>%
  group_by(symbol) %>% slice(1) %>% ungroup()
dyn_cols <- unname(dyn_labs)
dc <- pair_corr(dyn_wide, dyn_cols)
save_corrplot(
  dc$cm, dc$fdr, cond_colors,
  "Response rankings compared with each other\nSpearman correlation",
  "dynamic_spearman", 7.6, 6.8, 1.4
)

# --- one matrix: levels, responses, Ma, tissues ------------------------------
read_effect <- function(path, value, lab) {
  read_csv(file.path(root, path), show_col_types = FALSE) %>%
    filter(passes_all_filters %in% TRUE, is.finite(.data[[value]]), hgnc_symbol != "") %>%
    transmute(symbol = toupper(hgnc_symbol), dataset = lab, effect = .data[[value]]) %>%
    group_by(symbol) %>%
    slice_max(abs(effect), n = 1, with_ties = FALSE) %>%
    ungroup()
}
level_files <- c(
  "32C" = "longevity/data/longevity_genes_static_32C.csv",
  "37C" = "longevity/data/longevity_genes_static_37C.csv",
  "41C" = "longevity/data/longevity_genes_static_41C.csv",
  "6H level" = "longevity/data/longevity_genes_static_hypoxia_6h.csv",
  "24H level" = "longevity/data/longevity_genes_static_hypoxia_24h.csv"
)
resp_files <- c(
  "Cold response" = "longevity/results/dynamic_genes/longevity_genes_dynamic_32C_vs_37C.csv",
  "Heat response" = "longevity/results/dynamic_genes/longevity_genes_dynamic_41C_vs_37C.csv",
  "6H response" = "longevity/results/dynamic_genes/longevity_genes_dynamic_hypoxia_6h_vs_0h.csv",
  "24H response" = "longevity/results/dynamic_genes/longevity_genes_dynamic_hypoxia_24h_vs_0h.csv"
)
fibro_cmp <- bind_rows(
  imap_dfr(level_files, function(path, lab) read_effect(path, "longevity_log2diff", lab)),
  imap_dfr(resp_files, function(path, lab) read_effect(path, "longevity_logFC_difference", lab))
)
ma_cmp <- read_csv(file.path(root, "longevity/results/static_genes/ma_allgene_slopes_vs_37C.csv"),
                   show_col_types = FALSE) %>%
  transmute(symbol = toupper(symbol), dataset = "Ma MLres", effect = ma_slope_MLres) %>%
  filter(symbol != "", is.finite(effect)) %>%
  group_by(symbol) %>% slice(1) %>% ungroup()
sig_cmp <- readRDS("/Users/katiarenault/PhD/Databases/Signatures_4.7.rds")
sig_cmp <- imap_dfr(sig_cmp, function(d, nm) {
  m <- regmatches(nm, regexec("^Species: (Brain|Kidney|Liver) MLres$", nm))[[1]]
  if (length(m) < 2) return(NULL)
  tibble(symbol = toupper(trimws(as.character(d$Mouse_Gene))),
         dataset = paste(m[2], "MLres"),
         effect = as.numeric(d$Slope)) %>%
    filter(symbol != "", is.finite(effect)) %>%
    group_by(symbol) %>% slice(1) %>% ungroup()
})
cmp_order <- c(
  "32C", "37C", "41C", "6H level", "24H level",
  "Cold response", "Heat response", "6H response", "24H response",
  "Ma MLres", "Brain MLres", "Kidney MLres", "Liver MLres"
)
cmp_wide <- bind_rows(fibro_cmp, ma_cmp, sig_cmp) %>%
  filter(dataset %in% cmp_order) %>%
  pivot_wider(names_from = dataset, values_from = effect) %>%
  dplyr::select(all_of(cmp_order))
cmp <- pair_corr(cmp_wide, cmp_order)
cmp_colors <- c(
  "32C" = unname(treatment_colors["cold"]),
  "37C" = "#3F545E",
  "41C" = unname(treatment_colors["heat"]),
  "6H level" = unname(treatment_colors["hypoxia"]),
  "24H level" = "#67A9CF",
  "Cold response" = unname(treatment_colors["cold"]),
  "Heat response" = unname(treatment_colors["heat"]),
  "6H response" = unname(treatment_colors["hypoxia"]),
  "24H response" = "#67A9CF",
  "Ma MLres" = "#3D6B3A",
  "Brain MLres" = "#31a6ad",
  "Kidney MLres" = "#D99F6A",
  "Liver MLres" = "#9b383a"
)
n_ma <- cmp$nn["37C", "Ma MLres"]
save_corrplot(
  cmp$cm, cmp$fdr, cmp_colors,
  paste0(
    "Culture levels plus perturbation responses\n",
    sprintf("Ma correlations use the genes present in both columns. 37C versus Ma MLres, n = %s.",
            format(n_ma, big.mark = ","))
  ),
  "levels_plus_responses", 13.6, 12.8, 0.9, mar_top = 4.6,
  nn = cmp$nn
)
ma_pts <- cmp_wide %>%
  transmute(fibro = `37C`, ma = `Ma MLres`) %>%
  filter(is.finite(fibro), is.finite(ma))
ma_rho <- cor.test(ma_pts$fibro, ma_pts$ma, method = "spearman", exact = FALSE)
p_ma <- ggplot(ma_pts, aes(ma, fibro)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_point(color = "#3D6B3A", size = 1.35, alpha = 0.35) +
  geom_smooth(method = "lm", se = FALSE, color = "grey20", linewidth = 0.7) +
  labs(
    title = "37C fibroblast difference versus Ma residual lifespan",
    subtitle = sprintf(
      "Spearman rho = %.2f, n = %s genes present in both.",
      unname(ma_rho$estimate), format(nrow(ma_pts), big.mark = ",")
    ),
    x = "Ma MLres slope",
    y = "Fibroblast longevity difference at 37C",
    caption = "Genes absent from the 37C filter pass or from the Ma slope file are omitted."
  ) +
  theme_pub(15)
save_fig(p_ma, "ma_shared_scatter", 7.4, 6.6)
message("37C vs responses: ", paste(round(cmp$cm["37C", c("Cold response", "Heat response", "6H response", "24H response")], 3), collapse = ", "))
message("37C vs Ma/tissues: ", paste(round(cmp$cm["37C", c("Ma MLres", "Brain MLres", "Kidney MLres", "Liver MLres")], 3), collapse = ", "))
message("37C vs Ma n: ", n_ma)

# --- 03 call counts ----------------------------------------------------------
call_counts <- static %>%
  filter(called) %>%
  count(condition_lab, direction) %>%
  mutate(direction = factor(direction, levels = c("higher", "lower")))
p_counts <- ggplot(call_counts, aes(condition_lab, n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66,
           color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(width = 0.72),
            vjust = -0.3, size = 3.3) +
  scale_fill_manual(values = dir_colors, labels = c(higher = "Higher", lower = "Lower"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Longevity calls at each condition",
    subtitle = "fc1.1 calls. Higher means greater expression at the higher-longevity anchor.",
    x = NULL, y = "Genes"
  ) +
  theme_pub(16) + theme(legend.position = "bottom")
save_fig(p_counts, "call_counts", 8.2, 5.2)

# --- 04 persistence ----------------------------------------------------------
persist <- static %>%
  filter(called) %>%
  group_by(symbol) %>%
  summarise(
    n_cond = n_distinct(condition),
    direction = if (n_distinct(direction) == 1) first(direction) else "mixed",
    .groups = "drop"
  )
persist_n <- persist %>%
  count(n_cond, direction) %>%
  mutate(n_cond = factor(n_cond))
p_persist <- ggplot(persist_n, aes(n_cond, n, fill = direction)) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.4, color = "white") +
  scale_fill_manual(values = dir_colors, labels = c(higher = "Higher", lower = "Lower"), name = NULL) +
  labs(
    title = "How many conditions call the same gene",
    subtitle = "Every gene called in two or more conditions keeps the same sign.",
    x = "Conditions with an fc1.1 call", y = "Genes"
  ) +
  theme_pub(16) + theme(legend.position = "bottom")
save_fig(p_persist, "call_persistence", 7.4, 5.2)

# --- 05 unique calls vs the 37C sign ----------------------------------------
sign37 <- static %>% filter(condition == "37C") %>% dplyr::select(symbol, sign37 = diff)
ncall <- static %>% filter(called) %>% distinct(symbol, condition) %>% count(symbol, name = "n_cond")
unique_calls <- static %>%
  filter(called) %>%
  left_join(ncall, by = "symbol") %>%
  filter(n_cond == 1, condition != "37C") %>%
  left_join(sign37, by = "symbol") %>%
  mutate(
    rope = log2(1.1),
    relation = case_when(
      sign(diff) == sign(sign37) ~ "Same sign as 37C",
      abs(sign37) >= rope ~ "Opposite 37C",
      TRUE ~ "37C inside the margin"
    )
  )
uniq_n <- unique_calls %>% count(condition_lab, relation)
p_uniq <- ggplot(uniq_n, aes(condition_lab, n, fill = relation)) +
  geom_col(width = 0.68, color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.4, color = "grey10") +
  scale_fill_manual(
    values = c(
      "Same sign as 37C" = unname(treatment_colors["hypoxia"]),
      "37C inside the margin" = "grey80",
      "Opposite 37C" = unname(treatment_colors["heat"])
    ),
    name = NULL
  ) +
  labs(
    title = "Calls made in only one condition",
    subtitle = "Opposite means the 37C difference has the other sign and sits outside the fc1.1 margin.",
    x = NULL, y = "Genes called in that condition only"
  ) +
  theme_pub(16) + theme(legend.position = "bottom")
save_fig(p_uniq, "unique_calls_vs_37C", 8.0, 5.2)

# --- 06 07 pathway NES -------------------------------------------------------
static_paths <- c(
  "HALLMARK_E2F_TARGETS" = "E2F targets",
  "REACTOME_COLLAGEN_FORMATION" = "Collagen formation",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = "Oxidative phosphorylation",
  "GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING" = "Chaperone-mediated folding",
  "GOBP_BASE_EXCISION_REPAIR" = "Base-excision repair",
  "GOBP_ATP_BIOSYNTHETIC_PROCESS" = "ATP synthesis",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" = "Unfolded protein response"
)
dyn_paths <- c(
  "GOBP_RIBOSOME_BIOGENESIS" = "Ribosome biogenesis",
  "REACTOME_FATTY_ACID_METABOLISM" = "Fatty-acid metabolism",
  "REACTOME_TRANSLATION" = "Translation",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE" = "Unfolded protein response",
  "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT" = "ATP synthesis, electron transport",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT" = "Respiratory electron transport",
  "GOBP_OXIDATIVE_PHOSPHORYLATION" = "Oxidative phosphorylation"
)

nes_grid <- function(path, pathways, analyses, labs) {
  raw <- read_csv(path, show_col_types = FALSE) %>%
    filter(pathway %in% names(pathways), analysis %in% analyses) %>%
    dplyr::select(analysis, pathway, NES, padj)
  expand_grid(analysis = analyses, pathway = names(pathways)) %>%
    left_join(raw, by = c("analysis", "pathway")) %>%
    mutate(
      analysis_lab = factor(labs[analysis], levels = unname(labs)),
      pathway_lab = factor(unname(pathways[pathway]), levels = rev(unname(pathways))),
      star = get_sig_star(padj),
      label = if_else(is.finite(NES), paste0(sprintf("%+.1f", NES), star), "")
    )
}

heat_plot <- function(df, title, subtitle, caption) {
  lim <- max(abs(df$NES), na.rm = TRUE)
  ggplot(df, aes(analysis_lab, pathway_lab, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.7) +
    geom_text(aes(label = label, color = abs(NES) > 1.7), size = 4.6, show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey10"), na.value = "grey10") +
    scale_fill_gradient2(
      low = gradient_low, mid = "white", high = gradient_high,
      midpoint = 0, limits = c(-lim, lim), na.value = "grey92", name = "NES"
    ) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL, caption = caption) +
    theme_tile(16)
}

p_nes_s <- heat_plot(
  nes_grid(file.path(root, "longevity/results/static_genes/fgsea_msigdb_v2023_2_all.csv"),
           static_paths, static_conds, static_labs),
  "Longevity-associated expression level",
  "Positive NES: higher at the higher-longevity anchor.",
  paste0(sig_caption, ". Each column is one culture condition.")
)
p_nes_d <- heat_plot(
  nes_grid(file.path(root, "longevity/results/dynamic_genes/fgsea_msigdb_v2023_2_all.csv"),
           dyn_paths, dyn_conds, dyn_labs),
  "Each perturbation rearranges a different program",
  "Positive NES: a more positive fold change at the higher-longevity anchor. Blank tiles are not enriched.",
  paste0(sig_caption, ". Each column is a perturbation relative to baseline.")
)
save_fig(p_nes_s, "static_pathways", 9.6, 5.8)
save_fig(p_nes_d, "dynamic_pathways", 9.2, 5.4)

# --- 07 response signatures --------------------------------------------------
sig_keep <- tibble::tribble(
  ~analysis, ~pathway, ~pathway_lab,
  "32C_vs_37C", "REACTOME_FATTY_ACID_METABOLISM", "Fatty-acid metabolism",
  "32C_vs_37C", "GOBP_RIBOSOME_BIOGENESIS", "Ribosome biogenesis",
  "41C_vs_37C", "REACTOME_TRANSLATION", "Translation",
  "41C_vs_37C", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "Unfolded protein response",
  "41C_vs_37C", "GOBP_OXIDATIVE_PHOSPHORYLATION", "Oxidative phosphorylation",
  "hypoxia_6h_vs_0h", "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT", "Electron transport",
  "hypoxia_6h_vs_0h", "GOBP_NCRNA_PROCESSING", "ncRNA processing",
  "hypoxia_24h_vs_0h", "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT", "Electron transport",
  "hypoxia_24h_vs_0h", "GOBP_RIBOSOME_BIOGENESIS", "Ribosome biogenesis"
)
sig_df <- read_csv(
  file.path(root, "longevity/results/dynamic_genes/fgsea_msigdb_v2023_2_padj05.csv"),
  show_col_types = FALSE
) %>%
  inner_join(sig_keep, by = c("analysis", "pathway")) %>%
  mutate(
    condition_lab = factor(dyn_labs[analysis], levels = unname(dyn_labs)),
    star = get_sig_star(padj),
    label = paste0(sprintf("%+.1f", NES), star)
  ) %>%
  group_by(condition_lab) %>%
  mutate(pathway_lab = forcats::fct_reorder(pathway_lab, NES)) %>%
  ungroup()
p_sig <- ggplot(sig_df, aes(NES, pathway_lab, fill = NES)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.3) +
  geom_text(
    aes(
      x = NES + if_else(NES > 0, 0.12, -0.12),
      label = label,
      hjust = if_else(NES > 0, 0, 1)
    ),
    size = 3.1, color = "grey10"
  ) +
  scale_fill_gradient2(
    low = gradient_low, mid = "white", high = gradient_high,
    midpoint = 0, limits = c(-2.8, 2.8), oob = scales::squish, guide = "none"
  ) +
  scale_x_continuous(limits = c(-5.6, 5.2), breaks = seq(-2, 2, 2)) +
  facet_wrap(~condition_lab, nrow = 1, scales = "free_y") +
  labs(
    title = "Each perturbation has its own longevity-associated response",
    subtitle = "Positive NES: a more positive fold change at the higher-longevity anchor.",
    x = "Normalized enrichment score", y = NULL,
    caption = paste0(
      sig_caption, ". Electron transport is positive at 6 h and negative at 24 h. Full-fit estimates."
    )
  ) +
  theme_pub(15) +
  theme(panel.spacing.x = unit(1.1, "lines"))
save_fig(p_sig, "response_signatures", 13.6, 4.8)

# --- every dynamic FDR hit, no hand-picked list -----------------------------
pretty_pathway <- function(x) {
  raw <- dplyr::case_when(
    startsWith(x, "KEGG_MEDICUS_REFERENCE_") ~
      paste0("Reference: ", tolower(gsub("_", " ", sub("^KEGG_MEDICUS_REFERENCE_", "", x)))),
    startsWith(x, "KEGG_MEDICUS_VARIANT_") ~
      paste0("Variant: ", tolower(gsub("_", " ", sub(
        "^KEGG_MEDICUS_VARIANT_(MUTATION_CAUSED_ABERRANT_|MUTATION_INACTIVATED_|SCRAPIE_CONFORMATION_)",
        "", x
      )))),
    TRUE ~ stringr::str_to_sentence(tolower(gsub("_", " ", sub("^(HALLMARK|REACTOME|GOBP)_", "", x))))
  )
  raw <- gsub("ncrna", "ncRNA", raw, ignore.case = TRUE)
  raw <- gsub("rrna", "rRNA", raw, ignore.case = TRUE)
  raw <- gsub("(?<![A-Za-z])rna(?![A-Za-z])", "RNA", raw, ignore.case = TRUE, perl = TRUE)
  raw <- gsub("(?<![A-Za-z])dna(?![A-Za-z])", "DNA", raw, ignore.case = TRUE, perl = TRUE)
  raw <- gsub("atp", "ATP", raw, ignore.case = TRUE)
  raw <- gsub("nadh", "NADH", raw, ignore.case = TRUE)
  raw <- gsub("tca", "TCA", raw, ignore.case = TRUE)
  raw <- gsub("G2m", "G2/M", raw, fixed = TRUE)
  raw <- gsub("G2 m", "G2/M", raw, fixed = TRUE)
  raw <- gsub("g2 m", "G2/M", raw, fixed = TRUE)
  raw <- gsub("G1 s", "G1/S", raw, fixed = TRUE)
  raw <- gsub("Upr", "UPR", raw, fixed = TRUE)
  raw <- gsub("E2f", "E2F", raw, fixed = TRUE)
  for (a in list(
    c("snca", "SNCA"), c("tdp43", "TDP43"), c("pink1", "PINK1"),
    c("htt", "HTT"), c("ubqln2", "UBQLN2"), c("vcp", "VCP"),
    c("sod1", "SOD1"), c("prpsc", "PrPSc"), c("complex i", "complex I")
  )) raw <- gsub(a[[1]], a[[2]], raw, ignore.case = TRUE)
  raw
}
dyn_sig <- read_csv(
  file.path(root, "longevity/results/dynamic_genes/fgsea_msigdb_v2023_2_padj05.csv"),
  show_col_types = FALSE
) %>%
  filter(!startsWith(pathway, "KEGG_MEDICUS_")) %>%
  mutate(
    condition_lab = factor(dyn_labs[analysis], levels = unname(dyn_labs)),
    pathway_lab = pretty_pathway(pathway),
    label = paste0(sprintf("%+.1f", NES), get_sig_star(padj))
  )
share_class <- dyn_sig %>%
  group_by(pathway) %>%
  summarise(n_pos = sum(NES > 0), n_neg = sum(NES < 0), .groups = "drop") %>%
  mutate(share = case_when(
    n_pos + n_neg == 1 ~ "Unique",
    n_pos > 0 & n_neg > 0 ~ "Opposite direction",
    TRUE ~ "Same direction"
  ))
dyn_sig <- dyn_sig %>%
  left_join(share_class, by = "pathway") %>%
  mutate(share = factor(share, levels = c("Unique", "Same direction", "Opposite direction")))
message("Dynamic FDR sharing: ", paste(names(table(share_class$share)), table(share_class$share), sep = " ", collapse = "; "))
sig_panel <- function(df, title, show_legend = TRUE) {
  df <- df %>% mutate(pathway_lab = forcats::fct_reorder(pathway_lab, NES))
  ggplot(df, aes(NES, pathway_lab, fill = NES)) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70") +
    geom_col(aes(color = share), width = 0.72, linewidth = 1.15) +
    geom_text(aes(x = 2.55, label = label), hjust = 0, size = 4.2, color = "grey10") +
    scale_fill_gradient2(
      low = gradient_low, mid = "white", high = gradient_high,
      midpoint = 0, limits = c(-2.8, 2.8), oob = scales::squish, guide = "none"
    ) +
    scale_color_manual(
      values = c(
        "Unique" = "grey55",
        "Same direction" = "#1B7837",
        "Opposite direction" = "#762A83"
      ),
      limits = c("Unique", "Same direction", "Opposite direction"),
      labels = c("Unique", "Same direction", "Opposite direction"),
      drop = FALSE,
      name = "Across perturbations"
    ) +
    guides(color = guide_legend(
      override.aes = list(fill = "white", linewidth = 1.8),
      nrow = 1, byrow = TRUE
    )) +
    scale_x_continuous(limits = c(-3.2, 4.1), breaks = seq(-2, 2, 2)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_pub(15) +
    theme(
      axis.text.y = element_text(size = 12, color = "grey10"),
      axis.text.x = element_text(size = 12, color = "grey10"),
      plot.title = element_text(size = 15),
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 15, color = "grey10"),
      legend.title = element_text(size = 15, face = "bold", color = "grey10"),
      legend.key.width = unit(1.3, "cm"),
      legend.key.height = unit(0.55, "cm")
    )
}
h24 <- dyn_sig %>% filter(condition_lab == "Hypoxia 24H")
keep_24 <- bind_rows(
  slice_min(h24, NES, n = 5, with_ties = FALSE),
  slice_max(h24, NES, n = 5, with_ties = FALSE)
)
dyn_plot <- dyn_sig %>%
  filter(condition_lab != "Hypoxia 24H" | pathway %in% keep_24$pathway)
panel_titles <- c(
  "Cold 32C" = "Cold 32C",
  "Heat 41C" = "Heat 41C",
  "Hypoxia 6H" = "Hypoxia 6H",
  "Hypoxia 24H" = "Hypoxia 24H, strongest five each way"
)
p_panels <- lapply(unname(dyn_labs), function(lab) {
  sig_panel(filter(dyn_plot, condition_lab == lab), panel_titles[[lab]], show_legend = FALSE)
})
legend_df <- tibble(
  share = factor(
    c("Unique", "Same direction", "Opposite direction"),
    levels = c("Unique", "Same direction", "Opposite direction")
  )
)
p_legend <- ggplot(legend_df, aes(x = 1, y = 1, color = share)) +
  geom_point(shape = 22, size = 0, stroke = 0, fill = "white") +
  scale_color_manual(
    values = c(
      "Unique" = "grey55",
      "Same direction" = "#1B7837",
      "Opposite direction" = "#762A83"
    ),
    limits = c("Unique", "Same direction", "Opposite direction"),
    name = "Across perturbations"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 6, shape = 22, fill = "white", stroke = 1.5),
    nrow = 1, byrow = TRUE
  )) +
  coord_cartesian(ylim = c(2, 3), clip = "on") +
  theme_void(base_size = 15) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 15, color = "grey10"),
    legend.title = element_text(size = 15, face = "bold", color = "grey10"),
    legend.key.width = unit(1.3, "cm"),
    legend.key.height = unit(0.55, "cm")
  )
p_all <- wrap_plots(
  p_panels[[1]], p_panels[[2]], p_panels[[3]], p_panels[[4]], p_legend,
  design = "AB\nCD\nEE",
  heights = c(12, 10, 1.3)
) +
  plot_annotation(
    title = "Response pathways at FDR < 0.05",
    subtitle = "Cold, heat, and 6 h show every pathway. 24 h shows the five strongest in each direction.",
    caption = paste0(
      sig_caption,
      ". Bar outline: unique to one perturbation, shared with the same sign, or shared with the opposite sign. Sharing counts every 24 h pathway, including those not drawn. KEGG MEDICUS copies are omitted."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 20, hjust = 0.5, color = "grey10"),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30"),
      plot.caption = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )
save_fig(p_all, "all_sig_pathways", 18.5, 15.2)
p_all_24 <- sig_panel(filter(dyn_plot, condition_lab == "Hypoxia 24H"), NULL, show_legend = TRUE) +
  labs(
    title = "Strongest 24 h hypoxia pathways at FDR < 0.05",
    subtitle = "Five most negative and five most positive. KEGG MEDICUS copies are omitted.",
    x = "Normalized enrichment score",
    caption = paste0(
      sig_caption,
      ". Bar outline compares this pathway with cold, heat, and 6 h."
    )
  ) +
  theme(legend.position = "bottom")
save_fig(p_all_24, "all_sig_pathways_24h", 11.2, 7.2)

# --- static pathways with a stable direction ---------------------------------
st_all <- read_csv(
  file.path(root, "longevity/results/static_genes/fgsea_msigdb_v2023_2_all.csv"),
  show_col_types = FALSE
)
stable_ids <- st_all %>%
  group_by(pathway) %>%
  summarise(
    n_finite = sum(is.finite(NES)),
    n_pos = sum(NES > 0, na.rm = TRUE),
    n_neg = sum(NES < 0, na.rm = TRUE),
    n_sig = sum(padj < 0.05, na.rm = TRUE),
    med_NES = median(NES, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_finite >= 4, n_sig >= 2, n_pos == n_finite | n_neg == n_finite) %>%
  arrange(desc(med_NES))
message(
  "Stable pathways: ", nrow(stable_ids),
  ". Significant in all five conditions: ",
  st_all %>% group_by(pathway) %>% summarise(n5 = sum(padj < 0.05, na.rm = TRUE), .groups = "drop") %>%
    filter(n5 == 5) %>% nrow()
)
stable_paths <- setNames(pretty_pathway(stable_ids$pathway), stable_ids$pathway)
p_stable <- heat_plot(
  nes_grid(
    file.path(root, "longevity/results/static_genes/fgsea_msigdb_v2023_2_all.csv"),
    stable_paths, static_conds, static_labs
  ),
  "Pathways that keep their direction across conditions",
  "Same sign wherever a score exists, and FDR < 0.05 in at least two conditions.",
  paste0(sig_caption, ". Grey: no score. No pathway is significant in all five conditions.")
)
save_fig(p_stable, "consistent_pathways", 10.2, 5.6)

# --- 08 dynamic effect types -------------------------------------------------
type_group <- function(x) {
  case_when(
    x %in% c("increase_vs_decrease", "decrease_vs_increase") ~ "Opposite direction",
    x %in% c("greater_increase", "greater_decrease") ~ "Same direction, larger",
    x %in% c("smaller_increase", "smaller_decrease") ~ "Same direction, smaller",
    x %in% c("increase_vs_unchanged", "decrease_vs_unchanged",
             "unchanged_vs_increase", "unchanged_vs_decrease") ~ "One anchor moves",
    TRUE ~ NA_character_
  )
}
type_n <- dynamic %>%
  filter(called) %>%
  mutate(kind = type_group(effect)) %>%
  count(condition_lab, kind) %>%
  mutate(kind = factor(kind, levels = c(
    "Opposite direction", "One anchor moves",
    "Same direction, larger", "Same direction, smaller"
  )))
p_types <- ggplot(type_n, aes(condition_lab, n, fill = kind)) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.2, color = "white") +
  scale_fill_manual(
    values = c(
      "Opposite direction" = unname(treatment_colors["heat"]),
      "One anchor moves" = unname(treatment_colors["oxidative"]),
      "Same direction, larger" = unname(treatment_colors["hypoxia"]),
      "Same direction, smaller" = unname(treatment_colors["cold"])
    ),
    name = NULL, drop = FALSE
  ) +
  labs(
    title = "The called responses move in opposite directions",
    subtitle = "fc1.1 response calls. A smaller response of the same kind does not appear.",
    x = NULL, y = "Genes"
  ) +
  theme_pub(16) + theme(legend.position = "bottom")
save_fig(p_types, "dynamic_effect_types", 8.4, 5.6)

called37 <- static %>% filter(condition == "37C", called) %>% pull(symbol)
dyn_new <- dynamic %>%
  filter(called) %>%
  mutate(at_37 = symbol %in% called37) %>%
  count(condition_lab, at_37) %>%
  mutate(at_37 = if_else(at_37, "Also called at 37C", "Not called at 37C"))
p_new <- ggplot(dyn_new, aes(condition_lab, n, fill = at_37)) +
  geom_col(width = 0.68, color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.5, color = "white") +
  scale_fill_manual(
    values = c("Also called at 37C" = "#3F545E", "Not called at 37C" = unname(treatment_colors["heat"])),
    name = NULL
  ) +
  labs(
    title = "Response calls that baseline does not make",
    subtitle = "Dynamic fc1.1 calls split by whether 37C also calls the gene.",
    x = NULL, y = "Genes"
  ) +
  theme_pub(16) + theme(legend.position = "bottom")
save_fig(p_new, "dynamic_beyond_baseline", 8.0, 5.4)

# --- 10 gene-set medians -----------------------------------------------------
activation_genes <- c("COL1A1", "COL1A2", "POSTN", "ACTA2", "CCN2", "LOXL2", "ITGA11")
read_gmt_set <- function(set_name) {
  gmt <- "/Users/katiarenault/PhD/GSEA/c2.cgp.v2023.2.Hs.symbols.gmt"
  hit <- readLines(gmt)
  hit <- hit[startsWith(hit, paste0(set_name, "\t"))]
  toupper(strsplit(hit, "\t", fixed = TRUE)[[1]][-c(1, 2)])
}
set_defs <- list(
  Activation = activation_genes,
  "Basement membrane" = read_gmt_set("NABA_BASEMENT_MEMBRANES"),
  Collagens = read_gmt_set("NABA_COLLAGENS"),
  "Secreted factors" = read_gmt_set("NABA_SECRETED_FACTORS")
)
set_med <- imap_dfr(set_defs, function(genes, name) {
  static %>%
    group_by(condition_lab) %>%
    summarise(
      median = median(diff[symbol %in% genes], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(geneset = name)
}) %>%
  mutate(geneset = factor(geneset, levels = names(set_defs)))
p_sets <- ggplot(set_med, aes(condition_lab, median, fill = condition_lab)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%+.2f", median), vjust = ifelse(median < 0, 1.15, -0.35)), size = 3.0) +
  facet_wrap(~geneset, nrow = 1) +
  scale_fill_manual(values = cond_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.16))) +
  labs(
    title = "Longevity difference by gene set",
    subtitle = "Median log2 difference. Secreted factors stay near zero.",
    x = NULL, y = "Median longevity log2 difference"
  ) +
  theme_pub(15) + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_fig(p_sets, "geneset_medians", 11.2, 4.8)

act <- static %>%
  filter(symbol %in% activation_genes) %>%
  mutate(
    symbol = factor(symbol, levels = rev(activation_genes)),
    label = paste0(sprintf("%+.1f", diff), if_else(called, "*", ""))
  )
p_act <- ggplot(act, aes(condition_lab, symbol, fill = diff)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label, color = abs(diff) > 2), size = 3.4, show.legend = FALSE) +
  scale_fill_gradient2(
    low = gradient_low, mid = "white", high = gradient_high,
    midpoint = 0, limits = c(-5, 5), oob = scales::squish, name = "Longevity\nlog2 diff"
  ) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey10")) +
  labs(
    title = "Activated-fibroblast genes",
    subtitle = "Star: fc1.1 call. CCN2 is higher at 32C, 37C, and 41C.",
    x = NULL, y = NULL,
    caption = "COL1A1, COL1A2, POSTN, ACTA2, CCN2, ITGA11, LOXL2. LOX is absent."
  ) +
  theme_tile(16)
save_fig(p_act, "activation_genes", 8.4, 4.8)

# External concordance is the correlation heatmap
# longevity/results/static_genes/plots/corr_fibro_ma_sig47.png
# (fibroblast conditions, Ma ML and MLres, Signatures 4.7).

# --- 14 generic stress -------------------------------------------------------
generic <- read_csv(
  file.path(root, "generic_stress/results/mammalian_generic_stress/gene_human_mouse_generic.csv"),
  show_col_types = FALSE
) %>%
  filter(in_generic %in% TRUE, symbol != "") %>%
  transmute(symbol = toupper(symbol), generic_dir = if_else(direction == "up", "Generic up", "Generic down")) %>%
  distinct(symbol, .keep_all = TRUE)
gen_med <- static %>%
  inner_join(generic, by = "symbol") %>%
  group_by(condition_lab, generic_dir) %>%
  summarise(median = median(diff), n = n(), .groups = "drop")
p_gen <- ggplot(gen_med, aes(condition_lab, median, color = generic_dir, group = generic_dir)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.6) +
  scale_color_manual(
    values = c("Generic up" = unname(treatment_colors["heat"]),
               "Generic down" = unname(treatment_colors["hypoxia"])),
    name = NULL
  ) +
  labs(
    title = "Generic-stress genes across conditions",
    subtitle = "Median longevity difference. Stress-induced genes sit lower at baseline temperatures.",
    x = NULL, y = "Median longevity log2 difference"
  ) +
  theme_pub(16) + theme(legend.position = "bottom")
save_fig(p_gen, "generic_stress", 8.2, 5.2)

# --- 15 GenAge ---------------------------------------------------------------
genage <- read_csv(file.path(root, "longevity/results/static_genes/genage_set_wilcoxon.csv"),
                   show_col_types = FALSE) %>%
  mutate(
    class_lab = recode(class,
      mouse_pro = "Mouse pro", mouse_anti = "Mouse anti",
      other_pro = "Worm/fly/yeast pro", other_anti = "Worm/fly/yeast anti"
    ),
    class_lab = factor(class_lab, levels = c("Mouse pro", "Mouse anti",
                                            "Worm/fly/yeast pro", "Worm/fly/yeast anti")),
    condition_lab = factor(condition_lab, levels = unname(static_labs))
  )
p_gage <- ggplot(genage, aes(condition_lab, median, color = class_lab, group = class_lab)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_color_manual(
    values = c(
      "Mouse pro" = unname(treatment_colors["heat"]),
      "Mouse anti" = unname(treatment_colors["hypoxia"]),
      "Worm/fly/yeast pro" = "#E8A0A0",
      "Worm/fly/yeast anti" = "#9DC3E0"
    ),
    name = NULL
  ) +
  labs(
    title = "GenAge classes across conditions",
    subtitle = "Median longevity difference. No class separates from the rest of the genome.",
    x = NULL, y = "Median longevity log2 difference"
  ) +
  theme_pub(15) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_gage, "genage", 8.4, 5.4)

# --- 16 HMP and gene panel ---------------------------------------------------
hmp <- read_csv(file.path(root, "longevity/results/hmp/static_hmp.csv"), show_col_types = FALSE) %>%
  filter(in_set %in% TRUE) %>%
  mutate(
    direction = factor(direction, levels = c("higher", "lower")),
    symbol = reorder(symbol, effect_mean)
  )
p_hmp <- ggplot(hmp, aes(effect_mean, symbol, color = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_segment(aes(x = 0, xend = effect_mean, yend = symbol), linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = dir_colors,
    labels = c(higher = "Higher at higher LQ", lower = "Lower at higher LQ"),
    name = NULL
  ) +
  labs(
    title = "Genes that pass the cross-condition harmonic mean",
    subtitle = "BH-adjusted harmonic mean p < 0.05, at least 3 conditions, signs agree.",
    x = "Mean longevity log2 difference", y = NULL
  ) +
  theme_pub(16) + theme(legend.position = "bottom")
save_fig(p_hmp, "hmp_genes", 7.6, 5.0)

panel <- bind_rows(
  tibble(symbol = activation_genes, group = "Activation"),
  tibble(symbol = c("MAP3K4", "ACIN1", "ZDHHC2", "MEA1", "NDUFA2", "AGPAT5"), group = "Called in all 5"),
  tibble(symbol = as.character(hmp$symbol), group = "Harmonic mean")
) %>% distinct()
panel_df <- static %>%
  inner_join(panel, by = "symbol") %>%
  mutate(
    group = factor(group, levels = c("Activation", "Called in all 5", "Harmonic mean")),
    label = paste0(sprintf("%+.1f", diff), if_else(called, "*", ""))
  )
p_panel <- ggplot(panel_df, aes(condition_lab, symbol, fill = diff)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = label, color = abs(diff) > 2.2), size = 2.7, show.legend = FALSE) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = gradient_low, mid = "white", high = gradient_high,
    midpoint = 0, limits = c(-5, 5), oob = scales::squish, name = "Longevity\nlog2 diff"
  ) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey10")) +
  labs(
    title = "Genes carried through the summary",
    subtitle = "Star: fc1.1 call. A gene can appear in more than one block.",
    x = NULL, y = NULL
  ) +
  theme_tile(15) +
  theme(strip.text.y = element_text(angle = 0))
save_fig(p_panel, "gene_panel", 8.8, 8.4)

# --- 18 reactions ------------------------------------------------------------
rx_files <- c(
  "32C" = "longevity/results/reactions_static/longevity_reactions_static_32C.csv",
  "37C" = "longevity/results/reactions_static/longevity_reactions_static_37C.csv",
  "41C" = "longevity/results/reactions_static/longevity_reactions_static_41C.csv",
  "hypoxia_6h" = "longevity/results/reactions_static/longevity_reactions_static_hypoxia_6h.csv",
  "hypoxia_24h" = "longevity/results/reactions_static/longevity_reactions_static_hypoxia_24h.csv"
)
rx <- imap_dfr(rx_files, function(path, cond) {
  read_csv(file.path(root, path), show_col_types = FALSE) %>%
    filter(effecttype_drfp0_05 %in% c("higher_in_long_lived", "lower_in_long_lived")) %>%
    transmute(
      direction = if_else(effecttype_drfp0_05 == "higher_in_long_lived", "higher", "lower"),
      condition_lab = factor(static_labs[cond], levels = unname(static_labs))
    )
}) %>% count(condition_lab, direction)
p_rx <- ggplot(rx, aes(condition_lab, n, fill = direction)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66,
           color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_dodge(width = 0.72), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = dir_colors, labels = c(higher = "Higher", lower = "Lower"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Reaction calls at each condition",
    subtitle = "Relative flux potential, drfp 0.05. Higher means greater flux potential at the higher-longevity anchor.",
    x = NULL, y = "Reactions"
  ) +
  theme_pub(16) + theme(legend.position = "bottom")
save_fig(p_rx, "reaction_calls", 8.2, 5.2)

# --- 19 species PGLS leave-one-out -------------------------------------------
pgls <- read_csv(file.path(root, "longevity/results/static_genes/fibroblast_species_pgls.csv"),
                 show_col_types = FALSE)
full_fit <- pgls %>% filter(dropped == "none")
loo <- pgls %>%
  filter(dropped != "none") %>%
  mutate(
    score_lab = recode(score_name, activation = "Activation program",
                       basement_membrane = "Basement membrane"),
    species_lab = pretty_species(dropped)
  )
full_seg <- full_fit %>%
  mutate(score_lab = recode(score_name, activation = "Activation program",
                            basement_membrane = "Basement membrane"))
p_loo <- ggplot(loo, aes(slope, reorder(species_lab, slope))) +
  geom_vline(data = full_seg, aes(xintercept = slope),
             color = unname(treatment_colors["temperature"]), linewidth = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_point(size = 2.5, color = "#3F545E") +
  facet_wrap(~score_lab, scales = "free_x") +
  labs(
    title = "Leave-one-species-out slope",
    subtitle = "PGLS of the 37C signature score on residual lifespan. Red line: all 12 species.",
    x = "PGLS slope", y = "Species dropped"
  ) +
  theme_pub(15)
save_fig(p_loo, "species_loo", 9.0, 5.6)

talk_map <- c(
  "levels_plus_responses" = "levels_plus_responses",
  "ma_shared_scatter" = "ma_shared_scatter",
  "consistent_pathways" = "stable_pathways",
  "response_signatures" = "response_pathways",
  "dynamic_effect_types" = "response_direction"
)
old_talk <- list.files(talk, full.names = TRUE)
if (length(old_talk)) file.remove(old_talk)
for (src in names(talk_map)) {
  for (ext in c("png", "pdf")) {
    file.copy(
      file.path(fig, paste0(src, ".", ext)),
      file.path(talk, paste0(talk_map[[src]], ".", ext)),
      overwrite = TRUE
    )
  }
}
heatmap <- file.path(root, "longevity/results/static_genes/plots/corr_fibro_ma_sig47")
for (ext in c("png", "pdf")) {
  file.copy(
    paste0(heatmap, ".", ext),
    file.path(talk, paste0("levels_plus_raw_lifespan.", ext)),
    overwrite = TRUE
  )
}

source(file.path(root, "longevity/plots/plot_reaction_detail.R"))

message("Supplement written to ", fig)
message("Talk figures written to ", talk)
