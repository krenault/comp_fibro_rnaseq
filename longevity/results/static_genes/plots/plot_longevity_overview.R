#!/usr/bin/env Rscript
# Static longevity levels, dynamic responses, and overlap with the
# pan-cell generic-stress signature. Style matches plot_aesthetics.R.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

root <- "/Users/katiarenault/Documents/GitHub/comp_fibro_rnaseq"
source(file.path(root, "generic_stress/scripts/shared/plot_aesthetics.R"))

plot_dir <- file.path(root, "longevity/results/static_genes/plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, file, w, h) {
  ggsave(file.path(plot_dir, file), p, width = w, height = h, dpi = PUB_DPI, bg = "white")
  message("Wrote ", file)
}

static_conds <- c("32C", "37C", "41C", "hypoxia_6h", "hypoxia_24h")
static_labs <- c("32C" = "32C", "37C" = "37C", "41C" = "41C",
                 "hypoxia_6h" = "Hypoxia 6H", "hypoxia_24h" = "Hypoxia 24H")
dyn_conds <- c("32C_vs_37C", "41C_vs_37C", "hypoxia_6h_vs_0h", "hypoxia_24h_vs_0h")
dyn_labs <- c(
  "32C_vs_37C" = "Cold 32C",
  "41C_vs_37C" = "Heat 41C",
  "hypoxia_6h_vs_0h" = "Hypoxia 6H",
  "hypoxia_24h_vs_0h" = "Hypoxia 24H"
)

read_static <- function(cond) {
  read_csv(file.path(root, "longevity/data", paste0("longevity_genes_static_", cond, ".csv")),
           show_col_types = FALSE) %>%
    transmute(
      symbol = toupper(hgnc_symbol),
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

static <- map_dfr(static_conds, read_static)

read_dyn <- function(cond) {
  read_csv(file.path(root, "longevity/results/dynamic_genes",
                     paste0("longevity_genes_dynamic_", cond, ".csv")),
           show_col_types = FALSE) %>%
    transmute(
      symbol = toupper(hgnc_symbol),
      called = effecttype_fc1_1 != "",
      condition = cond
    ) %>%
    filter(symbol != "")
}
dynamic <- map_dfr(dyn_conds, read_dyn)

nes_of <- function(path, keep_path, keep_analysis, lab_map) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(pathway %in% keep_path, analysis %in% keep_analysis, !is.na(NES)) %>%
    mutate(
      pathway = factor(pathway, levels = keep_path),
      analysis = factor(analysis, levels = keep_analysis),
      analysis_lab = factor(lab_map[as.character(analysis)], levels = unname(lab_map)),
      star = get_sig_star(padj),
      label = paste0(sprintf("%+.1f", NES), star)
    )
}

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

static_nes <- nes_of(
  file.path(root, "longevity/results/static_genes/fgsea_msigdb_v2023_2_all.csv"),
  names(static_paths), static_conds, static_labs
) %>%
  mutate(pathway_lab = factor(static_paths[as.character(pathway)], levels = rev(unname(static_paths))))

dyn_nes <- nes_of(
  file.path(root, "longevity/results/dynamic_genes/fgsea_msigdb_v2023_2_all.csv"),
  names(dyn_paths), dyn_conds, dyn_labs
) %>%
  mutate(pathway_lab = factor(dyn_paths[as.character(pathway)], levels = rev(unname(dyn_paths))))

heat_plot <- function(df, title, subtitle, xlab) {
  lim <- max(abs(df$NES), na.rm = TRUE)
  ggplot(df, aes(analysis_lab, pathway_lab, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label, color = abs(NES) > 1.7), size = 3.3, show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey10")) +
    scale_fill_gradient2(
      low = gradient_low, mid = "white", high = gradient_high,
      midpoint = 0, limits = c(-lim, lim), name = "NES"
    ) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL,
         caption = paste0(sig_caption, ".  ", xlab)) +
    theme_pub(12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      legend.position = "right"
    )
}

p_static <- heat_plot(
  static_nes,
  "Longevity-associated expression level",
  "Static model. Positive NES: higher at the higher-longevity anchor.",
  "Each column is one culture condition."
)
p_dyn <- heat_plot(
  dyn_nes,
  "Longevity-associated stress response",
  "Dynamic model. Positive NES: a more positive fold-change at the higher-longevity anchor.",
  "Each column is a perturbation relative to its baseline."
)

# --- generic stress ----------------------------------------------------------
generic <- read_csv(
  file.path(root, "generic_stress/results/mammalian_generic_stress/gene_human_mouse_generic.csv"),
  show_col_types = FALSE
) %>%
  transmute(
    symbol = toupper(symbol),
    in_generic = in_generic %in% TRUE,
    generic_dir = if_else(in_generic %in% TRUE, direction, "not generic"),
    ec_mean = ec_mean
  ) %>%
  filter(symbol != "") %>%
  distinct(symbol, .keep_all = TRUE)

joined <- static %>%
  inner_join(generic, by = "symbol")

rho_tbl <- joined %>%
  group_by(condition) %>%
  summarise(
    rho = cor(diff, ec_mean, method = "spearman", use = "complete.obs"),
    n = sum(is.finite(diff) & is.finite(ec_mean)),
    .groups = "drop"
  ) %>%
  mutate(condition_lab = factor(static_labs[condition], levels = unname(static_labs)))

# Called longevity genes vs generic membership, 37C and each condition
call37 <- static %>% filter(condition == "37C", called) %>% pull(symbol)
overlap_rows <- list()
for (cond in static_conds) {
  calls <- static %>% filter(condition == cond, called)
  g <- generic %>% filter(in_generic)
  both <- intersect(calls$symbol, g$symbol)
  # sign agreement: generic up with higher longevity, generic down with lower
  m <- calls %>%
    inner_join(g, by = "symbol") %>%
    mutate(agree = (direction == "higher" & generic_dir == "up") |
                    (direction == "lower" & generic_dir == "down"))
  overlap_rows[[cond]] <- tibble(
    condition = cond,
    n_calls = nrow(calls),
    n_generic_calls = nrow(m),
    n_agree = sum(m$agree),
    frac_calls_generic = nrow(m) / nrow(calls)
  )
}
overlap_tbl <- bind_rows(overlap_rows)
write_csv(overlap_tbl, file.path(plot_dir, "generic_stress_overlap_calls.csv"))
write_csv(rho_tbl, file.path(plot_dir, "generic_stress_spearman.csv"))

# Wilcoxon: generic up / down vs everyone else, longevity diff
class_stats <- joined %>%
  mutate(grp = case_when(
    generic_dir == "up" ~ "Generic up",
    generic_dir == "down" ~ "Generic down",
    TRUE ~ "Other genes"
  )) %>%
  group_by(condition, grp) %>%
  summarise(
    n = n(),
    median = median(diff),
    mean = mean(diff),
    .groups = "drop"
  )

p_vals <- joined %>%
  mutate(grp = case_when(
    generic_dir == "up" ~ "Generic up",
    generic_dir == "down" ~ "Generic down",
    TRUE ~ "Other genes"
  )) %>%
  group_by(condition) %>%
  group_modify(~{
    ref <- .x$diff[.x$grp == "Other genes"]
    tibble(
      grp = c("Generic up", "Generic down"),
      p = c(
        wilcox.test(.x$diff[.x$grp == "Generic up"], ref)$p.value,
        wilcox.test(.x$diff[.x$grp == "Generic down"], ref)$p.value
      )
    )
  }) %>%
  ungroup()

message("\nSpearman longevity diff vs generic-stress effect:")
print(rho_tbl)
message("\nCall overlap:")
print(overlap_tbl)
message("\nClass medians:")
print(class_stats %>% filter(grp != "Other genes"))
message("\nWilcox vs other genes:")
print(p_vals)

# Dynamic calls that are generic
dyn_overlap <- dynamic %>%
  filter(called) %>%
  left_join(generic %>% select(symbol, in_generic, generic_dir), by = "symbol") %>%
  group_by(condition) %>%
  summarise(
    n_calls = n(),
    n_generic = sum(in_generic %in% TRUE),
    .groups = "drop"
  )
message("\nDynamic calls in generic signature:")
print(dyn_overlap)

p_rho <- ggplot(rho_tbl, aes(condition_lab, rho, fill = condition == "37C")) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", rho), vjust = ifelse(rho < 0, 1.3, -0.4)),
            fontface = "bold", size = 3.6) +
  scale_fill_manual(values = c(`TRUE` = unname(treatment_colors["heat"]),
                               `FALSE` = unname(treatment_colors["hypoxia"])),
                    guide = "none") +
  scale_y_continuous(limits = c(-0.12, 0.08), expand = c(0, 0)) +
  labs(
    title = "Longevity levels vs the generic-stress effect",
    subtitle = "Spearman correlation across shared genes. Near zero means the two rankings are different.",
    x = NULL, y = "Spearman rho"
  ) +
  theme_pub(13)

box_df <- joined %>%
  filter(condition == "37C") %>%
  mutate(grp = case_when(
    generic_dir == "up" ~ "Generic up",
    generic_dir == "down" ~ "Generic down",
    TRUE ~ "Other genes"
  ))
lab_n <- box_df %>% count(grp)
lab_n <- lab_n %>% mutate(grp_lab = paste0(grp, "\n(n = ", n, ")"))
box_df <- box_df %>%
  left_join(lab_n, by = "grp") %>%
  mutate(grp_lab = factor(grp_lab, levels = lab_n$grp_lab[match(
    c("Generic up", "Generic down", "Other genes"), lab_n$grp
  )]))

p_box <- ggplot(box_df, aes(grp_lab, diff, fill = grp)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_violin(scale = "width", color = NA, alpha = 0.9) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(
    `Generic up` = unname(treatment_colors["heat"]),
    `Generic down` = unname(treatment_colors["hypoxia"]),
    `Other genes` = "grey80"
  ), guide = "none") +
  coord_cartesian(ylim = c(-4, 4)) +
  labs(
    title = "37C longevity difference by stress class",
    subtitle = "Generic up: induced across heat, hypoxia, oxidative stress, and DNA damage.",
    x = NULL,
    y = "Longevity log2 difference"
  ) +
  theme_pub(13)

p_gen <- p_rho / p_box + plot_annotation(
  caption = "Generic signature: 960 pan-cell genes, BH-adjusted harmonic mean p < 0.05, sign agreement across classes."
)
# patchwork caption sits better on the combined plot via plot_annotation theme
p_gen <- p_gen & theme(plot.caption = element_text(size = 9, color = "grey40", hjust = 1))

save_fig(p_static, "longevity_static_nes.png", 9.4, 5.6)
save_fig(p_dyn, "longevity_dynamic_nes.png", 8.8, 5.4)
save_fig(p_gen, "longevity_vs_generic_stress.png", 8.4, 9.2)

# also pdf
ggsave(file.path(plot_dir, "longevity_static_nes.pdf"), p_static, width = 9.4, height = 5.6, bg = "white")
ggsave(file.path(plot_dir, "longevity_dynamic_nes.pdf"), p_dyn, width = 8.8, height = 5.4, bg = "white")
ggsave(file.path(plot_dir, "longevity_vs_generic_stress.pdf"), p_gen, width = 8.4, height = 9.2, bg = "white")
message("done")
