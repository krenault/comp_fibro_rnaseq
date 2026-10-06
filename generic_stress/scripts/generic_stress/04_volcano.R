#!/usr/bin/env Rscript
# Per-species Drive volcanoes. Labels = DESeq2 padj < 0.05 AND public generic.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
source("scripts/shared/plot_aesthetics.R")
source("scripts/shared/geo_helpers.R")
source("scripts/generic_stress/plot_helpers.R")

RES <- file.path(ROOT, "results/mammalian_generic_stress")
exp_levels <- c("Hypoxia 6H", "Hypoxia 24H", "Cold 32C", "Heat 41C")
species_levels <- c(
  "human", "gelada", "rat", "squirrel", "little_brown_bat",
  "honey_badger", "seal", "dolphin", "whale", "rhino",
  "dromedary_camel", "bactrian_camel"
)

mam <- read_csv(file.path(RES, "gene_human_mouse_generic.csv"), show_col_types = FALSE)
gen <- mam %>%
  filter(in_generic) %>%
  transmute(gene, gen_dir = direction, symbol = coalesce(symbol, gene))

drive_files <- list_drive_dge()
drive_files <- drive_files[!str_detect(drive_files, "glucose")]

parse_one <- function(f) {
  bn <- tools::file_path_sans_ext(basename(f))
  tt <- str_extract(bn, "hypoxia|temperature")
  cmp <- str_extract(bn, "(6H_vs_0|24H_vs_0|32C_vs_37C|41C_vs_37C)$")
  if (is.na(tt) || is.na(cmp)) return(NULL)
  sp <- str_remove(bn, paste0("_", tt, "_", cmp, "$"))
  exposure <- case_when(
    tt == "hypoxia" & cmp == "6H_vs_0" ~ "hypoxia_6H",
    tt == "hypoxia" & cmp == "24H_vs_0" ~ "hypoxia_24H",
    grepl("32C", cmp) ~ "cold",
    grepl("41C", cmp) ~ "heat",
    TRUE ~ NA_character_
  )
  d <- read_csv(f, show_col_types = FALSE)
  tibble(
    gene = strip(d$gene),
    lfc = as.numeric(d$log2FoldChange),
    p = as.numeric(d$pvalue),
    padj = as.numeric(d$padj),
    species = sp,
    exposure = exposure
  ) %>%
    filter(is.finite(lfc), str_starts(gene, "ENSG"), !is.na(exposure))
}

drive <- map_dfr(drive_files, parse_one) %>%
  filter(species %in% species_levels)
message("Drive rows: ", nrow(drive), " contrasts: ", n_distinct(paste(drive$species, drive$exposure)))

# -----------------------------------------------------------------------------
# One volcano file per species. Labels are that panel's significant generic DEGs.
# -----------------------------------------------------------------------------
plot_species <- function(sp) {
  d <- drive %>%
    filter(species == sp) %>%
    left_join(gen, by = "gene") %>%
    mutate(
      is_generic = !is.na(gen_dir),
      is_sig = is.finite(padj) & padj < 0.05,
      flag = is_generic & is_sig,
      exp_lab = factor(pretty_exp(exposure), levels = exp_levels),
      y = -log10(pmin(pmax(p, 1e-20), 1))
    )
  lab <- d %>%
    filter(flag, !is.na(symbol), nzchar(symbol), !str_starts(symbol, "ENSG")) %>%
    group_by(exp_lab) %>%
    mutate(score = abs(lfc) * y) %>%
    slice_max(order_by = score, n = 8, with_ties = FALSE) %>%
    ungroup()
  n_flag <- d %>%
    group_by(exp_lab) %>%
    summarise(n = sum(flag), .groups = "drop")
  p <- ggplot(d, aes(lfc, y)) +
    geom_point(
      data = filter(d, !is_generic),
      shape = 21, fill = "grey82", color = "grey70",
      size = 1.05, stroke = 0.12, alpha = 0.6
    ) +
    geom_point(
      data = filter(d, is_generic, !is_sig),
      shape = 21, fill = unname(treatment_colors["stress"]),
      color = "grey35", size = 1.35, stroke = 0.2, alpha = 0.55
    ) +
    geom_point(
      data = filter(d, flag),
      shape = 21, fill = unname(treatment_colors["stress"]),
      color = "grey10", size = 1.9, stroke = 0.28, alpha = 0.9
    ) +
    geom_text_repel(
      data = lab, aes(label = symbol), size = 2.7, fontface = "bold",
      min.segment.length = 0, box.padding = 0.18, max.overlaps = 18,
      segment.color = "grey40"
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    facet_wrap(~ exp_lab, ncol = 2) +
    coord_cartesian(xlim = c(-4.5, 4.5), ylim = c(0, 20)) +
    labs(
      title = paste0(pretty_species(sp), " Drive volcano"),
      subtitle = "Grey = not in the public 1,344. Purple = public generic. Labels = padj < 0.05 AND public generic (top 8 per panel by |LFC| × −log10(p)).",
      x = "Drive log2 fold change",
      y = expression(-log[10](italic(p)))
    ) +
    theme_pub(12)
  save_png(p, paste0("volcano_", sp, ".png"), 9.6, 7.4)
  n_flag
}

counts <- map_dfr(species_levels, function(sp) {
  message("volcano ", sp)
  n <- plot_species(sp)
  n %>% mutate(species = sp)
})
write_csv(counts, file.path(RES, "drive_volcano_generic_counts.csv"))

# -----------------------------------------------------------------------------
# Density of median Drive LFC (generic vs other)
# -----------------------------------------------------------------------------
med <- drive %>%
  group_by(gene, exposure) %>%
  summarise(
    lfc_med = median(lfc, na.rm = TRUE),
    n_sp = n_distinct(species),
    .groups = "drop"
  ) %>%
  filter(n_sp >= 5) %>%
  left_join(gen, by = "gene") %>%
  mutate(
    is_generic = !is.na(gen_dir),
    exp_lab = factor(pretty_exp(exposure), levels = exp_levels),
    set = if_else(is_generic, "Public generic", "Other Drive genes")
  )
p_dens <- ggplot(med, aes(lfc_med, color = set, fill = set)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_density(alpha = 0.18, linewidth = 0.8) +
  scale_color_manual(values = c(`Public generic` = unname(treatment_colors["stress"]),
                                `Other Drive genes` = "grey45"), name = NULL) +
  scale_fill_manual(values = c(`Public generic` = unname(treatment_colors["stress"]),
                               `Other Drive genes` = "grey70"), name = NULL) +
  facet_wrap(~ exp_lab, ncol = 2) +
  coord_cartesian(xlim = c(-2.2, 2.2)) +
  labs(
    title = "Public generic genes are not the extreme Drive DEGs",
    subtitle = "Density of median Drive LFC across species. If this were just CSR, purple would sit in the tails.",
    x = "Median Drive log2 fold change", y = "Density"
  ) +
  theme_pub(13) +
  theme(legend.position = "bottom")
save_png(p_dens, "drive_lfc_density.png", 9.2, 6.6)

cat("Per-species volcanoes written.\n")
