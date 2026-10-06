#!/usr/bin/env Rscript
# Public-generic vs Drive figures. Written to plots/mammalian_generic_stress/.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
source("scripts/shared/plot_aesthetics.R")
source("scripts/shared/geo_helpers.R")
source("scripts/generic_stress/plot_helpers.R")

RES <- file.path(ROOT, "results/mammalian_generic_stress")

species_levels <- c(
  "human", "gelada", "rat", "squirrel", "little_brown_bat",
  "honey_badger", "seal", "dolphin", "whale", "rhino",
  "dromedary_camel", "bactrian_camel"
)
exp_levels <- c("Hypoxia 6H", "Hypoxia 24H", "Cold 32C", "Heat 41C")
call_cols <- c(Induced = "#7CB577", Unclear = "#D0D4D8", Opposite = "#2166AC")
tile_theme <- theme_pub(13) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

honest <- read_csv(file.path(RES, "drive_gsea_honest_holdout.csv"), show_col_types = FALSE) %>%
  filter(!exposure %in% c("glucose_low", "glucose_high"), species != "rousette") %>%
  mutate(
    species_lab = factor(pretty_species(species), levels = pretty_species(species_levels)),
    exp_lab = factor(pretty_exp(exposure), levels = exp_levels),
    call = factor(
      case_when(
        NES_up > 0 & padj_up < 0.05 ~ "Induced",
        NES_up < 0 & padj_up < 0.05 ~ "Opposite",
        TRUE ~ "Unclear"
      ),
      levels = c("Induced", "Unclear", "Opposite")
    )
  )

hold <- read_csv(file.path(RES, "drive_gsea_class_holdout.csv"), show_col_types = FALSE) %>%
  filter(!exposure %in% c("glucose_low", "glucose_high"), species != "rousette")
asg <- read_csv(file.path(RES, "drive_conserved_assignment.csv"), show_col_types = FALSE) %>%
  filter(!exposure %in% c("glucose_low", "glucose_high"))
mam <- read_csv(file.path(RES, "gene_human_mouse_generic.csv"), show_col_types = FALSE)
gsea_pw <- read_csv(file.path(RES, "cell_meta_gsea_hallmark_reactome_gobp.csv"),
                    show_col_types = FALSE)
jac <- read_csv(file.path(RES, "cell_meta_celltype_jaccard.csv"), show_col_types = FALSE)
inv <- read_csv(file.path(RES, "cell_meta_contrasts.csv"), show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Discovery inventory
# -----------------------------------------------------------------------------
inv2 <- inv %>%
  count(stress_class, source, name = "n") %>%
  mutate(
    class_lab = factor(
      recode(stress_class, heat = "Heat", hypoxia = "Hypoxia",
             oxidative = "Oxidative (H2O2)", dna_damage = "DNA damage (UV/IR)"),
      levels = c("Heat", "Hypoxia", "Oxidative (H2O2)", "DNA damage (UV/IR)")
    )
  )
p_disc <- ggplot(inv2, aes(class_lab, n, fill = source)) +
  geom_col(width = 0.7, color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            fontface = "bold", size = 3.6, color = "grey10") +
  scale_fill_manual(values = c(ASTRA = unname(treatment_colors["stress"]),
                               GEO = unname(treatment_colors["other"])),
                    name = NULL) +
  labs(
    title = "Public discovery is four ASTRA classes only",
    subtitle = "108 contrasts / 59 GSE. ASTRA = heat, H2O2, hypoxia, UV. Extra GEO fills mouse + a few human fibroblasts.",
    x = NULL, y = "Contrasts"
  ) +
  theme_pub(13) +
  theme(legend.position = "bottom")
save_png(p_disc, "discovery_classes.png", 8.2, 5.2)

# -----------------------------------------------------------------------------
# Pathways
# -----------------------------------------------------------------------------
is_cycle <- function(nm) {
  str_detect(
    tolower(nm),
    "mitot|chromatid|chromosome|spindle|e2f|g2m|cell cycle|prometaphase|metaphase|meiot|nuclear division|formins"
  )
}
pw_h <- gsea_pw %>% filter(collection == "Hallmark", padj < 0.05)
pw_r <- gsea_pw %>% filter(collection == "Reactome", padj < 0.05) %>% arrange(padj)
pw_g <- gsea_pw %>% filter(collection == "GO BP", padj < 0.05) %>% arrange(padj)
pw <- bind_rows(
  pw_h,
  bind_rows(
    pw_r %>% filter(!is_cycle(name)) %>% slice_head(n = 5),
    pw_r %>% filter(is_cycle(name)) %>% slice_head(n = 2)
  ),
  bind_rows(
    pw_g %>% filter(!is_cycle(name)) %>% slice_head(n = 4),
    pw_g %>% filter(is_cycle(name)) %>% slice_head(n = 1)
  )
) %>%
  mutate(
    collection = factor(collection, levels = c("Hallmark", "Reactome", "GO BP")),
    stars = get_sig_star(padj),
    axis = if_else(
      duplicated(name) | duplicated(name, fromLast = TRUE),
      paste0(name, " · ", collection),
      name
    ),
    axis = str_wrap(axis, 42)
  ) %>%
  arrange(NES) %>%
  mutate(axis = fct_inorder(axis))
p_pw <- ggplot(pw, aes(NES, axis)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_segment(aes(x = 0, xend = NES, y = axis, yend = axis),
               color = "grey70", linewidth = 0.65) +
  geom_point(aes(fill = NES, shape = collection),
             color = "grey25", size = 3.4, stroke = 0.5) +
  geom_text(aes(label = stars, x = NES + 0.16 * sign(NES)), size = 2.7) +
  scale_fill_gradient2(low = gradient_low, mid = "white", high = gradient_high,
                       midpoint = 0, guide = "none") +
  scale_shape_manual(values = c(Hallmark = 21, Reactome = 22, `GO BP` = 23),
                     name = NULL) +
  labs(
    title = "Public generic CSR: mitosis off, p53 / hypoxia / ISR on",
    subtitle = "Hallmark, Reactome, and GO BP combined. Redundant cell-cycle terms collapsed.",
    x = "NES", y = NULL, caption = sig_caption
  ) +
  theme_pub(12) +
  theme(legend.position = "bottom")
save_png(p_pw, "public_pathways.png", 8.4, 8.8)

# -----------------------------------------------------------------------------
# Shared gene universe
# -----------------------------------------------------------------------------
drive_files <- list_drive_dge()
drive_files <- drive_files[!str_detect(drive_files, "glucose")]
drive_genes <- map(drive_files, function(f) {
  d <- read_csv(f, show_col_types = FALSE)
  g <- if ("gene" %in% names(d)) d$gene else d[[1]]
  unique(strip(g)[str_starts(strip(g), "ENSG")])
})
drive_union <- unique(unlist(drive_genes))
pub_all <- unique(mam$gene)
pub_gen <- unique(mam$gene[mam$in_generic])
univ <- tibble(
  set = factor(
    c("Public meta genes", "Public generic (1,344)",
      "In Drive (any species)", "Generic in Drive"),
    levels = c("Public meta genes", "Public generic (1,344)",
               "In Drive (any species)", "Generic in Drive")
  ),
  n = c(
    length(pub_all),
    length(pub_gen),
    length(intersect(pub_all, drive_union)),
    length(intersect(pub_gen, drive_union))
  ),
  note = c(
    "tested in public HMP table",
    "HMP ≥3/4 classes",
    sprintf("%.0f%% of public meta", 100 * length(intersect(pub_all, drive_union)) / length(pub_all)),
    sprintf("%.0f%% of generic list", 100 * length(intersect(pub_gen, drive_union)) / length(pub_gen))
  )
)
p_univ <- ggplot(univ, aes(set, n)) +
  geom_col(width = 0.62, fill = unname(treatment_colors["stress"]),
           color = "grey25", linewidth = 0.35) +
  geom_text(aes(label = paste0(comma(n), "\n", note)), vjust = -0.15,
            fontface = "bold", size = 3.3, lineheight = 1.05) +
  coord_cartesian(ylim = c(0, max(univ$n) * 1.22)) +
  labs(
    title = "Gene universes overlap — we are not losing the comparison to IDs",
    subtitle = "74% of public-tested genes and 88% of the 1,344 generic genes appear in at least one Drive file.",
    x = NULL, y = "Genes"
  ) +
  theme_pub(13)
save_png(p_univ, "gene_universe.png", 8.8, 5.6)

# -----------------------------------------------------------------------------
# Conserved Drive genes vs public generic
# -----------------------------------------------------------------------------
dir_map <- mam %>% filter(in_generic) %>% distinct(gene, direction)
asg2 <- asg %>%
  left_join(dir_map, by = "gene") %>%
  mutate(
    exp_lab = factor(pretty_exp(exposure), levels = exp_levels),
    bin = case_when(
      assignment == "generic_mammalian_stress" &
        ((direction == "up" & lfc_med > 0) | (direction == "down" & lfc_med < 0)) ~
        "Public generic, same direction",
      assignment == "generic_mammalian_stress" ~
        "Public generic, opposite direction",
      assignment == "matched_class_specific" ~ "Public class-specific",
      TRUE ~ "Drive residual"
    ),
    bin = factor(bin, levels = c(
      "Public generic, same direction",
      "Public generic, opposite direction",
      "Public class-specific",
      "Drive residual"
    ))
  )
asg_n <- asg2 %>% count(exp_lab, bin, name = "n") %>%
  group_by(exp_lab) %>%
  mutate(tot = sum(n), frac = n / tot) %>%
  ungroup()
bin_cols <- c(
  "Public generic, same direction" = unname(treatment_colors["stress"]),
  "Public generic, opposite direction" = unname(treatment_colors["heat"]),
  "Public class-specific" = unname(treatment_colors["oxidative"]),
  "Drive residual" = "#B0B7BE"
)
p_cons <- ggplot(asg_n, aes(exp_lab, n, fill = bin)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.3) +
  geom_text(
    data = asg_n %>% filter(bin == "Public generic, same direction"),
    aes(label = sprintf("%d same-dir\n(%.0f%% of conserved)", n, 100 * frac),
        y = tot + 80),
    size = 3.1, fontface = "bold", color = "grey20", lineheight = 0.95
  ) +
  scale_fill_manual(values = bin_cols, name = NULL) +
  labs(
    title = "Most conserved Drive genes are not the public generic",
    subtitle = "Conserved = ≥5 species, 80% agree, |median LFC| ≥ 0.25. Grey = Drive-only. Red = on the public list but opposite sign.",
    x = NULL, y = "Conserved genes"
  ) +
  theme_pub(13) +
  theme(legend.position = "bottom")
save_png(p_cons, "conserved_overlap.png", 9.2, 5.8)

# -----------------------------------------------------------------------------
# Honest holdout scorecard + class-composition strip
# -----------------------------------------------------------------------------
score <- honest %>% mutate(species_lab = fct_rev(species_lab))
p_score <- ggplot(score, aes(exp_lab, species_lab, fill = call)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = recode(as.character(call),
                               Induced = "YES", Unclear = "ns", Opposite = "NO")),
            fontface = "bold", size = 3.3) +
  scale_fill_manual(values = call_cols, name = NULL) +
  labs(
    title = "Does Drive induce the class-held-out public up-genes?",
    subtitle = "Coloured cells in the strip = public classes that enter the HMP gene set. Grey = held out. Cold uses all four (includes heat).",
    x = NULL, y = NULL
  ) +
  tile_theme
p_hold <- stack_strip(
  class_composition_strip(c("hypoxia_6H", "hypoxia_24H", "cold", "heat")),
  p_score
)
save_png(p_hold, "holdout_scorecard.png", 8.8, 7.6)

# -----------------------------------------------------------------------------
# Unmatched ox + DNA, every Drive exposure
# -----------------------------------------------------------------------------
unm <- hold %>%
  filter(pathway == "unmatched_up") %>%
  mutate(
    species_lab = factor(pretty_species(species),
                         levels = rev(pretty_species(species_levels))),
    exp_lab = factor(pretty_exp(exposure), levels = exp_levels),
    call = factor(
      case_when(
        NES > 0 & padj < 0.05 ~ "Induced",
        NES < 0 & padj < 0.05 ~ "Opposite",
        TRUE ~ "Unclear"
      ),
      levels = c("Induced", "Unclear", "Opposite")
    )
  )
p_unm_tiles <- ggplot(unm, aes(exp_lab, species_lab, fill = call)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = recode(as.character(call),
                               Induced = "YES", Unclear = "ns", Opposite = "NO")),
            fontface = "bold", size = 3.3) +
  scale_fill_manual(values = call_cols, name = NULL) +
  labs(
    title = "Same gene set on every column: public H2O2 + UV up-genes only",
    subtitle = "Heat and hypoxia are both held out. The strip is identical across Drive contrasts.",
    x = NULL, y = NULL
  ) +
  tile_theme
unm_fn <- function(exposure) holdout_included("unmatched")
p_unm <- stack_strip(
  class_composition_strip(c("hypoxia_6H", "hypoxia_24H", "cold", "heat"), unm_fn),
  p_unm_tiles
)
save_png(p_unm, "unmatched_ox_dna_scorecard.png", 8.8, 7.6)

# -----------------------------------------------------------------------------
# Cold forest
# -----------------------------------------------------------------------------
cold <- honest %>%
  filter(exposure == "cold") %>%
  mutate(species_lab = fct_reorder(species_lab, NES_up))
p_cold <- ggplot(cold, aes(NES_up, species_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_segment(aes(x = 0, xend = NES_up, y = species_lab, yend = species_lab),
               color = "grey70", linewidth = 0.65) +
  geom_point(aes(fill = call), shape = 21, color = "grey25", size = 3.8, stroke = 0.65) +
  scale_fill_manual(values = call_cols, name = NULL, drop = FALSE) +
  labs(
    title = sprintf("Cold: %d of 12 species induce the full public generic-up list",
                    sum(cold$call == "Induced")),
    subtitle = "GSEA NES of Drive cold (ranked −log10(p)×sign(LFC)) vs the 737 public generic-up genes. Cold is not in ASTRA, so the gene set includes public heat.",
    x = "GSEA NES (public generic-up genes)", y = NULL
  ) +
  theme_pub(13) +
  theme(legend.position = "bottom")
save_png(p_cold, "cold_forest.png", 8.0, 5.8)

# -----------------------------------------------------------------------------
# Human vs other mammals — NES of the holdout generic-up set
# -----------------------------------------------------------------------------
p_hum <- ggplot(honest, aes(exp_lab, NES_up)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_jitter(
    data = filter(honest, species != "human"),
    width = 0.14, height = 0, shape = 21, fill = "grey78",
    color = "grey30", size = 2.8, stroke = 0.45
  ) +
  geom_point(
    data = filter(honest, species == "human"),
    shape = 21, fill = unname(treatment_colors["stress"]),
    color = "grey20", size = 4.4, stroke = 0.7
  ) +
  labs(
    title = "Human Drive is not a better match than other mammals",
    subtitle = "Purple = human; grey = other 11 species. Y is GSEA NES of Drive (ranked −log10(p)×sign(LFC)) vs the held-out public generic-up genes. Strip = public classes in that gene set.",
    x = NULL, y = "GSEA NES vs held-out public generic-up genes"
  ) +
  theme_pub(12)
p_hum2 <- stack_strip(
  class_composition_strip(c("hypoxia_6H", "hypoxia_24H", "cold", "heat")),
  p_hum,
  heights = c(1.05, 3.6)
)
save_png(p_hum2, "human_vs_others_holdout_nes.png", 8.4, 6.8)

# -----------------------------------------------------------------------------
# REML / cell-type sensitivity
# -----------------------------------------------------------------------------
jac2 <- jac %>%
  filter(a == "dl_full", b != "dl_full") %>%
  mutate(
    lab = recode(
      b,
      reml_noc2c12 = "REML instead of DL\n(drop C2C12)",
      reml_nocancer = "REML, drop cancer lines\nand C2C12",
      reml_fibroblast = "REML, fibroblast-like\ncells only",
      rmamv_noc2c12 = "rma.mv REML\n(drop C2C12)"
    ),
    lab = fct_reorder(lab, jaccard)
  )
p_jac <- ggplot(jac2, aes(jaccard, lab)) +
  geom_col(width = 0.62, fill = unname(treatment_colors["stress"]),
           color = "grey25", linewidth = 0.35) +
  geom_text(aes(label = sprintf("Jaccard %.2f   %d of %d genes in this list are in the 1,344",
                                jaccard, n_shared, n_b)),
            hjust = -0.02, fontface = "bold", size = 3.2) +
  coord_cartesian(xlim = c(0, 1.55)) +
  labs(
    title = "The 1,344-gene list is stable in direction, not in membership",
    subtitle = "Each bar is an alternative generic list vs the published 1,344. Jaccard 1 = identical set. No genes flip sign under REML. Fibroblast-only hypoxia is only 3 GSE.",
    x = "Jaccard overlap with the 1,344-gene list", y = NULL
  ) +
  theme_pub(12)
save_png(p_jac, "reml_celltype_overlap.png", 9.4, 5.2)

writeLines(c(
  "# Generic-stress figures",
  "",
  "This folder is the paper set. No locked/ numbering.",
  "",
  "discovery_classes — public ASTRA/GEO contrasts by class",
  "public_pathways — Hallmark + Reactome + GO BP GSEA of the 1,344",
  "gene_universe — public vs Drive ID overlap",
  "conserved_overlap — conserved Drive genes vs public generic",
  "holdout_scorecard — Drive vs class-held-out public up-genes; strip = classes in the HMP set",
  "unmatched_ox_dna_scorecard — same H2O2+UV gene set on every Drive contrast",
  "cold_forest — GSEA NES of Drive cold vs the 737 public generic-up genes",
  "human_vs_others_holdout_nes — GSEA NES of each Drive contrast vs the held-out public generic-up set",
  "reml_celltype_overlap — Jaccard of alternative generic lists vs the 1,344",
  "volcano_<species> — Drive DEGs; labels = padj<0.05 AND public generic",
  "drive_lfc_density — median Drive LFC, generic vs other genes",
  "species_hmp_counts — per-species heat+cold+hypoxia HMP counts",
  "species_hmp_vs_public — numerator = genes in public 1,344; denominator = that species HMP list",
  "species_hmp_conservation — how many Drive species call each gene",
  "species_hmp_conserved_ora — ORA of the ≥5-species core"
), file.path(PLOT_DIR, "README.md"))

cat("Figures written to ", PLOT_DIR, "\n", sep = "")
