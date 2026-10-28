#!/usr/bin/env Rscript
# GenAge versus the fibroblast longevity ranking.
# Model-organism genes are directional. Human GenAge is not.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

root <- "/Users/katiarenault/Documents/GitHub/comp_fibro_rnaseq"
source(file.path(root, "generic_stress/scripts/shared/plot_aesthetics.R"))
hagr <- "/Users/katiarenault/Documents/GitHub/comparative_reprogramming/data/external/hagr"
out_dir <- file.path(root, "longevity/results/static_genes")
plot_dir <- file.path(out_dir, "plots")

save_fig <- function(p, file, w, h) {
  ggsave(file.path(plot_dir, paste0(file, ".png")), p, width = w, height = h, dpi = PUB_DPI, bg = "white")
  ggsave(file.path(plot_dir, paste0(file, ".pdf")), p, width = w, height = h, bg = "white")
  message("Wrote ", file)
}

activation_genes <- c("COL1A1", "COL1A2", "POSTN", "ACTA2", "CCN2", "LOXL2", "ITGA11")
static_conds <- c("32C", "37C", "41C", "hypoxia_6h", "hypoxia_24h")
static_labs <- c("32C" = "32C", "37C" = "37C", "41C" = "41C",
                 "hypoxia_6h" = "Hypoxia 6H", "hypoxia_24h" = "Hypoxia 24H")
taxa <- c(
  `Caenorhabditis elegans` = "6239",
  `Drosophila melanogaster` = "7227",
  `Mus musculus` = "10090",
  `Saccharomyces cerevisiae` = "559292"
)

read_static <- function(cond) {
  read_csv(file.path(root, "longevity/data", paste0("longevity_genes_static_", cond, ".csv")),
           show_col_types = FALSE) %>%
    transmute(
      symbol = toupper(hgnc_symbol),
      diff = longevity_log2diff,
      called = effecttype_fc1_1 %in% c("higher_in_long_lived", "lower_in_long_lived"),
      condition = cond
    ) %>%
    filter(symbol != "", is.finite(diff))
}
static <- map_dfr(static_conds, read_static)

# --- GenAge models -> human, best-scoring orthologs only ---
gm <- read_csv(file.path(hagr, "genage_models.csv"), show_col_types = FALSE) %>%
  filter(organism %in% names(taxa)) %>%
  transmute(
    taxon = unname(taxa[organism]),
    organism,
    model_symbol = tolower(symbol),
    influence = `longevity influence`,
    group = if_else(organism == "Mus musculus", "mouse", "other")
  )

ali <- fread(file.path(hagr, "alliance_orthology.tsv.gz"),
             select = c("Gene1Symbol", "Gene1SpeciesTaxonID", "Gene2Symbol",
                        "Gene2SpeciesTaxonID", "IsBestScore"))
ali <- ali[IsBestScore == "Yes"]
ali[, t1 := sub(".*:", "", Gene1SpeciesTaxonID)]
ali[, t2 := sub(".*:", "", Gene2SpeciesTaxonID)]
to_human <- rbind(
  ali[t2 == "9606", .(taxon = t1, model_symbol = tolower(Gene1Symbol), Symbol = toupper(Gene2Symbol))],
  ali[t1 == "9606", .(taxon = t2, model_symbol = tolower(Gene2Symbol), Symbol = toupper(Gene1Symbol))]
)
to_human <- unique(to_human[taxon %in% unname(taxa) & Symbol != ""])

mapped <- as.data.table(gm)[to_human, on = c("taxon", "model_symbol"), nomatch = NULL]
mapped <- unique(mapped, by = c("group", "organism", "model_symbol", "Symbol", "influence"))

evidence <- mapped %>%
  as_tibble() %>%
  group_by(Symbol, group) %>%
  summarise(
    n_pro = n_distinct(paste(organism, model_symbol)[influence == "Pro-Longevity"]),
    n_anti = n_distinct(paste(organism, model_symbol)[influence == "Anti-Longevity"]),
    organisms = paste(sort(unique(organism)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(direction = case_when(
    n_pro > 0 & n_anti == 0 ~ "pro",
    n_anti > 0 & n_pro == 0 ~ "anti",
    n_pro > 0 & n_anti > 0 ~ "conflict",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(direction))

human_genes <- toupper(read_csv(file.path(hagr, "genage_human.csv"), show_col_types = FALSE)$symbol)

class_of <- evidence %>%
  filter(direction %in% c("pro", "anti")) %>%
  transmute(symbol = Symbol, class = paste(group, direction, sep = "_"))

message(sprintf(
  "Mapped model records %d; human genes with a direction: mouse pro %d, mouse anti %d, other pro %d, other anti %d, conflicts %d",
  nrow(mapped),
  sum(evidence$group == "mouse" & evidence$direction == "pro"),
  sum(evidence$group == "mouse" & evidence$direction == "anti"),
  sum(evidence$group == "other" & evidence$direction == "pro"),
  sum(evidence$group == "other" & evidence$direction == "anti"),
  sum(evidence$direction == "conflict")
))
message("Mouse GenAge genes: ", sum(gm$group == "mouse"),
        "; mapped to a human symbol: ", n_distinct(mapped$Symbol[mapped$group == "mouse"]))

# --- Calls replicated across conditions ---
call_n <- static %>%
  filter(called) %>%
  distinct(symbol, condition) %>%
  count(symbol, name = "n_conditions")
replicated <- call_n %>% filter(n_conditions >= 2) %>% pull(symbol)
called_37 <- static %>% filter(condition == "37C", called) %>% pull(symbol)
replicated_37 <- intersect(replicated, called_37)

# --- Tests ---
classes <- c("mouse_pro", "mouse_anti", "other_pro", "other_anti")
class_genes <- setNames(lapply(classes, function(cl) class_of$symbol[class_of$class == cl]), classes)

test_shift <- function(df, genes, alternative) {
  x <- df$diff[df$symbol %in% genes]
  y <- df$diff[!df$symbol %in% genes]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 3) {
    return(tibble(n = length(x), median = median(x), mean = mean(x),
                  median_background = median(y), p = NA_real_))
  }
  wt <- wilcox.test(x, y, alternative = alternative, exact = FALSE)
  tibble(n = length(x), median = median(x), mean = mean(x),
         median_background = median(y), p = unname(wt$p.value))
}

alt_of <- c(mouse_pro = "greater", mouse_anti = "less", other_pro = "greater", other_anti = "less")

set_stats <- map_dfr(static_conds, function(cond) {
  df <- filter(static, condition == cond)
  map_dfr(classes, function(cl) {
    test_shift(df, class_genes[[cl]], alt_of[[cl]]) %>%
      mutate(condition = cond, class = cl, alternative = alt_of[[cl]])
  })
}) %>%
  group_by(condition) %>%
  mutate(p_bh = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(star = get_sig_star(p_bh), condition_lab = static_labs[condition])

# Human list: no predicted sign. Signed shift, and whether |diff| is larger.
human_stats <- map_dfr(static_conds, function(cond) {
  df <- filter(static, condition == cond)
  signed <- test_shift(df, human_genes, "two.sided")
  x <- abs(df$diff[df$symbol %in% human_genes])
  y <- abs(df$diff[!df$symbol %in% human_genes])
  wt <- wilcox.test(x, y, alternative = "greater", exact = FALSE)
  hi <- df$diff > quantile(df$diff, 0.9)
  lo <- df$diff < quantile(df$diff, 0.1)
  inn <- df$symbol %in% human_genes
  tail_tab <- matrix(c(sum(inn & (hi | lo)), sum(inn & !(hi | lo)),
                       sum(!inn & (hi | lo)), sum(!inn & !(hi | lo))), nrow = 2)
  ft <- fisher.test(tail_tab, alternative = "greater")
  signed %>%
    mutate(condition = cond, p_abs = unname(wt$p.value),
           n_in_outer_decile = sum(inn & (hi | lo)),
           fisher_or = unname(ft$estimate), fisher_p = unname(ft$p.value))
})

# Membership of the fibroblast sets
member_class <- function(genes) {
  tibble(symbol = genes) %>%
    mutate(
      human = symbol %in% human_genes,
      mouse_pro = symbol %in% class_genes$mouse_pro,
      mouse_anti = symbol %in% class_genes$mouse_anti,
      other_pro = symbol %in% class_genes$other_pro,
      other_anti = symbol %in% class_genes$other_anti
    )
}
activation_hit <- member_class(activation_genes)
replicated_hit <- member_class(replicated_37)

d37 <- static %>% filter(condition == "37C") %>% dplyr::select(symbol, diff, called)
gene_table <- bind_rows(
  class_of,
  tibble(symbol = intersect(human_genes, d37$symbol), class = "human")
) %>%
  distinct() %>%
  left_join(d37, by = "symbol") %>%
  filter(!is.na(diff))

write_csv(evidence, file.path(out_dir, "genage_model_evidence.csv"))
write_csv(set_stats, file.path(out_dir, "genage_set_wilcoxon.csv"))
write_csv(human_stats, file.path(out_dir, "genage_human_tests.csv"))
write_csv(gene_table, file.path(out_dir, "genage_genes_in_ranking.csv"))
write_csv(activation_hit, file.path(out_dir, "genage_activation_membership.csv"))

# --- Plot ---
class_labs <- c(
  mouse_pro = "Mouse pro",
  mouse_anti = "Mouse anti",
  other_pro = "Worm/fly/yeast pro",
  other_anti = "Worm/fly/yeast anti"
)
class_colors <- c(
  "Mouse pro" = unname(treatment_colors["heat"]),
  "Mouse anti" = unname(treatment_colors["hypoxia"]),
  "Worm/fly/yeast pro" = "#E8A0A0",
  "Worm/fly/yeast anti" = "#9DC3E0",
  "Other genes" = "grey75"
)

plot_df <- d37 %>%
  left_join(class_of, by = "symbol") %>%
  mutate(class_lab = unname(class_labs[class])) %>%
  mutate(class_lab = if_else(is.na(class_lab), "Other genes", class_lab)) %>%
  mutate(class_lab = factor(class_lab, levels = c("Mouse pro", "Mouse anti",
                                                  "Worm/fly/yeast pro", "Worm/fly/yeast anti",
                                                  "Other genes")))
# A gene in both mouse and other would have been two rows from the join.
# class_of is one row per symbol-class, so the join duplicates those genes. Keep them
# in each class and drop the duplicate "other" assignment by filtering NA class only once.
plot_df <- plot_df %>% filter(!(class_lab == "Other genes" & symbol %in% class_of$symbol))
lab_n <- plot_df %>% count(class_lab) %>% mutate(lab = paste0(class_lab, "\n(n = ", n, ")"))
plot_df <- plot_df %>% left_join(lab_n %>% dplyr::select(class_lab, lab), by = "class_lab")
plot_df$lab <- factor(plot_df$lab, levels = lab_n$lab)

p_dist <- ggplot(plot_df, aes(lab, diff, fill = class_lab)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_violin(color = NA, width = 0.85, alpha = 0.9) +
  geom_boxplot(width = 0.16, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  coord_cartesian(ylim = c(-4, 4)) +
  scale_fill_manual(values = class_colors, guide = "none") +
  labs(
    title = "37C longevity difference by GenAge class",
    subtitle = "Pro-longevity: activity extends lifespan. Anti-longevity: activity limits it.",
    x = NULL, y = "Longevity log2 difference"
  ) +
  theme_pub(base_size = 12)

med_df <- set_stats %>%
  mutate(
    class_lab = factor(unname(class_labs[class]), levels = unname(class_labs)),
    condition_lab = factor(condition_lab, levels = unname(static_labs)),
    label = paste0(sprintf("%+.2f", median), star)
  )
p_med <- ggplot(med_df, aes(condition_lab, median, color = class_lab, group = class_lab)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.4) +
  scale_color_manual(values = class_colors[unname(class_labs)], name = NULL) +
  labs(
    title = "Same classes in every condition",
    subtitle = "Stars: BH-adjusted one-sided Wilcoxon within that condition. Positive median: higher at the higher-longevity anchor.",
    x = NULL, y = "Median longevity log2 difference"
  ) +
  theme_pub(base_size = 12) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

save_fig(p_dist / p_med + plot_layout(heights = c(1.15, 1)), "genage_longevity", 9.4, 8.6)

message("\n37C directional tests:")
print(set_stats %>% filter(condition == "37C") %>% dplyr::select(class, n, median, p, p_bh))
message("\nAll conditions:")
print(set_stats %>% dplyr::select(class, condition, n, median, p, p_bh) %>% arrange(class, condition))
message("\nHuman GenAge:")
print(human_stats)
message("\nActivation genes in GenAge:")
print(activation_hit)
message("\nReplicated 37C calls in GenAge: ",
        sum(replicated_hit$human | replicated_hit$mouse_pro | replicated_hit$mouse_anti |
              replicated_hit$other_pro | replicated_hit$other_anti),
        " of ", nrow(replicated_hit))
print(replicated_hit %>% filter(human | mouse_pro | mouse_anti | other_pro | other_anti))
message("done")
