#!/usr/bin/env Rscript
# Direction-aware harmonic mean p across longevity conditions.
# Same combiner as generic_stress/share/hmp_combine.R (Wilson 2019 via
# harmonicmeanp::p.hmp, equal weights, multilevel = FALSE): BH after the
# harmonic mean, at least 3 tests, at least 80% sign agreement.
#
# The longevity tables have no raw p-values. Each condition contributes a
# two-sided normal tail from its 95% HPD: z = estimate / ((upper-lower)/(2*1.96)).
# A condition enters only when passes_all_filters is TRUE. Missing conditions
# are left out. L is the number of conditions in that analysis.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

root <- "/Users/katiarenault/Documents/GitHub/comp_fibro_rnaseq"
source(file.path(root, "generic_stress/scripts/shared/plot_aesthetics.R"))
source(file.path(root, "generic_stress/share/hmp_combine.R"))
if (!requireNamespace("harmonicmeanp", quietly = TRUE)) stop("harmonicmeanp is required")

out_dir <- file.path(root, "longevity/results/hmp")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

static_files <- c(
  "32C" = "longevity/data/longevity_genes_static_32C.csv",
  "37C" = "longevity/data/longevity_genes_static_37C.csv",
  "41C" = "longevity/data/longevity_genes_static_41C.csv",
  "hypoxia_6h" = "longevity/data/longevity_genes_static_hypoxia_6h.csv",
  "hypoxia_24h" = "longevity/data/longevity_genes_static_hypoxia_24h.csv"
)
dynamic_files <- c(
  "32C_vs_37C" = "longevity/results/dynamic_genes/longevity_genes_dynamic_32C_vs_37C.csv",
  "41C_vs_37C" = "longevity/results/dynamic_genes/longevity_genes_dynamic_41C_vs_37C.csv",
  "hypoxia_6h_vs_0h" = "longevity/results/dynamic_genes/longevity_genes_dynamic_hypoxia_6h_vs_0h.csv",
  "hypoxia_24h_vs_0h" = "longevity/results/dynamic_genes/longevity_genes_dynamic_hypoxia_24h_vs_0h.csv"
)

tail_rows <- function(path, test, effect_col) {
  d <- read_csv(file.path(root, path), show_col_types = FALSE)
  effect <- d[[effect_col]]
  se <- (d$longevity_upper - d$longevity_lower) / (2 * 1.96)
  z <- effect / se
  p <- 2 * pnorm(-abs(z))
  ok <- d$passes_all_filters %in% TRUE & is.finite(effect) & is.finite(se) & se > 0 & is.finite(p)
  tibble(
    gene = d$gene[ok],
    symbol = toupper(d$hgnc_symbol[ok]),
    p = pmax(p[ok], 1e-300),
    effect = effect[ok],
    test = test,
    called_fc11 = d$effecttype_fc1_1[ok] %in% c(
      "higher_in_long_lived", "lower_in_long_lived",
      "greater_increase", "smaller_increase", "greater_decrease", "smaller_decrease",
      "increase_vs_decrease", "decrease_vs_increase",
      "increase_vs_unchanged", "decrease_vs_unchanged",
      "unchanged_vs_increase", "unchanged_vs_decrease"
    )
  )
}

run_hmp <- function(files, effect_col, L) {
  tests <- imap_dfr(files, ~ tail_rows(.x, .y, effect_col))
  tests$symbol[tests$symbol == ""] <- NA_character_
  res <- hmp_combine(tests, id = "gene", p = "p", effect = "effect", test = "test", L = L,
                     min_n = 3L, sign_agree = 0.80, padj_cutoff = 0.05)
  sym <- tests %>%
    filter(!is.na(symbol)) %>%
    distinct(gene, symbol) %>%
    group_by(gene) %>%
    slice(1) %>%
    ungroup()
  calls <- tests %>%
    filter(called_fc11) %>%
    distinct(gene, test) %>%
    count(gene, name = "n_fc11")
  res %>%
    left_join(sym, by = "gene") %>%
    left_join(calls, by = "gene") %>%
    mutate(
      n_fc11 = replace_na(n_fc11, 0L),
      direction = case_when(
        in_set & hmp_dir > 0 ~ "higher",
        in_set & hmp_dir < 0 ~ "lower",
        TRUE ~ "none"
      )
    )
}

static_hmp <- run_hmp(static_files, "longevity_log2diff", L = 5L)
dynamic_hmp <- run_hmp(dynamic_files, "longevity_logFC_difference", L = 4L)

write_csv(static_hmp, file.path(out_dir, "static_hmp.csv"))
write_csv(dynamic_hmp, file.path(out_dir, "dynamic_hmp.csv"))

summarise_hmp <- function(x, analysis) {
  called <- x %>% filter(in_set)
  tibble(
    analysis = analysis,
    n_genes = nrow(x),
    n_hmp = nrow(called),
    n_higher = sum(called$direction == "higher"),
    n_lower = sum(called$direction == "lower"),
    n_unanimous = sum(called$frac_agree == 1),
    n_also_fc11_once = sum(called$n_fc11 >= 1),
    n_also_fc11_twice = sum(called$n_fc11 >= 2),
    median_abs_effect = median(abs(called$effect_mean))
  )
}

summary_tbl <- bind_rows(
  summarise_hmp(static_hmp, "static_5_conditions"),
  summarise_hmp(dynamic_hmp, "dynamic_4_contrasts")
)
write_csv(summary_tbl, file.path(out_dir, "hmp_summary.csv"))

# fc1.1 genes that the HMP does not keep
static_fc <- static_hmp %>% filter(n_fc11 >= 2)
message("Static HMP")
print(summary_tbl %>% filter(analysis == "static_5_conditions"))
message("fc1.1 in >=2 conditions: ", nrow(static_fc),
        "; of those in HMP: ", sum(static_fc$in_set))
message("Dynamic HMP")
print(summary_tbl %>% filter(analysis == "dynamic_4_contrasts"))
message("Dynamic fc1.1 in >=2 contrasts: ", sum(dynamic_hmp$n_fc11 >= 2),
        "; of those in HMP: ", sum(dynamic_hmp$n_fc11 >= 2 & dynamic_hmp$in_set))

activation <- c("COL1A1", "COL1A2", "POSTN", "ACTA2", "CCN2", "LOXL2", "ITGA11")
message("\nActivation genes, static HMP:")
print(static_hmp %>% filter(symbol %in% activation) %>%
        dplyr::select(symbol, hmp_p, hmp_padj, hmp_dir, n_test, frac_agree, effect_mean, in_set, n_fc11))

# Plot: the static genes that pass, labeled
plot_df <- static_hmp %>%
  filter(in_set) %>%
  mutate(
    direction = factor(direction, levels = c("higher", "lower")),
    symbol = reorder(symbol, effect_mean)
  )
p <- ggplot(plot_df, aes(effect_mean, symbol, color = direction)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_segment(aes(x = 0, xend = effect_mean, yend = symbol), linewidth = 0.6) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c(higher = unname(treatment_colors["heat"]),
                                lower = unname(treatment_colors["hypoxia"])),
                     labels = c(higher = "Higher at higher LQ", lower = "Lower at higher LQ"),
                     name = NULL) +
  labs(
    title = "Cross-condition harmonic mean",
    subtitle = "Seven genes pass. Sign agrees in every condition that was tested.",
    x = "Mean longevity log2 difference",
    y = NULL,
    caption = "BH-adjusted harmonic mean p < 0.05, at least 3 of 5 conditions, at least 80% sign agreement."
  ) +
  theme_pub(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "static_hmp_effects.png"), p, width = 7.2, height = 4.8, dpi = PUB_DPI, bg = "white")
ggsave(file.path(out_dir, "static_hmp_effects.pdf"), p, width = 7.2, height = 4.8, bg = "white")
message("Wrote static_hmp_effects")
message("done")
