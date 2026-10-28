#!/usr/bin/env Rscript
# Activated-fibroblast state: longevity ranking, mammal cross-tissue
# PGLS, and a species score with Pagel's lambda and leave-one-species-out.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ape)
  library(caper)
  library(nlme)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

root <- "/Users/katiarenault/Documents/GitHub/comp_fibro_rnaseq"
source(file.path(root, "generic_stress/scripts/shared/plot_aesthetics.R"))

plot_dir <- file.path(root, "longevity/results/static_genes/plots")
out_dir <- file.path(root, "longevity/results/static_genes")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, file, w, h) {
  ggsave(file.path(plot_dir, paste0(file, ".png")), p,
         width = w, height = h, dpi = PUB_DPI, bg = "white")
  ggsave(file.path(plot_dir, paste0(file, ".pdf")), p,
         width = w, height = h, bg = "white")
  message("Wrote ", file)
}

# Collagen I, periostin, ACTA2, CCN2, integrin a11.
# LOX is absent from the ortholog tables; LOXL2 is the crosslinking gene present.
activation_genes <- c("COL1A1", "COL1A2", "POSTN", "ACTA2", "CCN2", "LOXL2", "ITGA11")

static_conds <- c("32C", "37C", "41C", "hypoxia_6h", "hypoxia_24h")
static_labs <- c(
  "32C" = "32C", "37C" = "37C", "41C" = "41C",
  "hypoxia_6h" = "Hypoxia 6H", "hypoxia_24h" = "Hypoxia 24H"
)
cond_colors <- c(
  "32C" = unname(treatment_colors["cold"]),
  "37C" = "#3F545E",
  "41C" = unname(treatment_colors["heat"]),
  "Hypoxia 6H" = unname(treatment_colors["hypoxia"]),
  "Hypoxia 24H" = "#67A9CF"
)

read_gmt_set <- function(path, set_name) {
  lines <- readLines(path)
  hit <- lines[startsWith(lines, paste0(set_name, "\t"))]
  if (length(hit) != 1) stop("GMT set not found: ", set_name)
  parts <- strsplit(hit, "\t", fixed = TRUE)[[1]]
  toupper(parts[-c(1, 2)])
}

gmt <- "/Users/katiarenault/PhD/GSEA/c2.cgp.v2023.2.Hs.symbols.gmt"
bm_genes <- read_gmt_set(gmt, "NABA_BASEMENT_MEMBRANES")
col_genes <- read_gmt_set(gmt, "NABA_COLLAGENS")
sec_genes <- read_gmt_set(gmt, "NABA_SECRETED_FACTORS")

read_static <- function(cond) {
  read_csv(
    file.path(root, "longevity/data", paste0("longevity_genes_static_", cond, ".csv")),
    show_col_types = FALSE
  ) %>%
    transmute(
      gene,
      symbol = toupper(hgnc_symbol),
      diff = longevity_log2diff,
      called = effecttype_fc1_1 %in% c("higher_in_long_lived", "lower_in_long_lived"),
      condition = cond
    ) %>%
    filter(symbol != "", is.finite(diff))
}

static <- map_dfr(static_conds, read_static)
symbol_to_gene <- static %>%
  distinct(symbol, gene) %>%
  group_by(symbol) %>%
  slice(1) %>%
  ungroup()

set_test <- function(df, genes) {
  x <- df$diff[df$symbol %in% genes]
  y <- df$diff[!df$symbol %in% genes]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  wt_low <- wilcox.test(x, y, alternative = "less", exact = FALSE)
  wt_two <- wilcox.test(x, y, alternative = "two.sided", exact = FALSE)
  tibble(
    n = length(x),
    median = median(x),
    mean = mean(x),
    median_background = median(y),
    p_lower = unname(wt_low$p.value),
    p_two_sided = unname(wt_two$p.value)
  )
}

set_defs <- list(
  activation = activation_genes,
  basement_membrane = bm_genes,
  collagens = col_genes,
  secreted_factors = sec_genes
)

set_stats <- imap_dfr(set_defs, function(genes, set_name) {
  map_dfr(static_conds, function(cond) {
    set_test(filter(static, condition == cond), genes) %>%
      mutate(condition = cond, geneset = set_name)
  })
}) %>%
  group_by(geneset) %>%
  mutate(p_lower_bh = p.adjust(p_lower, method = "BH")) %>%
  ungroup() %>%
  mutate(
    condition_lab = static_labs[condition],
    star = get_sig_star(p_lower_bh)
  )

gene_stats <- static %>%
  filter(symbol %in% activation_genes) %>%
  mutate(
    condition_lab = factor(static_labs[condition], levels = unname(static_labs)),
    symbol = factor(symbol, levels = rev(activation_genes)),
    label = paste0(sprintf("%+.1f", diff), ifelse(called, "*", ""))
  )

write_csv(gene_stats %>% dplyr::select(condition, symbol, diff, called),
          file.path(out_dir, "fibroblast_activation_genes.csv"))
write_csv(set_stats, file.path(out_dir, "fibroblast_set_wilcoxon.csv"))

# --- Mammal cross-tissue PGLS (BBM), positive beta = higher with MLres ---
bbm <- read_csv(
  "/Users/katiarenault/Documents/GitHub/birds_bats_and_mammals/results/mammals/cross_tissue_pgls_meta_nobatch/gene_meta_symbol_MLres.csv",
  show_col_types = FALSE
) %>%
  transmute(
    symbol = toupper(Symbol),
    beta = beta_IVW,
    fdr = FDR_IVW,
    hit = hit_strict == 1
  ) %>%
  filter(symbol != "", is.finite(beta))

bbm_set_test <- function(genes, set_name) {
  inn <- bbm$symbol %in% genes
  x <- bbm$beta[inn]
  y <- bbm$beta[!inn]
  wt_low <- wilcox.test(x, y, alternative = "less", exact = FALSE)
  wt_two <- wilcox.test(x, y, alternative = "two.sided", exact = FALSE)
  a <- sum(inn & bbm$hit)
  b <- sum(inn & !bbm$hit)
  c <- sum(!inn & bbm$hit)
  d <- sum(!inn & !bbm$hit)
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  tibble(
    geneset = set_name,
    n = length(x),
    n_hit = a,
    median_beta = median(x),
    mean_beta = mean(x),
    frac_negative = mean(x < 0),
    p_lower = unname(wt_low$p.value),
    p_two_sided = unname(wt_two$p.value),
    fisher_or = unname(ft$estimate),
    fisher_p = unname(ft$p.value)
  )
}

bbm_sets <- map2_dfr(set_defs, names(set_defs), bbm_set_test) %>%
  mutate(p_lower_bh = p.adjust(p_lower, method = "BH"))
write_csv(bbm_sets, file.path(out_dir, "fibroblast_vs_bbm_sets.csv"))

our37 <- static %>%
  filter(condition == "37C") %>%
  dplyr::select(symbol, diff37 = diff)

bbm_genes <- bbm %>%
  inner_join(our37, by = "symbol") %>%
  mutate(
    in_activation = symbol %in% activation_genes,
    in_bm = symbol %in% bm_genes
  ) %>%
  filter(in_activation | in_bm)

write_csv(
  bbm %>% filter(symbol %in% unique(unlist(set_defs))),
  file.path(out_dir, "fibroblast_vs_bbm_genes.csv")
)

spear_bm <- bbm_genes %>%
  filter(in_bm) %>%
  summarise(
    n = n(),
    rho = cor(diff37, beta, method = "spearman"),
    p = cor.test(diff37, beta, method = "spearman", exact = FALSE)$p.value
  )

# --- Species means at 37C / 8 mM / normoxia ---
tpm_path <- "/Users/katiarenault/Desktop/PhD/Flexible_homeostasis/STAR_featurecounts/TPM_Matthew/TPM_13-species.one2one.csv"
header <- names(fread(tpm_path, nrows = 0))
base_cols <- header[grepl("_37C_8mM_0$", header)]
species_order <- c(
  "honey_badger", "little_brown_bat", "dromedary_camel", "bactrian_camel",
  "rousette", "dolphin", "whale", "gelada", "human", "squirrel", "rhino",
  "seal", "rat"
)
species_order <- species_order[order(nchar(species_order), decreasing = TRUE)]
sp_of <- function(col) {
  hit <- species_order[startsWith(col, paste0(species_order, "_"))]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}
base_sp <- setNames(vapply(base_cols, sp_of, character(1)), base_cols)

score_symbols <- unique(c(activation_genes, bm_genes))
score_map <- symbol_to_gene %>% filter(symbol %in% score_symbols)
tpm <- fread(tpm_path, select = c("GeneID", base_cols))
tpm <- tpm[GeneID %in% score_map$gene]
tpm_long <- melt(tpm, id.vars = "GeneID", variable.name = "sample", value.name = "tpm")
tpm_long[, species := base_sp[as.character(sample)]]
tpm_long[, log2tpm := log2(tpm + 1)]
species_mean <- tpm_long[, .(log2tpm = mean(log2tpm, na.rm = TRUE), n_indiv = uniqueN(sample)),
                         by = .(GeneID, species)]
species_mean <- merge(
  species_mean,
  as.data.table(score_map),
  by.x = "GeneID", by.y = "gene",
  all.x = TRUE
)

pheno_raw <- read_csv(
  "/Users/katiarenault/Desktop/PhD/Flexible_homeostasis/Lucas_method/phenotype_metadata_additional.csv",
  show_col_types = FALSE
) %>%
  dplyr::select(1:3)
names(pheno_raw) <- c("species_label", "ml", "mass")
pheno_key <- c(
  "whale" = "whale", "Seal" = "seal", "Squirrel" = "squirrel",
  "Bactrian Camel" = "bactrian_camel", "Dolphin" = "dolphin",
  "Dromedary Camel" = "dromedary_camel", "Gelada" = "gelada",
  "Honey Badger" = "honey_badger", "Human" = "human",
  "Little Brown Bat" = "little_brown_bat", "Rat" = "rat", "Rhino" = "rhino"
)
pheno <- pheno_raw %>%
  mutate(
    species = unname(pheno_key[species_label]),
    ml = as.numeric(ml),
    mass = as.numeric(str_replace_all(as.character(mass), ",", ""))
  ) %>%
  filter(!is.na(species))

tree <- read.tree("/Users/katiarenault/Desktop/PhD/Flexible_homeostasis/Lucas_method/species_tree.nwk")
if (!is.rooted(tree)) tree <- root(tree, outgroup = "little_brown_bat", resolve.root = TRUE)

panel_species <- intersect(tree$tip.label, pheno$species)

signature_score <- function(symbols, species_keep) {
  sub <- species_mean[symbol %in% symbols & species %in% species_keep]
  wide <- data.table::dcast(sub, species ~ symbol, value.var = "log2tpm")
  sp <- wide$species
  gene_cols <- setdiff(names(wide), "species")
  z <- scale(as.matrix(wide[, gene_cols, with = FALSE]))
  score <- rowMeans(z, na.rm = TRUE)
  tibble(species = sp, score = as.numeric(score), n_genes = length(gene_cols))
}

add_resid <- function(df) {
  meta <- pheno %>% filter(species %in% df$species)
  fit <- lm(log10(ml) ~ log10(mass), data = meta)
  meta$ml_resid <- unname(resid(fit))
  meta$log10_ml <- log10(meta$ml)
  meta$log10_mass <- log10(meta$mass)
  df %>% left_join(meta %>% dplyr::select(species, ml, mass, ml_resid, log10_ml, log10_mass), by = "species")
}

fit_row <- function(n, lambda, intercept, slope, se, t, p, model) {
  tibble(n = n, lambda = lambda, intercept = intercept, slope = slope,
         se = se, t = t, p = p, model = model)
}

fit_pgls <- function(df) {
  df <- as.data.frame(df)
  df <- df[is.finite(df$score) & is.finite(df$ml_resid), ]
  df$species <- as.character(df$species)
  phy <- keep.tip(tree, df$species)
  df <- df[match(phy$tip.label, df$species), ]
  rownames(df) <- df$species
  cd <- comparative.data(phy, df, species, vcv = TRUE, na.omit = FALSE, warn.dropped = FALSE)

  caper_row <- function(lambda_arg, model_name) {
    fit <- pgls(score ~ ml_resid, data = cd, lambda = lambda_arg,
                bounds = list(lambda = c(1e-6, 1)))
    coefs <- summary(fit)$coefficients
    cn <- colnames(coefs)
    se_col <- cn[grepl("Std", cn)][1]
    t_col <- cn[grepl("^t", cn)][1]
    p_col <- cn[grepl("Pr", cn)][1]
    lam <- if (is.numeric(lambda_arg)) lambda_arg else unname(fit$param["lambda"])
    fit_row(nrow(df), lam,
            unname(coefs["(Intercept)", "Estimate"]),
            unname(coefs["ml_resid", "Estimate"]),
            unname(coefs["ml_resid", se_col]),
            unname(coefs["ml_resid", t_col]),
            unname(coefs["ml_resid", p_col]),
            model_name)
  }

  ml_try <- try(caper_row("ML", "pgls_lambda_ML"), silent = TRUE)
  if (!inherits(ml_try, "try-error")) return(ml_try)

  df$sp <- factor(df$species, levels = phy$tip.label)
  nlme_try <- try(gls(
    score ~ ml_resid, data = df,
    correlation = corPagel(0.6, phy = phy, form = ~ sp),
    method = "ML",
    control = glsControl(msMaxIter = 80, tolerance = 1e-3, singular.ok = TRUE)
  ), silent = TRUE)
  if (!inherits(nlme_try, "try-error")) {
    tt <- summary(nlme_try)$tTable
    lam <- unname(coef(nlme_try$modelStruct$corStruct, unconstrained = FALSE))
    return(fit_row(
      nrow(df), lam, tt["(Intercept)", "Value"], tt["ml_resid", "Value"],
      tt["ml_resid", "Std.Error"], tt["ml_resid", "t-value"],
      tt["ml_resid", "p-value"], "gls_lambda_ML"
    ))
  }

  br_try <- try(caper_row(1, "pgls_brownian"), silent = TRUE)
  if (!inherits(br_try, "try-error")) return(br_try)
  caper_row(0, "pgls_lambda0")
}

run_score_models <- function(symbols, score_name) {
  full <- signature_score(symbols, panel_species) %>% add_resid()
  full_fit <- fit_pgls(full) %>% mutate(dropped = "none", score_name = score_name)
  loo <- map_dfr(panel_species, function(sp) {
    kept <- setdiff(panel_species, sp)
    signature_score(symbols, kept) %>%
      add_resid() %>%
      fit_pgls() %>%
      mutate(dropped = sp, score_name = score_name)
  })
  list(points = full %>% mutate(score_name = score_name), fits = bind_rows(full_fit, loo))
}

act_model <- run_score_models(activation_genes, "activation")
bm_present <- intersect(bm_genes, unique(species_mean$symbol))
bm_model <- run_score_models(bm_present, "basement_membrane")

score_points <- bind_rows(act_model$points, bm_model$points)
score_fits <- bind_rows(act_model$fits, bm_model$fits)
write_csv(score_points, file.path(out_dir, "fibroblast_species_scores.csv"))
write_csv(score_fits, file.path(out_dir, "fibroblast_species_pgls.csv"))

indiv_n <- species_mean %>%
  filter(species %in% panel_species) %>%
  distinct(species, n_indiv)
message("Individuals per species at baseline:")
print(indiv_n)

# --- Plots ---
p_genes <- ggplot(gene_stats, aes(condition_lab, symbol, fill = diff)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = label, color = abs(diff) > 2), size = 3.4, show.legend = FALSE) +
  scale_fill_gradient2(
    low = gradient_low, mid = "white", high = gradient_high,
    midpoint = 0, limits = c(-5, 5), oob = squish, name = "Longevity\nlog2 diff"
  ) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey10")) +
  scale_x_discrete(expand = expansion(add = 0.05)) +
  labs(
    title = "Activated-fibroblast genes",
    subtitle = "Positive: higher at the higher-longevity anchor. Star: fc1.1 call.",
    x = NULL, y = NULL,
    caption = "COL1A1, COL1A2, POSTN, ACTA2, CCN2, ITGA11, and LOXL2. LOX is absent."
  ) +
  theme_pub(base_size = 13) +
  theme(legend.position = "right", panel.grid = element_blank())

set_plot_df <- set_stats %>%
  mutate(
    geneset = recode(
      geneset,
      activation = "Activation",
      basement_membrane = "Basement membrane",
      collagens = "Collagens",
      secreted_factors = "Secreted factors"
    ),
    geneset = factor(geneset, levels = c("Activation", "Basement membrane", "Collagens", "Secreted factors")),
    condition_lab = factor(condition_lab, levels = unname(static_labs)),
    label = paste0(sprintf("%+.2f", median), star)
  )

p_sets <- ggplot(set_plot_df, aes(condition_lab, median, fill = condition_lab)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = label, vjust = ifelse(median < 0, 1.2, -0.3)), size = 3.1) +
  facet_wrap(~geneset, nrow = 1) +
  scale_fill_manual(values = cond_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.18))) +
  labs(
    title = "Longevity difference by gene set",
    subtitle = "Median log2 difference. Stars: BH-adjusted one-sided Wilcoxon, set lower than other genes.",
    x = NULL, y = "Median longevity log2 difference"
  ) +
  theme_pub(base_size = 12) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

bbm_plot <- bbm_genes %>%
  mutate(
    group = case_when(
      in_activation ~ "Activation",
      in_bm ~ "Basement membrane",
      TRUE ~ "Other"
    )
  )

p_bbm <- ggplot(bbm_plot, aes(diff37, beta)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
  geom_point(data = filter(bbm_plot, !in_activation), color = "#67A9CF", size = 2.4, alpha = 0.9) +
  geom_point(data = filter(bbm_plot, in_activation), color = unname(treatment_colors["heat"]), size = 2.8) +
  geom_text_repel(
    data = filter(bbm_plot, in_activation),
    aes(label = symbol),
    size = 3.3, color = "grey10", seed = 1, min.segment.length = 0,
    box.padding = 0.35, max.overlaps = 20
  ) +
  labs(
    title = "Same genes in mammal tissues",
    subtitle = sprintf(
      "Basement-membrane genes, Spearman rho = %.2f (n = %d, p = %.3f). Red: activation genes.",
      spear_bm$rho, spear_bm$n, spear_bm$p
    ),
    x = "Fibroblast longevity difference at 37C",
    y = "Mammal cross-tissue PGLS beta",
    caption = "Positive beta: higher expression with residual lifespan across brain, liver, and kidney."
  ) +
  theme_pub(base_size = 13)

score_lab <- c(activation = "Activation program", basement_membrane = "Basement membrane")
full_fits <- score_fits %>% filter(dropped == "none")
point_df <- score_points %>%
  mutate(
    score_name = factor(score_lab[score_name], levels = unname(score_lab)),
    species_lab = pretty_species(species)
  )
line_df <- full_fits %>%
  mutate(score_name = factor(score_lab[score_name], levels = unname(score_lab)))

p_scatter <- ggplot(point_df, aes(ml_resid, score)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey75") +
  geom_abline(data = line_df, aes(intercept = intercept, slope = slope),
              color = unname(treatment_colors["temperature"]), linewidth = 0.7) +
  geom_point(size = 2.6, color = "#3F545E") +
  geom_text_repel(aes(label = species_lab), size = 3.1, seed = 2,
                  box.padding = 0.4, max.overlaps = 20, color = "grey20") +
  facet_wrap(~score_name, scales = "free_y") +
  labs(
    title = "Species score versus residual lifespan",
    subtitle = "Line is the PGLS fit (Pagel's lambda). Score is the mean z-scored log2(TPM+1) at 37C.",
    x = "Residual log10 lifespan | log10 mass",
    y = "Signature score"
  ) +
  theme_pub(base_size = 12)

loo_df <- score_fits %>%
  filter(dropped != "none") %>%
  mutate(
    score_name = factor(score_lab[score_name], levels = unname(score_lab)),
    species_lab = pretty_species(dropped),
    species_lab = reorder(species_lab, slope)
  )
full_seg <- full_fits %>%
  mutate(score_name = factor(score_lab[score_name], levels = unname(score_lab)))

p_loo <- ggplot(loo_df, aes(slope, species_lab)) +
  geom_vline(data = full_seg, aes(xintercept = slope),
             color = unname(treatment_colors["temperature"]), linewidth = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_point(size = 2.4, color = "#3F545E") +
  facet_wrap(~score_name, scales = "free_x") +
  labs(
    title = "Leave-one-species-out slope",
    subtitle = "Each point refits PGLS, z-scores, and the mass residual without that species. Red line: full fit.",
    x = "PGLS slope", y = "Species dropped"
  ) +
  theme_pub(base_size = 12)

p_species <- p_scatter / p_loo + plot_layout(heights = c(1.05, 1))

save_fig(p_genes, "fibroblast_activation_genes", 8.2, 4.6)
save_fig(p_sets, "fibroblast_set_shift", 11.2, 4.8)
save_fig(p_bbm, "fibroblast_vs_mammals", 7.6, 6.2)
save_fig(p_species, "fibroblast_species_pgls", 9.2, 9.4)

message("\nSet medians:")
print(set_stats %>% dplyr::select(geneset, condition, n, median, p_lower, p_lower_bh) %>% arrange(geneset, condition))
message("\nBBM genome-wide median beta: ", median(bbm$beta, na.rm = TRUE))
message("\nBBM sets:")
print(bbm_sets)
message("\nPGLS:")
print(score_fits %>% filter(dropped == "none"))
message("\nLOO slope range:")
print(score_fits %>% filter(dropped != "none") %>% group_by(score_name) %>%
        summarise(min_slope = min(slope), max_slope = max(slope),
                  n_pos = sum(slope > 0), n_p05 = sum(p < 0.05), .groups = "drop"))
message("done")
