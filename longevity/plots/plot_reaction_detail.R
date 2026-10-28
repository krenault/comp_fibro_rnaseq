# Named reaction figures. Exchange reactions, a single metabolite opening to nothing, are omitted.
# rFP is relative flux potential. Positive means higher potential at the higher-longevity anchor.

if (!exists("root")) {
  root <- "/Users/katiarenault/Documents/GitHub/comp_fibro_rnaseq"
  source(file.path(root, "generic_stress/scripts/shared/plot_aesthetics.R"))
  fig <- file.path(root, "longevity/plots/supplement")
  dir.create(fig, recursive = TRUE, showWarnings = FALSE)
  save_fig <- function(p, file, w, h) {
    ggsave(file.path(fig, paste0(file, ".png")), p, width = w, height = h, dpi = PUB_DPI, bg = "white")
    ggsave(file.path(fig, paste0(file, ".pdf")), p, width = w, height = h, bg = "white")
    message("Wrote ", file)
  }
  static_labs <- c(
    "32C" = "32C", "37C" = "37C", "41C" = "41C",
    "hypoxia_6h" = "Hypoxia 6H", "hypoxia_24h" = "Hypoxia 24H"
  )
  dyn_labs <- c(
    "32C_vs_37C" = "Cold 32C", "41C_vs_37C" = "Heat 41C",
    "hypoxia_6h_vs_0h" = "Hypoxia 6H", "hypoxia_24h_vs_0h" = "Hypoxia 24H"
  )
}

suppressPackageStartupMessages(library(tidyverse))

gem <- read_tsv(
  file.path(root, "longevity/data/external/human_gem/reactions.tsv"),
  show_col_types = FALSE
)

cofactors <- c(
  "H+", "H2O", "NAD+", "NADH", "NADP+", "NADPH", "ATP", "ADP", "AMP", "Pi", "PPi",
  "CoA", "FAD", "FADH2", "O2", "SAM", "SAH", "NH3", "CO2", "H2O2", "GTP", "GDP",
  "ubiquinone", "ubiquinol"
)
comp_lab <- c(
  c = "cytosol", m = "mito", x = "peroxisome", r = "ER", n = "nucleus",
  g = "golgi", l = "lysosome", e = "extra", s = "extra", i = "IMS"
)

reaction_name <- function(equation) {
  equation <- ifelse(is.na(equation), "", equation)
  vapply(equation, function(eq) {
    if (!nzchar(eq)) return("unmapped")
    rev <- grepl("<=>", eq, fixed = TRUE)
    sides <- strsplit(eq, "<=>|=>", perl = TRUE)[[1]]
    parts <- lapply(sides, function(side) {
      bits <- trimws(strsplit(side, " + ", fixed = TRUE)[[1]])
      bits <- bits[nzchar(bits)]
      met <- sub("^\\d+\\s+", "", bits)
      comp <- sub("^.*\\[([a-z])\\]$", "\\1", met)
      comp[!grepl("\\[[a-z]\\]$", met)] <- NA_character_
      met <- sub("\\[[a-z]\\]$", "", met)
      keep <- !(met %in% cofactors)
      if (!any(keep) && length(met)) keep <- TRUE
      list(met = met[keep], comp = comp[keep], cof = met[!keep & met %in% cofactors])
    })
    left <- parts[[1]]$met
    right <- if (length(parts) > 1) parts[[2]]$met else character()
    if (!length(left) && !length(right)) return("cofactor only")
    arrow <- if (rev) " <-> " else " to "
    nm <- paste0(paste(left, collapse = " + "), arrow, paste(right, collapse = " + "))
    comps <- unique(na.omit(unlist(lapply(parts, `[[`, "comp"))))
    if (length(comps) == 1 && comps %in% names(comp_lab)) {
      nm <- paste0(nm, ", ", comp_lab[[comps]])
    }
    cofs <- unique(unlist(lapply(parts, `[[`, "cof")))
    tag <- c(
      if (any(c("NAD+", "NADH") %in% cofs)) "NAD",
      if (any(c("NADP+", "NADPH") %in% cofs)) "NADP",
      if (any(c("FAD", "FADH2") %in% cofs)) "FAD"
    )
    if (length(tag)) nm <- paste0(nm, ", ", paste(tag, collapse = "/"))
    if (nchar(nm) > 62) nm <- paste0(substr(nm, 1, 59), "...")
    nm
  }, character(1))
}

is_exchange <- function(equation) {
  !is.na(equation) & grepl("^1 [^\\+]+\\[[a-z]\\] <=>\\s*$", equation)
}

name_table <- gem %>%
  distinct(reaction, .keep_all = TRUE) %>%
  mutate(
    exchange = is_exchange(equation),
    rxn = reaction_name(equation),
    rxn = ifelse(rxn == "unmapped", reaction, rxn)
  )
name_table$rxn <- make.unique(name_table$rxn, sep = " ")

annotate_names <- function(df) {
  df %>% left_join(name_table %>% select(reaction, equation, exchange, rxn), by = "reaction") %>%
    mutate(
      exchange = coalesce(exchange, FALSE),
      rxn = coalesce(rxn, reaction)
    )
}

static_files <- c(
  "32C" = "longevity/results/reactions_static/longevity_reactions_static_32C.csv",
  "37C" = "longevity/results/reactions_static/longevity_reactions_static_37C.csv",
  "41C" = "longevity/results/reactions_static/longevity_reactions_static_41C.csv",
  "hypoxia_6h" = "longevity/results/reactions_static/longevity_reactions_static_hypoxia_6h.csv",
  "hypoxia_24h" = "longevity/results/reactions_static/longevity_reactions_static_hypoxia_24h.csv"
)
static_rx <- imap_dfr(static_files, function(path, cond) {
  read_csv(file.path(root, path), show_col_types = FALSE) %>%
    transmute(
      reaction, condition = cond,
      diff = longevity_rfp_difference,
      called = effecttype_drfp0_05 %in% c("higher_in_long_lived", "lower_in_long_lived"),
      call = effecttype_drfp0_05
    )
}) %>% annotate_names()

kept <- static_rx %>%
  filter(!exchange, called) %>%
  count(reaction, rxn) %>%
  filter(n >= 3)
shared <- static_rx %>%
  filter(reaction %in% kept$reaction) %>%
  mutate(
    condition_lab = factor(static_labs[condition], levels = unname(static_labs)),
    group = case_when(
      grepl("thymine|uracil|dihydro|ureido", rxn, ignore.case = TRUE) ~ "Pyrimidine",
      grepl("retinol|retinal|retin", rxn, ignore.case = TRUE) ~ "Retinoid",
      grepl("sarcosine|dimethylglycine|formaldehyde", rxn, ignore.case = TRUE) ~ "One carbon",
      grepl("cholest|calci|vitamin", rxn, ignore.case = TRUE) ~ "Sterol",
      grepl("histone|lysine", rxn, ignore.case = TRUE) ~ "Methylation",
      TRUE ~ "Other"
    ),
    group = factor(group, levels = c(
      "Pyrimidine", "Retinoid", "One carbon", "Sterol", "Methylation", "Other"
    ))
  )
shared <- shared %>%
  mutate(rxn = forcats::fct_reorder(rxn, diff, .fun = median, .desc = TRUE))
message("Shared internal reactions called in >=3 conditions: ", n_distinct(shared$reaction))

p_shared <- ggplot(shared, aes(condition_lab, rxn, fill = diff)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_point(data = filter(shared, called), size = 1.1, color = "grey10") +
  scale_fill_gradient2(
    low = gradient_low, mid = "white", high = gradient_high,
    midpoint = 0, limits = c(-0.8, 0.8), oob = scales::squish,
    name = "rFP difference"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Reactions called in three or more conditions",
    subtitle = "Exchange reactions removed. Dot: drfp 0.05 call. Positive: higher flux potential at the higher-longevity anchor.",
    x = NULL, y = NULL
  ) +
  theme_pub(11) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    axis.text.y = element_text(size = 7),
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom",
    panel.grid = element_blank()
  )
save_fig(p_shared, "reactions_shared", 13.2, 0.26 * n_distinct(shared$reaction) + 2.6)

dyn_files <- c(
  "32C_vs_37C" = "longevity/results/reactions_dynamic/longevity_reactions_dynamic_32C_vs_37C.csv",
  "41C_vs_37C" = "longevity/results/reactions_dynamic/longevity_reactions_dynamic_41C_vs_37C.csv",
  "hypoxia_6h_vs_0h" = "longevity/results/reactions_dynamic/longevity_reactions_dynamic_hypoxia_6h_vs_0h.csv",
  "hypoxia_24h_vs_0h" = "longevity/results/reactions_dynamic/longevity_reactions_dynamic_hypoxia_24h_vs_0h.csv"
)
dyn_rx <- imap_dfr(dyn_files, function(path, cond) {
  read_csv(file.path(root, path), show_col_types = FALSE) %>%
    filter(!is.na(effecttype_drfp0_05), effecttype_drfp0_05 != "") %>%
    transmute(
      reaction, condition = cond,
      diff = longevity_delta_rfp_difference,
      call = effecttype_drfp0_05
    )
}) %>% annotate_names() %>%
  filter(!exchange) %>%
  mutate(condition_lab = factor(dyn_labs[condition], levels = unname(dyn_labs)))

response_panel <- function(df, title, subtitle) {
  df <- df %>% mutate(rxn = forcats::fct_reorder(rxn, diff))
  ggplot(df, aes(diff, rxn, fill = diff)) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
    geom_col(width = 0.72, color = "grey25", linewidth = 0.3) +
    scale_fill_gradient2(
      low = gradient_low, mid = "white", high = gradient_high,
      midpoint = 0, limits = c(-0.6, 0.6), oob = scales::squish, guide = "none"
    ) +
    facet_grid(condition_lab ~ ., scales = "free_y", space = "free_y") +
    labs(title = title, subtitle = subtitle, x = "Response difference in rFP", y = NULL) +
    theme_pub(11) +
    theme(axis.text.y = element_text(size = 7), strip.text.y = element_text(angle = 0))
}

dyn_small <- dyn_rx %>% filter(condition != "hypoxia_24h_vs_0h")
dyn_24 <- dyn_rx %>% filter(condition == "hypoxia_24h_vs_0h")
message("Internal response calls: ", paste(names(table(dyn_rx$condition_lab)), table(dyn_rx$condition_lab), sep = " ", collapse = "; "))

p_dyn <- response_panel(
  dyn_small,
  "Reaction responses at cold, heat, and 6 h",
  "Exchange reactions removed. Positive: the higher-longevity anchor has the more positive change in flux potential."
)
save_fig(p_dyn, "reaction_responses", 12.4, 0.32 * nrow(dyn_small) + 1.8)

p_dyn24 <- response_panel(
  dyn_24,
  "Reaction responses at 24 h hypoxia",
  "Exchange reactions removed. Positive: the higher-longevity anchor has the more positive change in flux potential."
)
save_fig(p_dyn24, "reaction_responses_24h", 12.4, 0.22 * nrow(dyn_24) + 1.8)

focus_ids <- c(
  "MAR04344", "MAR04345", "MAR04346",
  "MAR04470", "MAR04471", "MAR04472", "MAR04736",
  "MAR06630", "MAR06631", "MAR06632", "MAR06633", "MAR06673"
)
focus_lab <- c(
  MAR04344 = "Uracil / dihydrouracil, NAD, cytosol",
  MAR04345 = "Uracil / dihydrouracil, NADP, cytosol",
  MAR04346 = "Dihydrouracil to ureidopropionate",
  MAR04470 = "Thymine to dihydrothymine, NADH, cytosol",
  MAR04471 = "Thymine to dihydrothymine, NADH, mito",
  MAR04472 = "Thymine to dihydrothymine, NADPH, cytosol",
  MAR04736 = "Dihydrothymine to ureidoisobutyrate",
  MAR06630 = "Retinol to retinal, NADP, cytosol",
  MAR06631 = "Retinol to retinal, NAD, cytosol",
  MAR06632 = "Retinol to retinal, NADP, ER",
  MAR06633 = "Retinol to retinal, NAD, ER",
  MAR06673 = "13-cis-retinal <-> retinal, ER"
)
focus <- static_rx %>%
  filter(reaction %in% focus_ids) %>%
  mutate(
    condition_lab = factor(static_labs[condition], levels = unname(static_labs)),
    rxn = factor(focus_lab[reaction], levels = rev(unname(focus_lab)))
  )
p_focus <- ggplot(focus, aes(condition_lab, rxn, fill = diff)) +
  geom_tile(aes(color = called, linewidth = called)) +
  geom_text(aes(label = sprintf("%+.2f", diff)), size = 2.8, color = "grey10") +
  scale_color_manual(values = c(`TRUE` = "grey10", `FALSE` = "white"), guide = "none") +
  scale_linewidth_manual(values = c(`TRUE` = 1.15, `FALSE` = 0.5), guide = "none") +
  scale_fill_gradient2(
    low = gradient_low, mid = "white", high = gradient_high,
    midpoint = 0, limits = c(-0.8, 0.8), oob = scales::squish,
    name = "rFP difference"
  ) +
  labs(
    title = "Pyrimidine reduction and retinol oxidation",
    subtitle = "Black border: drfp 0.05 call. Positive: higher flux potential at the higher-longevity anchor.",
    x = NULL, y = NULL
  ) +
  theme_pub(12) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1),
    legend.position = "bottom",
    panel.grid = element_blank()
  )
save_fig(p_focus, "pyrimidine_retinol", 8.6, 6.4)
