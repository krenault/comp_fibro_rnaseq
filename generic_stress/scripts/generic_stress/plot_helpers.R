# Plot output for the generic-stress figures (no numbered locked/ subdir).
PLOT_DIR <- file.path(ROOT, "plots/mammalian_generic_stress")
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

save_png <- function(p, file, w = 8, h = 5.4) {
  ggplot2::ggsave(file.path(PLOT_DIR, file), p, width = w, height = h,
                  dpi = PUB_DPI, bg = "white")
  message("Wrote ", file)
}

CLASS_COLS <- c(
  heat       = unname(treatment_colors["heat"]),
  hypoxia    = unname(treatment_colors["hypoxia"]),
  oxidative  = unname(treatment_colors["oxidative"]),
  dna_damage = unname(treatment_colors["radiation"])
)
CLASS_LABS <- c(
  heat = "Heat", hypoxia = "Hypoxia",
  oxidative = "H2O2", dna_damage = "UV / IR"
)

pretty_exp <- function(x) {
  dplyr::recode(
    x,
    hypoxia_6H = "Hypoxia 6H", hypoxia_24H = "Hypoxia 24H",
    cold = "Cold 32C", heat = "Heat 41C", unmatched = "Unmatched",
    .default = x
  )
}

# Public classes that enter the HMP gene set for each Drive contrast.
# Heat holds out public heat. Hypoxia holds out public hypoxia.
# Cold is not in ASTRA, so it is scored on the full four-class list (includes heat).
holdout_included <- function(exposure) {
  switch(
    as.character(exposure),
    heat         = c("hypoxia", "oxidative", "dna_damage"),
    hypoxia_6H   = c("heat", "oxidative", "dna_damage"),
    hypoxia_24H  = c("heat", "oxidative", "dna_damage"),
    cold         = c("heat", "hypoxia", "oxidative", "dna_damage"),
    unmatched    = c("oxidative", "dna_damage"),
    stop("unknown exposure: ", exposure)
  )
}

holdout_gene_set_label <- function(exposure) {
  switch(
    as.character(exposure),
    heat        = "Drive heat vs public hypoxia + H2O2 + UV",
    hypoxia_6H  = "Drive hypoxia vs public heat + H2O2 + UV",
    hypoxia_24H = "Drive hypoxia vs public heat + H2O2 + UV",
    cold        = "Drive cold vs public heat + hypoxia + H2O2 + UV",
    unmatched   = "Drive vs public H2O2 + UV only",
    as.character(exposure)
  )
}

class_composition_strip <- function(exposures, included_fn = holdout_included) {
  exposures <- as.character(unique(exposures))
  grid <- tidyr::expand_grid(
    exposure = exposures,
    class = names(CLASS_LABS)
  ) %>%
    dplyr::mutate(
      on = purrr::map2_lgl(exposure, class, ~ .y %in% included_fn(.x)),
      exp_lab = factor(pretty_exp(exposure), levels = pretty_exp(exposures)),
      class_lab = factor(CLASS_LABS[class], levels = rev(unname(CLASS_LABS))),
      fill = dplyr::if_else(on, CLASS_COLS[class], "#F0F0F0"),
      lab = dplyr::if_else(on, "in", "out")
    )
  ggplot2::ggplot(grid, ggplot2::aes(exp_lab, class_lab, fill = fill)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = lab, color = on),
                       fontface = "bold", size = 2.7) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_manual(values = c(`TRUE` = "grey15", `FALSE` = "grey55"),
                                guide = "none") +
    ggplot2::labs(x = NULL, y = "Public classes\nin HMP set") +
    theme_pub(11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 15, 0, 10)
    )
}

stack_strip <- function(strip, panel, heights = c(1.05, 4.4)) {
  strip + panel + patchwork::plot_layout(heights = heights, ncol = 1)
}
