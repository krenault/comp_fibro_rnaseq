# Shared publication aesthetics (matched to birds_bats_and_mammals/scripts/shared/plot_aesthetics.R)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 required for plot_aesthetics.R")
}

# Treatment / panel colors (Drive Flexible_homeostasis)
treatment_colors <- c(
  hypoxia    = "#2166AC",
  temperature = "#B2182B",
  glucose    = "#7CB577",
  stress     = "#7671A3",
  cold       = "#6B8E9F",
  heat       = "#CC5454",
  oxidative  = "#D99F6A",
  radiation  = "#5E4B7A",
  desert     = "#C76F84",
  other      = "#3F545E"
)

gradient_low  <- "#2166AC"
gradient_high <- "#B2182B"
title_color   <- "#9b383a"

PUB_BASE <- 16L
PUB_DPI  <- 300L

theme_pub <- function(base_size = PUB_BASE, title_color_override = NULL) {
  tc <- if (is.null(title_color_override)) "grey10" else title_color_override
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = base_size + 4, hjust = 0.5, color = tc
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 2, hjust = 0.5, color = "grey30"
      ),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 2),
      legend.text = ggplot2::element_text(size = base_size - 2),
      legend.title = ggplot2::element_text(size = base_size - 1, face = "bold"),
      plot.caption = ggplot2::element_text(
        size = base_size - 4, hjust = 1, color = "grey40"
      ),
      strip.text = ggplot2::element_text(face = "bold", size = base_size - 1),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 15, 10, 10)
    )
}

theme_cmp <- function(bs = 12) {
  ggplot2::theme_bw(base_size = bs) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey30")
    )
}

sig_caption <- "* FDR/p < 0.05, ** < 0.01, *** < 0.001"

get_sig_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

pretty_ours <- function(x) {
  dplyr::recode(
    x,
    hypoxia_6H_vs_0 = "Hypoxia 6H",
    hypoxia_24H_vs_0 = "Hypoxia 24H",
    temperature_32C_vs_37C = "Cold 32C",
    temperature_41C_vs_37C = "Heat 41C",
    glucose_2.5mM_vs_8mM = "Glucose 2.5 mM",
    glucose_30mM_vs_8mM = "Glucose 30 mM",
    .default = x
  )
}

pretty_species <- function(x) {
  stringr::str_replace_all(x, "_", " ") %>%
    stringr::str_to_title() %>%
    stringr::str_replace("Bactrian Camel", "Bactrian camel") %>%
    stringr::str_replace("Dromedary Camel", "Dromedary camel")
}

div_fill <- function(...) {
  ggplot2::scale_fill_gradient2(
    low = gradient_low, mid = "white", high = gradient_high,
    midpoint = 0, limits = c(-1, 1), oob = scales::squish, ...
  )
}
