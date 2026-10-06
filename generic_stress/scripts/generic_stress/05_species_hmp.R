#!/usr/bin/env Rscript
# 05 — per-species Drive cross-perturbation HMP, conservation, enrichment.
# Same combination rule as the public generic (BH HMP < 0.05, ≥3 arms, ≥80% sign).
# Primary insults: heat, cold, hypoxia (24H preferred). Glucose dropped.

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(fgsea)
})

.file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
.here <- if (length(.file)) dirname(normalizePath(.file[[1]])) else getwd()
source(file.path(.here, "..", "shared", "root.R"))
source("scripts/shared/plot_aesthetics.R")
source("scripts/shared/geo_helpers.R")
source("scripts/shared/hmp_twosided.R")
source("scripts/generic_stress/plot_helpers.R")

RES <- file.path(ROOT, "results/mammalian_generic_stress")
OUT <- file.path(RES, "drive_species_hmp")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

SIGN_AGREE <- 0.80
HMP_PADJ <- 0.05
MIN_N <- 3L
N_SP_CONS <- 5L

csv3 <- file.path(OUT, "drive_species_hmp_three_insult.csv")

if (file.exists(csv3)) {
  message("Loading existing species HMP table...")
  h3 <- fread(csv3)
} else {
  source("scripts/shared/hmp_twosided.R")
  drive_files <- list_drive_dge()
  drive_files <- drive_files[!grepl("glucose", drive_files)]
  parse_one <- function(f) {
    bn <- tools::file_path_sans_ext(basename(f))
    tt <- stringr::str_extract(bn, "hypoxia|temperature")
    cmp <- stringr::str_extract(bn, "(6H_vs_0|24H_vs_0|32C_vs_37C|41C_vs_37C)$")
    if (is.na(tt) || is.na(cmp)) return(NULL)
    sp <- sub(paste0("_", tt, "_", cmp, "$"), "", bn)
    exposure <- dplyr::case_when(
      tt == "hypoxia" & cmp == "6H_vs_0" ~ "hypoxia_6H",
      tt == "hypoxia" & cmp == "24H_vs_0" ~ "hypoxia_24H",
      grepl("32C", cmp) ~ "cold",
      grepl("41C", cmp) ~ "heat",
      TRUE ~ NA_character_
    )
    d <- data.table::fread(f, showProgress = FALSE)
    data.table(
      gene = strip(d$gene), lfc = as.numeric(d$log2FoldChange),
      p = as.numeric(d$pvalue), species = sp, exposure = exposure
    )[is.finite(lfc) & grepl("^ENSG", gene) & !is.na(exposure)]
  }
  drive <- rbindlist(lapply(drive_files, parse_one))
  drive <- drive[species != "rousette"]
  hyp_pref <- data.table::copy(drive[exposure == "hypoxia_24H"])
  hyp_alt <- data.table::copy(
    drive[exposure == "hypoxia_6H" &
            !paste(species, gene) %in% hyp_pref[, paste(species, gene)]]
  )
  three <- rbindlist(list(
    data.table::copy(drive[exposure %in% c("heat", "cold")]),
    hyp_pref[, exposure := "hypoxia"],
    hyp_alt[, exposure := "hypoxia"]
  ))
  hmp_by_species <- function(dt, L) {
    rbindlist(lapply(split(dt, by = "species"), function(spdt) {
      wide <- dcast(spdt, gene ~ exposure, value.var = c("p", "lfc"))
      p_cols <- grep("^p_", names(wide), value = TRUE)
      lfc_cols <- grep("^lfc_", names(wide), value = TRUE)
      rows <- lapply(seq_len(nrow(wide)), function(i) {
        pv <- unlist(wide[i, ..p_cols], use.names = FALSE)
        ev <- unlist(wide[i, ..lfc_cols], use.names = FALSE)
        h <- hmp_twosided(pv, ev, L = L)
        data.table(
          gene = wide$gene[i], hmp_p = h$p, hmp_dir = h$dir,
          n_arm = h$n, frac_agree = h$frac_agree,
          lfc_mean = mean(ev[is.finite(ev)])
        )
      })
      out <- rbindlist(rows)
      out[, species := unique(spdt$species)]
      out[, hmp_padj := p.adjust(hmp_p, method = "BH")]
      out[, in_cross := is.finite(hmp_padj) & hmp_padj < HMP_PADJ &
            n_arm >= MIN_N & frac_agree >= SIGN_AGREE]
      out[, direction := fcase(
        in_cross & hmp_dir > 0, "up",
        in_cross & hmp_dir < 0, "down",
        default = "none"
      )]
      out
    }))
  }
  message("Computing species HMP...")
  h3 <- hmp_by_species(three, L = 3L)
  h3[, panel := "three_insult"]
  fwrite(h3, csv3)
}

pub <- fread(file.path(RES, "gene_human_mouse_generic.csv"))
pub_gen <- pub[in_generic == TRUE, .(gene, pub_dir = direction)]

overlap_one <- function(h) {
  h[in_cross == TRUE, {
    m <- merge(.SD, pub_gen, by = "gene", all.x = TRUE)
    n <- .N
    n_pub <- sum(gene %in% pub_gen$gene)
    .(
      n_cross = n,
      n_up = sum(direction == "up"),
      n_down = sum(direction == "down"),
      n_in_public_generic = n_pub,
      n_same_dir_public = sum(!is.na(m$pub_dir) & m$pub_dir == direction),
      frac_of_cross_in_public = n_pub / n,
      jaccard_public = length(intersect(gene, pub_gen$gene)) /
        length(union(gene, pub_gen$gene))
    )
  }, by = species]
}
sum3 <- overlap_one(h3)
fwrite(sum3, file.path(OUT, "drive_species_hmp_overlap_three_insult.csv"))

sp_lev <- c(
  "human", "gelada", "rat", "squirrel", "little_brown_bat",
  "honey_badger", "seal", "dolphin", "whale", "rhino",
  "dromedary_camel", "bactrian_camel"
)
sum3b <- sum3[species %in% sp_lev]
sum3b[, species_lab := factor(pretty_species(species),
                              levels = pretty_species(sp_lev))]
n_long <- melt(
  sum3b,
  id.vars = "species_lab",
  measure.vars = c("n_up", "n_down"),
  variable.name = "direction", value.name = "n"
)
n_long[, direction := fifelse(direction == "n_up", "Up", "Down")]
p_n <- ggplot(n_long, aes(species_lab, n, fill = direction)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.3) +
  scale_fill_manual(values = c(Up = unname(gradient_high),
                               Down = unname(gradient_low)), name = NULL) +
  labs(
    title = "Drive cross-perturbation genes, per species",
    subtitle = "HMP of heat, cold, and hypoxia. Same rule as the public generic: BH HMP < 0.05, >=3 arms, >=80% sign.",
    x = NULL, y = "Genes"
  ) +
  theme_pub(13) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        legend.position = "bottom")
save_png(p_n, "species_hmp_counts.png", 9.2, 5.6)

p_ov <- ggplot(sum3b, aes(species_lab, frac_of_cross_in_public)) +
  geom_col(width = 0.72, fill = unname(treatment_colors["stress"]),
           color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d / %d\n(%.0f%%)",
                                n_in_public_generic, n_cross,
                                100 * frac_of_cross_in_public)),
            vjust = -0.15, size = 2.9, fontface = "bold", lineheight = 0.95) +
  coord_cartesian(ylim = c(0, 0.28)) +
  labs(
    title = "Most species-cross-perturbation genes are not the public generic",
    subtitle = "Each bar = (genes in the public 1,344) / (genes called by that species' heat+cold+hypoxia HMP). Human 204/3612 is 6%.",
    x = NULL, y = "Public 1,344 / species HMP list"
  ) +
  theme_pub(12) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_png(p_ov, "species_hmp_vs_public.png", 9.4, 5.8)

# -----------------------------------------------------------------------------
# Conservation across Drive species
# -----------------------------------------------------------------------------
x <- h3[in_cross == TRUE]
cons <- x[, .(
  n_species = uniqueN(species),
  n_up = sum(direction == "up"),
  n_down = sum(direction == "down"),
  mean_lfc = mean(lfc_mean, na.rm = TRUE)
), by = gene]
cons[, frac_maj := pmax(n_up, n_down) / n_species]
cons[, direction := fifelse(n_up >= n_down, "up", "down")]
cons[, conserved5 := n_species >= N_SP_CONS & frac_maj >= SIGN_AGREE]
cons[, in_public_generic := gene %in% pub_gen$gene]
cons <- merge(cons, pub_gen, by = "gene", all.x = TRUE)
fwrite(cons, file.path(OUT, "drive_species_hmp_conservation.csv"))

cat("\nConservation (in_cross, 12 Drive species):\n")
print(table(cons$n_species))
for (k in c(3L, 5L, 8L, 10L, 12L)) {
  s <- cons[n_species >= k & frac_maj >= SIGN_AGREE]
  cat(sprintf(
    "  >=%d species, >=80%% same dir: %d  (up %d / down %d); public generic %d\n",
    k, nrow(s), sum(s$direction == "up"), sum(s$direction == "down"),
    sum(s$in_public_generic)
  ))
}

n_sp <- as.integer(uniqueN(x$species))
hist_df <- cons[, .N, by = n_species]
p_cons <- ggplot(hist_df, aes(n_species, N)) +
  geom_col(width = 0.8, fill = unname(treatment_colors["stress"]),
           color = "grey25", linewidth = 0.3) +
  geom_text(aes(label = N), vjust = -0.3, fontface = "bold", size = 3.2) +
  scale_x_continuous(breaks = 1:n_sp) +
  labs(
    title = "Most Drive cross-perturbation genes are species-restricted",
    subtitle = sprintf(
      "%d genes in >=5 species with >=80%% direction. None in all %d species. %d of the >=5 set are in the public 1,344.",
      sum(cons$conserved5), n_sp, sum(cons$conserved5 & cons$in_public_generic)
    ),
    x = "Drive species calling the gene", y = "Genes"
  ) +
  theme_pub(13)
save_png(p_cons, "species_hmp_conservation.png", 8.4, 5.4)

# -----------------------------------------------------------------------------
# ORA of the conserved ≥5 set (enough genes for Hallmark / Reactome / GO BP)
# -----------------------------------------------------------------------------
msig_one <- function(collection, sub = NULL) {
  tryCatch({
    if (is.null(sub)) msigdbr::msigdbr(species = "Homo sapiens", collection = collection)
    else msigdbr::msigdbr(species = "Homo sapiens", collection = collection, subcollection = sub)
  }, error = function(e) {
    if (is.null(sub)) msigdbr::msigdbr(species = "Homo sapiens", category = collection)
    else msigdbr::msigdbr(species = "Homo sapiens", category = collection, subcategory = sub)
  })
}
to_list <- function(tab) {
  split(tab$ensembl_gene, tab$gs_name)
}
pw <- c(
  lapply(to_list(msig_one("H")), function(g) unique(g)),
  lapply(to_list(tryCatch(msig_one("C2", "REACTOME"),
                          error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"))),
         unique),
  lapply(to_list(tryCatch(msig_one("C5", "GO:BP"),
                          error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"))),
         unique)
)
names(pw) <- ifelse(duplicated(names(pw)), paste0(names(pw), "_dup"), names(pw))
universe <- unique(h3$gene)

ora_dir <- function(genes, label) {
  genes <- unique(genes)
  if (length(genes) < 15) {
    message("Skip ORA ", label, " n=", length(genes))
    return(data.table())
  }
  message("ORA ", label, " n=", length(genes))
  fg <- as.data.table(fgsea::fora(pw, genes = genes, universe = universe,
                                  minSize = 10, maxSize = 500))
  fg[, set := label]
  fg[, overlap := lengths(overlapGenes)]
  fg[, overlapGenes := NULL]
  fg[order(padj)]
}

up5 <- cons[conserved5 == TRUE & direction == "up", gene]
dn5 <- cons[conserved5 == TRUE & direction == "down", gene]
ora <- rbindlist(list(
  ora_dir(up5, "conserved_up"),
  ora_dir(dn5, "conserved_down")
), fill = TRUE)
fwrite(ora, file.path(OUT, "drive_species_hmp_conserved_ora.csv"))
cat("ORA padj<0.1:", ora[padj < 0.1, .N], "\n")

pretty_pw <- function(nm) {
  nm <- sub("^(HALLMARK|REACTOME|GOBP)_", "", nm)
  nm <- gsub("_", " ", nm)
  stringr::str_to_sentence(tolower(nm))
}
is_cycle <- function(nm) {
  str_detect(tolower(nm), "mitot|chromatid|chromosome|spindle|e2f|g2m|cell cycle|prometaphase|metaphase|meiot")
}
ora_plot <- ora[padj < 0.1][order(padj)]
ora_plot[, collection := fcase(
  startsWith(pathway, "HALLMARK"), "Hallmark",
  startsWith(pathway, "REACTOME"), "Reactome",
  default = "GO BP"
)]
keep <- rbindlist(list(
  ora_plot[set == "conserved_up"][order(padj)][1:8],
  ora_plot[set == "conserved_down" & collection == "Hallmark"][order(padj)][1:6],
  ora_plot[set == "conserved_down" & collection == "Reactome" & !is_cycle(pathway)][order(padj)][1:4],
  ora_plot[set == "conserved_down" & collection == "GO BP" & !is_cycle(pathway)][order(padj)][1:3],
  ora_plot[set == "conserved_down" & is_cycle(pathway)][order(padj)][1]
), fill = TRUE)
keep <- unique(keep, by = c("set", "pathway"))
keep <- keep[!is.na(pathway)]
keep[, name := pretty_pw(pathway)]
keep[, stars := get_sig_star(padj)]
keep[, signed := fifelse(set == "conserved_up", -log10(pmax(padj, 1e-20)),
                         log10(pmax(padj, 1e-20)))]
keep[, name := fct_reorder(str_wrap(name, 42), signed)]

p_ora <- ggplot(keep, aes(signed, name, fill = set, shape = collection)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(size = 3.3, color = "grey25") +
  geom_text(aes(label = stars, x = signed - 0.85 * sign(signed)), size = 3.0) +
  scale_fill_manual(
    values = c(conserved_up = unname(gradient_high),
               conserved_down = unname(gradient_low)),
    labels = c(conserved_up = "Up in ≥5 species", conserved_down = "Down in ≥5 species"),
    name = NULL
  ) +
  scale_shape_manual(values = c(Hallmark = 21, Reactome = 22, `GO BP` = 23), name = NULL) +
  labs(
    title = "Hypoxia / glycolysis up, sterol biosynthesis down",
    subtitle = sprintf("ORA of Drive cross-perturbation genes conserved in >=5 species (%d up, %d down).",
                       length(up5), length(dn5)),
    x = "signed −log10(FDR)  (down ← 0 → up)", y = NULL,
    caption = sig_caption
  ) +
  theme_pub(12) +
  theme(legend.position = "bottom")
save_png(p_ora, "species_hmp_conserved_ora.png", 9.0, 7.2)

message("Wrote conservation + ORA to ", OUT)
