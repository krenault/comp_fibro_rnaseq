#!/usr/bin/env Rscript
# =============================================================================
# READING COPY. Not used by the analysis.
#
# The file the scripts actually source is share/hmp_combine.R.
# This copy is the same math, with the comments I would want if I were
# reading it for the first time. Run the example at the bottom if you want
# to watch one gene go through the steps:
#
#   Rscript share/hmp_combine_annotated.R
#
# Why the p-value is two-sided
# ---------------------------
# A DESeq2 or mixed-model p-value is already two-sided. It says the gene
# moved. It does not say whether it moved up or down.
#
# Split each p in half. The half on the observed side is evidence for that
# direction. The other half is evidence for the opposite direction.
#
# Combine the up-halves. Combine the down-halves. Keep the smaller one.
# Then double it. If the winning side was 0.01, report 0.02.
#
# The doubling is the whole reason it is still called two-sided. Both
# directions were tried, and the winner was kept. Without the double, every
# gene would look twice as significant as the test allows. Cap the result at 1.
#
# The sign rule is a separate filter, and it is stricter. A gene that is
# strongly up in two stresses and strongly down in two can still have a
# small combined p, because each side listens to its own small halves.
# The 80% rule removes that gene. With 3 or 4 stresses, 80% means every
# stress that was measured points the same way.
#
# Which calculator
# ----------------
# harmonicmeanp::p.hmp, version that produced the tables already on disk
# (package 3.0.1 when those tables were checked). Equal weight on each
# test that is actually present. multilevel = FALSE.
#
# multilevel = FALSE is the important setting. In that mode the package
# IGNORES L. L is still an argument, and the analysis scripts still pass
# L = 4 (public classes) or L = 3 (FH experiment heat/cold/hypoxia), but the
# combined p-value depends only on the p-values you hand it. A gene
# measured in 3 classes is combined as 3 tests, not as 3 tests inside a
# family of 4. Do not flip multilevel to TRUE: that would change the
# p-values and the gene list.
#
# There used to be a Landau-tail approximation for machines without the
# package. It is gone. If harmonicmeanp is not installed, this stops.
# The 894-gene list matches p.hmp, not the approximation.
#
# The three rules that turn a combined p into "in the set"
# --------------------------------------------------------
# Applied AFTER the p-value, and after Benjamini-Hochberg across genes:
#
#   1. BH-adjusted HMP < 0.05
#   2. at least 3 tests survived (finite p, finite effect, p > 0)
#   3. at least 80% of those tests share the winning sign
#
# With 3 or 4 tests, rule 3 is unanimity:
#   3 of 4 = 75%, which fails 80%
#   2 of 3 = 67%, which fails 80%
#   4 of 5 = 80%, which would pass, but this analysis does not have 5
# So a public generic gene is up in every class it was measured in, or
# down in every class it was measured in, and it was measured in at least 3.
#
# What this function does NOT do
# ------------------------------
# It does not fit the mixed models. For the public list, each "test" is
# already one class (heat, hypoxia, oxidative, DNA damage): one p and one
# effect from metafor, after many studies were combined. This function
# only combines those four numbers.
#
# For one species in the FH experiments, each "test" is that species' own DESeq2 contrast
# (heat 41 vs 37, cold 32 vs 37, and one hypoxia). Same function, different
# rows. That table is drive_species_hmp_three_insult.csv.
# =============================================================================

# One gene, one call.
#
# p         two-sided p-values, one per test
# sign_ec   signed effects, same order as p. Only the sign is used here.
#           The size of the effect does not enter the p-value.
# L         accepted so older calls still work. Not used by p.hmp when
#           multilevel = FALSE. See the header.
hmp_twosided <- function(p, sign_ec, L = length(p)) {
  # Drop a test if its p or its effect is missing, infinite, or if p is 0.
  # p = 0 would make the harmonic mean explode (it divides by p). DESeq2
  # sometimes reports p = 0 for very strong effects; those tests are left out
  # rather than treated as infinitely significant.
  ok <- is.finite(p) & is.finite(sign_ec) & p > 0
  p <- p[ok]
  s <- sign(sign_ec[ok])

  # One test is not a combination. The caller will not call this gene a hit,
  # because n will be 0 or 1, and the minimum is 3.
  if (length(p) < 2L) {
    return(list(p = NA_real_, dir = 0, n = length(p), frac_agree = NA_real_))
  }

  # An effect of exactly 0 has no sign. Counting it as up is a corner case.
  # Real log2 fold-changes are almost never bitwise zero.
  s[s == 0] <- 1

  # Turn each two-sided p into a one-sided p for "up" and a one-sided p for "down".
  #
  # If the gene went UP in this test (s > 0):
  #   evidence for up   = p/2          (the small tail)
  #   evidence for down = 1 - p/2      (almost 1: this test argues against down)
  # If the gene went DOWN, those two are swapped.
  #
  # The pmin/pmax lines only keep the number inside (0, 1) so log and 1/p
  # stay defined. They are not a scientific threshold.
  one_pos <- ifelse(s > 0, pmin(pmax(p / 2, 1e-300), 1 - 1e-16),
                    pmin(pmax(1 - p / 2, 1e-300), 1 - 1e-16))
  one_neg <- ifelse(s < 0, pmin(pmax(p / 2, 1e-300), 1 - 1e-16),
                    pmin(pmax(1 - p / 2, 1e-300), 1 - 1e-16))

  # Harmonic mean p of one side.
  #
  # Intuition: the harmonic mean is dragged down by the smallest p-values.
  # Three moderate p-values do less than one very small p-value and two
  # ordinary ones. That is the opposite of a plain average, where every
  # test pulls equally. Wilson 2019 turns that harmonic mean into a
  # p-value that still means "probability under the null".
  #
  # Equal weights: each test that is present counts the same. A class with
  # 20 studies does not count more here than a class with 3 studies. That
  # study-size difference already happened inside the mixed model, before
  # this function saw the number.
  #
  # multilevel = FALSE: stand-alone combination. L is ignored. We then
  # Benjamini-Hochberg the resulting p-values across genes. We are not using
  # Wilson's multilevel family-wise error procedure.
  hmp_one <- function(ps) {
    if (!requireNamespace("harmonicmeanp", quietly = TRUE)) {
      stop(
        "Package harmonicmeanp is required. ",
        "This analysis does not use a Landau-tail approximation."
      )
    }
    as.numeric(harmonicmeanp::p.hmp(
      ps, w = rep(1 / length(ps), length(ps)), L = L, multilevel = FALSE
    ))
  }

  pp <- hmp_one(one_pos)   # combined evidence that the shared direction is up
  pn <- hmp_one(one_neg)   # combined evidence that the shared direction is down
  if (!is.finite(pp)) pp <- 1
  if (!is.finite(pn)) pn <- 1

  # Keep the side with the smaller combined p. Double it to make a two-sided p.
  # frac_agree is the fraction of tests whose sign matches that winning side.
  # It is not a p-value. The 80% rule is applied later, by the caller.
  if (pp <= pn) {
    list(p = min(1, 2 * pp), dir = 1, n = length(p), frac_agree = mean(s > 0))
  } else {
    list(p = min(1, 2 * pn), dir = -1, n = length(p), frac_agree = mean(s < 0))
  }
}

# Whole table.
#
# x is one row per test, not one row per gene. Typical columns:
#   gene, class, p, effect
# for the public analysis, class is heat / hypoxia / oxidative / dna_damage
# and effect is the class mixed-model estimate.
#
# For a FH experiment species the same shape works with
#   gene, exposure, pvalue, log2FoldChange
# and the column-name guesses below.
#
# padj_cutoff, min_n, sign_agree are the three rules. Defaults are the
# ones used for both the public list and the species list.
hmp_combine <- function(x, id = NULL, p = NULL, effect = NULL, test = NULL,
                        L = NULL, padj_cutoff = 0.05, min_n = 3, sign_agree = 0.80) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  # Column names are guessed so a DESeq2 table and a class-model table can
  # both be passed without renaming. An explicit argument overrides the guess.
  if (is.null(id)) id <- .pick_col(x, c("gene", "ensembl", "id", "symbol"))
  if (is.null(p)) p <- .pick_col(x, c("p", "pvalue", "p_value", "pval"))
  if (is.null(effect)) {
    effect <- .pick_col(x, c("effect", "ec", "log2fc", "log2FoldChange", "logFC", "estimate", "lfc"))
  }
  if (is.null(test)) {
    hit <- intersect(c("test", "class", "signature", "exposure"), names(x))
    test <- if (length(hit)) hit[[1]] else NULL
  }

  # L is stored and passed through. It does not change p.hmp while
  # multilevel = FALSE. Default: number of distinct test labels, or the
  # most tests any gene has.
  if (is.null(L)) {
    L <- if (!is.null(test)) length(unique(x[[test]][!is.na(x[[test]])])) else max(table(x[[id]]))
  }

  ids <- unique(x[[id]])
  rows <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    sub <- x[x[[id]] == ids[[i]], , drop = FALSE]
    h <- hmp_twosided(sub[[p]], sub[[effect]], L = L)
    ev <- sub[[effect]]
    rows[[i]] <- data.frame(
      hmp_p = h$p,             # two-sided combined p, before BH
      hmp_dir = h$dir,         # +1 up, -1 down, 0 if not combined
      n_test = h$n,            # how many tests survived the finite/p>0 filter
      frac_agree = h$frac_agree,
      effect_mean = mean(ev[is.finite(ev)]),  # descriptive only; not used in the call
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out <- cbind(setNames(data.frame(ids, stringsAsFactors = FALSE), id), out)
  rownames(out) <- NULL

  # Benjamini-Hochberg across genes in this table.
  # Public analysis: one BH across all genes.
  # Species analysis: 05_species_hmp.R does this BH separately per species.
  # If you paste every species into one CSV and call hmp_combine once, the
  # BH will be across species-gene rows and will NOT match that file.
  out$hmp_padj <- p.adjust(out$hmp_p, method = "BH")

  out$in_set <- is.finite(out$hmp_padj) & out$hmp_padj < padj_cutoff &
    out$n_test >= min_n & !is.na(out$frac_agree) & out$frac_agree >= sign_agree
  out$direction <- ifelse(out$in_set & out$hmp_dir > 0, "up",
                          ifelse(out$in_set & out$hmp_dir < 0, "down", "none"))
  out
}

.pick_col <- function(x, choices) {
  hit <- intersect(choices, names(x))
  if (!length(hit)) {
    stop("Need a column named one of: ", paste(choices, collapse = ", "))
  }
  hit[[1]]
}

# -----------------------------------------------------------------------------
# Worked example. Four genes, four classes. This is the public rule (L = 4
# is passed and ignored). Read the comments next to each gene, then the
# printed table.
#
# GENE_UP
#   all four effects positive, all four p-values small
#   winning side = up, frac_agree = 4/4 = 1
#   passes BH, passes n >= 3, passes 80%
#
# GENE_DOWN
#   same pattern, signs negative
#
# GENE_SPLIT
#   very small p-values, but two up and two down
#   the winning side can still have a tiny combined p, because each side
#   sees two small p/2 values and two near-1 values, and the harmonic mean
#   listens to the small ones
#   frac_agree = 2/4 = 0.5, so the 80% rule throws it out
#   this is the whole point of "direction-aware" plus the sign rule:
#   a small p is not enough
#
# GENE_THREE_OF_FOUR
#   three up, one down. frac_agree = 3/4 = 0.75, which is under 80%
#   not a hit, even though most classes agree
#
# GENE_ONLY_TWO
#   both up and both tiny, but n_test = 2, under the minimum of 3
# -----------------------------------------------------------------------------
if (length(grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)) &&
    grepl("hmp_combine_annotated\\.R$",
          sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]))) {
  example <- data.frame(
    gene = c(
      rep("GENE_UP", 4),
      rep("GENE_DOWN", 4),
      rep("GENE_SPLIT", 4),
      rep("GENE_THREE_OF_FOUR", 4),
      rep("GENE_ONLY_TWO", 2)
    ),
    class = c(
      rep(c("heat", "hypoxia", "oxidative", "dna_damage"), 4),
      c("heat", "hypoxia")
    ),
    p = c(
      1e-6, 1e-5, 1e-4, 1e-3,
      1e-4, 1e-3, 0.01, 0.02,
      1e-6, 1e-6, 1e-4, 1e-3,
      1e-6, 1e-5, 1e-4, 1e-3,
      1e-8, 1e-8
    ),
    effect = c(
      0.8, 0.5, 0.4, 0.2,
      -0.6, -0.4, -0.3, -0.2,
      0.9, -0.9, 0.4, -0.3,
      0.8, 0.5, 0.4, -0.2,
      1.0, 1.0
    ),
    stringsAsFactors = FALSE
  )
  res <- hmp_combine(example, L = 4)
  print(res, digits = 4)
  message(
    "Hits should be GENE_UP and GENE_DOWN only. ",
    "GENE_SPLIT fails the sign rule, GENE_THREE_OF_FOUR is 75% agreement, ",
    "GENE_ONLY_TWO has only two tests."
  )
}
