#!/usr/bin/env Rscript
# Direction-aware two-sided harmonic mean p (Wilson, PNAS 2019; Yoon).
# Same rule as the generic-stress signature:
#   BH HMP < 0.05, at least 3 tests, at least 80% sign agreement.
#
# One row per test. Columns (names are detected):
#   gene          id
#   p             two-sided p
#   effect        signed effect (log2FC, class EC, ...)
#   class         optional test label. L defaults to the number of labels.
#
#   Rscript share/hmp_combine.R share/example_tests.csv share/example_hmp.csv
#   Rscript share/hmp_combine.R class_p.csv hmp.csv --L 4 --min-n 3 --sign-agree 0.8
#
# Requires the harmonicmeanp package. That is the only p-value this file
# computes. A long commented copy of the same math is
# share/hmp_combine_annotated.R.

hmp_twosided <- function(p, sign_ec, L = length(p)) {
  ok <- is.finite(p) & is.finite(sign_ec) & p > 0
  p <- p[ok]
  s <- sign(sign_ec[ok])
  if (length(p) < 2L) {
    return(list(p = NA_real_, dir = 0, n = length(p), frac_agree = NA_real_))
  }
  s[s == 0] <- 1
  one_pos <- ifelse(s > 0, pmin(pmax(p / 2, 1e-300), 1 - 1e-16),
                    pmin(pmax(1 - p / 2, 1e-300), 1 - 1e-16))
  one_neg <- ifelse(s < 0, pmin(pmax(p / 2, 1e-300), 1 - 1e-16),
                    pmin(pmax(1 - p / 2, 1e-300), 1 - 1e-16))
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
  pp <- hmp_one(one_pos)
  pn <- hmp_one(one_neg)
  if (!is.finite(pp)) pp <- 1
  if (!is.finite(pn)) pn <- 1
  if (pp <= pn) {
    list(p = min(1, 2 * pp), dir = 1, n = length(p), frac_agree = mean(s > 0))
  } else {
    list(p = min(1, 2 * pn), dir = -1, n = length(p), frac_agree = mean(s < 0))
  }
}

.pick_col <- function(x, choices) {
  hit <- intersect(choices, names(x))
  if (!length(hit)) {
    stop("Need a column named one of: ", paste(choices, collapse = ", "))
  }
  hit[[1]]
}

hmp_combine <- function(x, id = NULL, p = NULL, effect = NULL, test = NULL,
                        L = NULL, padj_cutoff = 0.05, min_n = 3, sign_agree = 0.80) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (is.null(id)) id <- .pick_col(x, c("gene", "ensembl", "id", "symbol"))
  if (is.null(p)) p <- .pick_col(x, c("p", "pvalue", "p_value", "pval"))
  if (is.null(effect)) {
    effect <- .pick_col(x, c("effect", "ec", "log2fc", "log2FoldChange", "logFC", "estimate", "lfc"))
  }
  if (is.null(test)) {
    hit <- intersect(c("test", "class", "signature", "exposure"), names(x))
    test <- if (length(hit)) hit[[1]] else NULL
  }
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
      hmp_p = h$p, hmp_dir = h$dir, n_test = h$n, frac_agree = h$frac_agree,
      effect_mean = mean(ev[is.finite(ev)]),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out <- cbind(setNames(data.frame(ids, stringsAsFactors = FALSE), id), out)
  rownames(out) <- NULL
  out$hmp_padj <- p.adjust(out$hmp_p, method = "BH")
  out$in_set <- is.finite(out$hmp_padj) & out$hmp_padj < padj_cutoff &
    out$n_test >= min_n & !is.na(out$frac_agree) & out$frac_agree >= sign_agree
  out$direction <- ifelse(out$in_set & out$hmp_dir > 0, "up",
                          ifelse(out$in_set & out$hmp_dir < 0, "down", "none"))
  out
}

.hmp_parse_args <- function(args) {
  opt <- list(
    input = NULL, output = NULL, L = NULL, min_n = 3, sign_agree = 0.80,
    padj = 0.05, id = NULL, p = NULL, effect = NULL, test = NULL
  )
  positional <- character()
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    take <- function() {
      i <<- i + 1L
      if (i > length(args)) stop("Missing value after ", a)
      args[[i]]
    }
    if (a == "--L") opt$L <- as.integer(take())
    else if (a == "--min-n") opt$min_n <- as.integer(take())
    else if (a == "--sign-agree") opt$sign_agree <- as.numeric(take())
    else if (a == "--padj") opt$padj <- as.numeric(take())
    else if (a == "--id") opt$id <- take()
    else if (a == "--p") opt$p <- take()
    else if (a == "--effect") opt$effect <- take()
    else if (a == "--test") opt$test <- take()
    else if (startsWith(a, "-")) stop("Unknown option: ", a)
    else positional <- c(positional, a)
    i <- i + 1L
  }
  if (length(positional) >= 1L) opt$input <- positional[[1]]
  if (length(positional) >= 2L) opt$output <- positional[[2]]
  opt
}

if (length(grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)) &&
    grepl("hmp_combine\\.R$", sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]))) {
  opt <- .hmp_parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(opt$input)) {
    cat(
      "Usage: Rscript hmp_combine.R tests.csv [hmp.csv]",
      "       [--L 4] [--min-n 3] [--sign-agree 0.8] [--padj 0.05]\n",
      "Input: one row per test, columns gene, p, effect, and optional class.\n",
      "Example: Rscript share/hmp_combine.R share/example_tests.csv share/example_hmp.csv --L 4\n",
      sep = ""
    )
    quit(status = 0)
  }
  tab <- utils::read.csv(opt$input, stringsAsFactors = FALSE)
  res <- hmp_combine(
    tab, id = opt$id, p = opt$p, effect = opt$effect, test = opt$test,
    L = opt$L, padj_cutoff = opt$padj, min_n = opt$min_n, sign_agree = opt$sign_agree
  )
  if (is.null(opt$output)) {
    utils::write.csv(res, row.names = FALSE)
  } else {
    utils::write.csv(res, opt$output, row.names = FALSE)
    message("Wrote ", nrow(res), " genes, ", sum(res$in_set), " in the set, to ", opt$output)
  }
}
