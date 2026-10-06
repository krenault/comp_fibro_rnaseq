# Direction-aware harmonic mean p (Yoon / Wilson 2019; Cell 2023 aging meta).
# Used by the public generic signature and by per-species Drive HMP.

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
    if (requireNamespace("harmonicmeanp", quietly = TRUE)) {
      as.numeric(harmonicmeanp::p.hmp(
        ps, w = rep(1 / length(ps), length(ps)), L = L, multilevel = FALSE
      ))
    } else {
      h <- length(ps) / sum(1 / ps)
      mu <- log(L) + 1 - digamma(1)
      sig <- pi / 2
      z <- (1 / h - mu) / sig
      min(1, max(h, stats::pnorm(z, lower.tail = FALSE) * 1.5))
    }
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
