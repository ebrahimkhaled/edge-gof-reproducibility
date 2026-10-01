## edge_external.R -- the directed test in external mode, as Section 7 of EDGE paper 3 uses it.
##
## When a published model is checked on new patients with its coefficients frozen, nothing is estimated
## from those patients: Omega = I exactly, and no score equation absorbs the overall level, so a constant
## column joins the basis and the statistic has d + 1 degrees of freedom. ebrahim.gof 2.8.0's
## edge.gof(y, predicted_probs = p) computes the d-column statistic with Omega = I and calls that
## reference conservative, because it is written for the case where the fit is unknown. This function
## is the external mode the paper reports; it is written to be folded into def.gof() as external = TRUE.
##
##   source("edge_external.R")
##   edge_external(y, p, G = 10, basis = "poly3")     # y outcomes, p frozen predicted probabilities
##
## Returns the statistic, its degrees of freedom, the p-value and the grouped residuals r_g.
## ---- archive paths (inserted by make_archive.py; the convention is in README.md, "How to run") -------------
EDGE_ARCHIVE_ROOT <- local({
  r <- Sys.getenv("EDGE_ARCHIVE_ROOT")
  if (!nzchar(r)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    r <- if (length(f)) file.path(dirname(normalizePath(f[1], winslash = "/")), "..", "..") else getwd()
  }
  r <- normalizePath(r, winslash = "/", mustWork = FALSE)
  Sys.setenv(EDGE_ARCHIVE_ROOT = r)        # so worker processes started from here resolve the same root
  r
})
edge_path <- function(...) file.path(EDGE_ARCHIVE_ROOT, ...)
edge_battery <- function(...) {           # the author's battery/ folder: analysis/ -> results/analysis, rest -> results/blocks
  if (...length() == 0L) return(edge_path("results", "blocks"))
  p <- file.path(...)
  ifelse(p == "analysis" | startsWith(p, "analysis/"), edge_path("results", p), edge_path("results", "blocks", p))
}
edge_out <- function(...) { d <- edge_path("output", ...); dir.create(d, showWarnings = FALSE, recursive = TRUE); d }
## ---------------------------------------------------------------------------------------------------------------
edge_external <- function(y, p, G = 10, basis = c("poly3", "stukel", "sym")) {
  basis <- match.arg(basis)
  stopifnot(length(y) == length(p), all(y %in% c(0, 1)), all(p > 0 & p < 1))
  if (identical(G, "auto")) G <- max(10L, as.integer(ceiling(length(y) / 25)))
  ## equal-frequency groups on the frozen risk, by rank, as def.gof() forms them
  g <- as.integer(pmin(ceiling(rank(p, ties.method = "first") / (length(p) / G)), G))
  M <- rowsum(cbind(y, p, p * (1 - p), 1), g, reorder = TRUE)
  r <- (M[, 1] - M[, 2]) / sqrt(M[, 3])
  pb <- M[, 2] / M[, 4]
  eta <- stats::qlogis(pb)
  Z <- switch(basis,
    poly3  = cbind(1, pb, pb^2, pb^3),                                        # the default basis + constant
    stukel = cbind(1, eta, eta^2 * (eta >= 0), -eta^2 * (eta < 0)),
    sym    = cbind(1, eta * abs(eta)))
  Q <- qr(Z)
  Z <- Z[, Q$pivot[seq_len(Q$rank)], drop = FALSE]
  S <- sum(qr.fitted(qr(Z), r)^2)                                             # Omega = I: chi-squared on ncol(Z)
  list(statistic = S, df = ncol(Z), p_value = stats::pchisq(S, ncol(Z), lower.tail = FALSE),
       groups = nrow(M), residuals = as.numeric(r), mean_risk = as.numeric(pb))
}

## self-check against the paper's Section 7 (run_I_bigdata_p3.R), when run as a script
if (sys.nframe() == 0L) {
  setwd(edge_path("code/simulations"))
  source("_cohort_p3.R")
  CO <- cohort_p3(getwd())
  fit <- stats::glm(CO$f, data = CO$Ddev, family = stats::binomial())
  p <- as.numeric(stats::predict(fit, newdata = CO$Dval, type = "response"))
  V <- read.csv(edge_path("results/cohort/runI_p3_val.csv"))
  for (G in V$G) {
    e <- edge_external(CO$Dval$y, p, G)
    ref <- V[V$G == G, ]
    cat(sprintf("G=%5d  S=%8.3f (paper %8.3f)  df=%d  p=%.3g (paper %.3g)\n", G, e$statistic, ref$poly_stat,
                e$df, e$p_value, ref$poly_p))
    stopifnot(abs(e$statistic - ref$poly_stat) < 1e-8)
  }
  cat("edge_external reproduces Section 7 at every partition\n")
}
