## battery_review_tests2.R -- follow-up to battery_review_tests.R: is the 3e-6 gap between the harness score form and the
## anova Rao reference (null_base_n1000, replicate 1, rule G) the default glm convergence, or the formula?
## Applies the harness functions to a tightly converged fit and compares again. Writes battery/_review/review_tests2.log.
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

SIMDIR <- edge_path("code/simulations")
OUT <- edge_battery("_review")
source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
sink(file.path(OUT, "review_tests2.log"), split = TRUE)
Cells <- battery_cells()
ce <- as.list(Cells[Cells$block == "1b" & Cells$cell == "null_base_n1000", ])
tight <- stats::glm.control(epsilon = 1e-14, maxit = 200)
RNGkind("L'Ecuyer-CMRG"); set.seed(ce$seed_base + 1); dat <- bt_data(ce)
fd <- suppressWarnings(stats::glm(dat$f, data = dat$d, family = stats::binomial()))
f0 <- suppressWarnings(stats::glm(dat$f, data = dat$d, family = stats::binomial(), control = tight))
ft <- suppressWarnings(stats::glm(dat$f, data = dat$d, family = stats::binomial(), control = tight, start = stats::coef(f0)))
cat(sprintf("default fit: %d iterations, max |coef - tight coef| %.2e\n", fd$iter, max(abs(stats::coef(fd) - stats::coef(ft)))))
mkfq <- function(fit) { eta <- as.numeric(fit$linear.predictors); ph <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6)
  dmu <- fit$family$mu.eta(eta); X <- stats::model.matrix(fit); w <- dmu^2 / (ph * (1 - ph))
  list(fit = fit, y = as.numeric(fit$y), eta = eta, ph = ph, p_raw = as.numeric(fitted(fit)), dmu = dmu, X = X, A = crossprod(X, w * X), n = length(fit$y)) }
G <- 40
for (lab in c("default", "tight")) {
  fq <- mkfq(if (lab == "default") fd else ft); gs <- bt_groups(fq, G)
  for (b in BT_BASES) {
    Z <- bt_basis(gs$pbar, b); sc <- bt_edge_score(gs, fq$A, Z)
    d2 <- dat$d; for (j in seq_len(ncol(Z))) d2[[paste0("zz", j)]] <- Z[gs$grp, j]
    m1 <- suppressWarnings(stats::glm(stats::update(dat$f, stats::as.formula(paste(". ~ . +", paste0("zz", seq_len(ncol(Z)), collapse = "+")))),
                                      data = d2, family = stats::binomial(), control = tight))
    a <- stats::anova(ft, m1, test = "Rao")
    cat(sprintf("  %-7s fit, %-5s: harness S %.10f p %.10f | anova Rao S %.10f p %.10f | |dp| %.1e\n", lab, b, sc$S, sc$p, a$Rao[2],
                a[["Pr(>Chi)"]][2], abs(sc$p - a[["Pr(>Chi)"]][2])))
  }
}
sink()
