## robust_influence_check.R -- what one record can do to each calibration test.
## One fixed data set; one record's covariate is moved along a grid while its outcome stays as drawn. For every value we
## record each test's p-value, the fitted coefficients, and, for EDGE, the contaminated group's residual, variance and the
## bound of ROBUSTNESS_bounded_influence.md Section 2(b). Output: theory/robust_influence_check.csv
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
suppressMessages({ library(ebrahim.gof); library(givitiR) })
OUT <- edge_path("declarations/robust_influence_check.csv")
tight <- glm.control(epsilon = 1e-12, maxit = 100)
options(width = 200)

set.seed(20260916)
n <- 1000
x0 <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
y <- rbinom(n, 1, plogis(0.6 * x0 + 0.5 * d))
i0 <- which(y == 0 & abs(x0) < 0.5)[1]                     # a middling record whose outcome is 0
GRID <- c(x0[i0], 2, 3, 4, 6, 8, 12, 16, 24, 40)           # the corrupted covariate value

grp_of <- function(ph, G) pmin(ceiling(rank(ph, ties.method = "first") / (length(ph) / G)), G)

one <- function(xv) {
  x <- x0; x[i0] <- xv
  fit <- suppressWarnings(glm(y ~ x + d, family = binomial(), control = tight))
  pr <- as.numeric(fitted(fit)); e <- as.numeric(predict(fit, type = "link"))
  za <- 0.5 * e^2 * (pr >= 0.5); zb <- -0.5 * e^2 * (pr < 0.5); zs <- e * abs(e)
  rao <- function(extra) {
    f1 <- suppressWarnings(glm(y ~ x + d + extra, family = binomial(), start = c(coef(fit), rep(0, NCOL(extra))), control = tight))
    tryCatch(anova(fit, f1, test = "Rao")[2, "Pr(>Chi)"], error = function(e) NA_real_)
  }
  lr <- function(extra) {
    f1 <- suppressWarnings(glm(y ~ x + d + extra, family = binomial(), start = c(coef(fit), rep(0, NCOL(extra))), control = tight))
    tryCatch(anova(fit, f1, test = "LRT")[2, "Pr(>Chi)"], error = function(e) NA_real_)
  }
  g0 <- suppressWarnings(glm(y ~ e, family = binomial()))
  g1 <- suppressWarnings(glm(y ~ e + I(e^2) + I(e^3), family = binomial()))
  pv <- function(f) tryCatch(f$p_value, error = function(e) NA_real_)
  ## the contaminated record's group, as def.gof builds it (G = max(10, round(n/25)))
  G <- max(10, round(n / 25)); g <- grp_of(pmin(pmax(pr, 1e-6), 1 - 1e-6), G)
  ii <- which(g == g[i0])
  Vg <- sum(pr[ii] * (1 - pr[ii])); rg <- (sum(y[ii]) - sum(pr[ii])) / sqrt(Vg)
  data.frame(x_corrupt = xv, eta_corrupt = e[i0], p_corrupt = pr[i0],
             b_x = unname(coef(fit)[2]), b_d = unname(coef(fit)[3]),
             Stk.joint = rao(cbind(za, zb)), Stk.sym1 = rao(zs), Stk.LR = lr(cbind(za, zb)),
             Cubic.LR = tryCatch(anova(g0, g1, test = "LRT")[2, "Pr(>Chi)"], error = function(e) NA_real_),
             GiViTI = tryCatch(suppressWarnings(givitiCalibrationTest(y, pr, devel = "internal")$p.value), error = function(e) NA_real_),
             EDGE.poly3.u = pv(def.gof(fit, G = "auto")), EDGE.sym.u = pv(def.gof(fit, G = "auto", basis = "sym")),
             EDGE.poly3.sc = pv(def.gof(fit, G = "auto", weights = "score")),
             EDGE.sym.sc = pv(def.gof(fit, G = "auto", basis = "sym", weights = "score")),
             HL.G10 = tryCatch(run.all.gof(fit, tests = "HL", G = 10)$p_value[1], error = function(e) NA_real_),
             group_size = length(ii), V_g = Vg, r_g = rg,
             bound = (1 + abs(rg) / (8 * Vg)) / sqrt(Vg))
}

R <- do.call(rbind, lapply(GRID, function(v) suppressWarnings(one(v))))
R$delta_r_g <- abs(R$r_g - R$r_g[1])
write.csv(R, OUT, row.names = FALSE)

cat("one record's covariate is moved; its outcome stays 0. n = 1000, G = 40, the model is correct for the other 999 records.\n\n")
cat("p-values:\n")
print(format(R[, c("x_corrupt", "eta_corrupt", "Stk.joint", "Stk.sym1", "Stk.LR", "Cubic.LR", "GiViTI",
                   "EDGE.poly3.u", "EDGE.sym.u", "EDGE.poly3.sc", "EDGE.sym.sc", "HL.G10")], digits = 3), row.names = FALSE)
cat("\nthe fit, and the contaminated group (its residual against the bound of Section 2b):\n")
print(format(R[, c("x_corrupt", "b_x", "b_d", "group_size", "V_g", "r_g", "delta_r_g", "bound")], digits = 4), row.names = FALSE)
cat("\nwritten:", OUT, "\n")
