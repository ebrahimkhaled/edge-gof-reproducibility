## robust_influence_rivals.R -- the influence curve of robust_influence_check.R, repeated for Liu et al.'s (2024) projection
## test and for BAGofT, with EDGE, Stukel and GiViTI on the same data. One fixed data set (n = 500, so the slow rivals finish);
## one record's covariate is moved along a grid while its outcome stays 0. Exploratory. Output: theory/robust_influence_rivals.csv
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
suppressMessages({ library(ebrahim.gof); library(givitiR); library(BAGofT) })
source(edge_path("code/simulations/_proj_test.R"))
OUT <- edge_path("declarations/robust_influence_rivals.csv")
tight <- glm.control(epsilon = 1e-12, maxit = 100)
options(width = 200)

set.seed(20260917)
n <- 500
x0 <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
y <- rbinom(n, 1, plogis(0.6 * x0 + 0.5 * d))
i0 <- which(y == 0 & abs(x0) < 0.5)[1]
GRID <- c(x0[i0], 4, 8, 12, 16, 24, 40)

rows <- list()
for (k in seq_along(GRID)) {
  x <- x0; x[i0] <- GRID[k]
  dat <- data.frame(y = y, x = x, d = d)
  fit <- suppressWarnings(glm(y ~ x + d, family = binomial(), data = dat, control = tight))
  pr <- as.numeric(fitted(fit)); e <- as.numeric(predict(fit, type = "link"))
  za <- 0.5 * e^2 * (pr >= 0.5); zb <- -0.5 * e^2 * (pr < 0.5)
  f1 <- suppressWarnings(glm(y ~ x + d + za + zb, family = binomial(), data = dat, start = c(coef(fit), 0, 0), control = tight))
  pv <- function(f) tryCatch(f$p_value, error = function(e) NA_real_)

  set.seed(1000 + k)
  t0 <- proc.time()[["elapsed"]]
  p_proj <- tryCatch(proj_pvalue(y, model.matrix(fit), B = 250)$p_value, error = function(e) NA_real_)
  s_proj <- proc.time()[["elapsed"]] - t0

  set.seed(2000 + k)
  t0 <- proc.time()[["elapsed"]]
  p_bag <- tryCatch({ invisible(capture.output(b <- suppressMessages(suppressWarnings(
                        BAGofT(testModel = testGlmBi(formula = y ~ x + d, link = "logit"), data = dat)))))
                      b$p.value }, error = function(e) NA_real_)
  s_bag <- proc.time()[["elapsed"]] - t0

  rows[[k]] <- data.frame(x_corrupt = GRID[k], eta_corrupt = e[i0], b_x = unname(coef(fit)[2]),
    Stk.joint = tryCatch(anova(fit, f1, test = "Rao")[2, "Pr(>Chi)"], error = function(e) NA_real_),
    GiViTI = tryCatch(suppressWarnings(givitiCalibrationTest(y, pr, devel = "internal")$p.value), error = function(e) NA_real_),
    EDGE.poly3.u = pv(def.gof(fit, G = "auto")), EDGE.sym.u = pv(def.gof(fit, G = "auto", basis = "sym")),
    HL.G10 = tryCatch(run.all.gof(fit, tests = "HL", G = 10)$p_value[1], error = function(e) NA_real_),
    projection = p_proj, BAGofT = p_bag, sec_projection = round(s_proj, 1), sec_BAGofT = round(s_bag, 1))
  cat(sprintf("x = %5.1f  proj %.3f (%.0fs)  BAGofT %.3f (%.0fs)\n", GRID[k], p_proj, s_proj, p_bag, s_bag))
}
R <- do.call(rbind, rows)
write.csv(R, OUT, row.names = FALSE)
cat("\none record's covariate is moved; its outcome stays 0; n = 500; p-values:\n")
print(format(R, digits = 3), row.names = FALSE)
cat("\nwritten:", OUT, "\n")
