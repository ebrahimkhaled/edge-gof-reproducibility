## robust_influence_check2.R -- the bound of ROBUSTNESS_bounded_influence.md 2(b), checked two ways.
## (i) membership fixed: the grouping is taken from the clean fit and kept, so only one member's prediction changes;
## (ii) membership as the test really builds it: the corrupted record moves to the top group and displaces one record
##      per boundary, so two groups change by at most one record each.
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
suppressMessages({ library(ebrahim.gof) })
OUT <- edge_path("declarations/robust_influence_bound.csv")
tight <- glm.control(epsilon = 1e-12, maxit = 100)
options(width = 190)

set.seed(20260916)
n <- 1000
x0 <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
y <- rbinom(n, 1, plogis(0.6 * x0 + 0.5 * d))
i0 <- which(y == 0 & abs(x0) < 0.5)[1]
GRID <- c(x0[i0], 2, 3, 4, 6, 8, 12, 16, 24, 40)
G <- max(10, round(n / 25))
grp_of <- function(ph) pmin(ceiling(rank(ph, ties.method = "first") / (length(ph) / G)), G)

fitp <- function(xv) {
  x <- x0; x[i0] <- xv
  fit <- suppressWarnings(glm(y ~ x + d, family = binomial(), control = tight))
  pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6)
}
res_of <- function(pr, g, gi) { ii <- which(g == gi); V <- sum(pr[ii] * (1 - pr[ii])); c(V = V, r = (sum(y[ii]) - sum(pr[ii])) / sqrt(V)) }

pr0 <- fitp(GRID[1]); g0 <- grp_of(pr0); gi0 <- g0[i0]                      # the clean grouping and the record's clean group
base_fixed <- res_of(pr0, g0, gi0)

R <- do.call(rbind, lapply(GRID, function(xv) {
  pr <- fitp(xv)
  ## (i) membership fixed at the clean grouping
  f <- res_of(pr, g0, gi0)
  bound_fixed <- (1 + abs(base_fixed[["r"]]) / (8 * base_fixed[["V"]])) / sqrt(base_fixed[["V"]])
  ## (ii) the grouping the test really uses
  g <- grp_of(pr); moved <- sum(g != g0)
  a <- res_of(pr, g, gi0)                                                   # the record's old group, re-formed
  top <- res_of(pr, g, G); top0 <- res_of(pr0, g0, G)                       # the top group, where the corrupted record lands
  data.frame(x_corrupt = xv, moved_records = moved,
             V_fixed = f[["V"]], r_fixed = f[["r"]], delta_fixed = abs(f[["r"]] - base_fixed[["r"]]), bound_fixed = bound_fixed,
             V_regrouped = a[["V"]], r_regrouped = a[["r"]], delta_regrouped = abs(a[["r"]] - base_fixed[["r"]]),
             V_top = top[["V"]], r_top = top[["r"]], delta_top = abs(top[["r"]] - top0[["r"]]),
             bound_top = 2 / sqrt(top[["V"]]))
}))
write.csv(R, OUT, row.names = FALSE)
cat("(i) membership fixed at the clean grouping: only the corrupted record's prediction changes inside its group\n")
print(format(R[, c("x_corrupt", "V_fixed", "r_fixed", "delta_fixed", "bound_fixed")], digits = 4), row.names = FALSE)
cat("\n(ii) the grouping the test really uses: the corrupted record leaves its group and lands at the top\n")
print(format(R[, c("x_corrupt", "moved_records", "V_regrouped", "r_regrouped", "delta_regrouped", "V_top", "r_top", "delta_top", "bound_top")], digits = 4), row.names = FALSE)
cat("\nwritten:", OUT, "\n")
