## pkg280_verify_num_warm.R -- joint Stukel and the DEF score form against anova(..., test = "Rao") on WARM-RESTARTED
## fits, using the INSTALLED ebrahim.gof 2.8.0. Verification of the package, not a paper result.
##
## anova.glmlist computes Rao from m1$residuals and m1$weights. glm.fit returns weights computed before its last
## IRLS step, and it stops on a deviance change, so that last step is about sqrt(epsilon * deviance): 4e-3 at the
## default epsilon, 1e-5 at epsilon = 1e-13. Refitting from the converged coefficients makes the last step ~1e-10,
## so anova's Rao becomes the exact u'I^-1 u and can be compared with the package at a tight tolerance.
##
## Helpers and fixtures are taken from pkg280_verify_num_exact.R (only the named assignments are evaluated).
## Run: Rscript pkg280_verify_num_warm.R > ../paper_EDGE/theory/pkg280_verify_num_warm.log 2>&1
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

suppressPackageStartupMessages(library(ebrahim.gof))
PKGV <- as.character(packageVersion("ebrahim.gof"))
stopifnot(PKGV == "2.8.0")
OUT <- edge_path("declarations")
options(width = 200, warn = 1)
cat("ebrahim.gof", PKGV, "(installed) |", R.version.string, "|", format(Sys.time()), "\n")

SRC <- edge_path("code/simulations/pkg280_verify_num_exact.R")
data("gof_demo", package = "ebrahim.gof")
TAKE <- c("mk", "find_fit", "augment", "score_qr", "stukel_Z", "stk", "Gof", "step_cov", "count_warn", "quiet",
          "demo_wrong", "A", "one_group_half", "BD")
for (e in parse(SRC, keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% TAKE)
    eval(e, globalenv())
stopifnot(all(TAKE %in% ls(globalenv())))

## tight fit, then one more fit started at its coefficients
warm <- function(fit) {
  env <- environment(formula(fit))
  cl <- getCall(fit); cl$control <- glm.control(epsilon = 1e-13, maxit = 100)
  ft <- suppressWarnings(eval(cl, env))
  cl$start <- unname(coef(ft))
  suppressWarnings(eval(cl, env))
}

ROWS <- list()
rec <- function(part, design, quantity, pkg, ref)
  ROWS[[length(ROWS) + 1]] <<- data.frame(part = part, design = design, quantity = quantity, pkg = pkg, ref = ref,
                                          abs_diff = abs(pkg - ref), rel_diff = abs(pkg - ref) / max(abs(ref), 1e-12),
                                          stringsAsFactors = FALSE)

cat("\n==== fixtures of (a): joint Stukel vs anova Rao, warm restart ====\n")
for (nm in names(A)) {
  fw <- warm(A[[nm]]); j <- stk(fw, "joint")
  a <- anova(fw, augment(fw, stukel_Z(fw)), test = "Rao")
  cat(sprintf("  %-10s joint %.12f df %g | anova %.12f df %g | |diff| %.2e\n", nm, j$Statistic, j$df, a$Rao[2], a$Df[2],
              abs(j$Statistic - a$Rao[2])))
  rec("a", nm, "joint stat", j$Statistic, a$Rao[2]); rec("a", nm, "joint df", j$df, a$Df[2])
}

cat("\n==== designs of (b): score form vs anova Rao on the step covariates, warm restart ====\n")
for (d in BD) {
  fw <- warm(d$fit)
  for (b in c("sym", "poly3", "poly2", "stukel")) {
    r <- quiet(def.gof(fw, G = d$G, basis = b, weights = "score"))
    a <- anova(fw, augment(fw, step_cov(fw, d$G, b)), test = "Rao")
    cat(sprintf("  %-44s %-6s S %.12f df %g | anova %.12f df %g | |diff| %.2e\n", d$name, b, r$Test_Statistic, r$df,
                a$Rao[2], a$Df[2], abs(r$Test_Statistic - a$Rao[2])))
    rec("b", d$name, paste(b, "score stat"), r$Test_Statistic, a$Rao[2]); rec("b", d$name, paste(b, "score df"), r$df, a$Df[2])
  }
}

cat("\n==== random data sets from the size cells (500 per cell), warm restart ====\n")
gen <- function(n, c0, s) {
  dat <- data.frame(xa = runif(n, -3, 3), db = rbinom(n, 1, 0.5))
  dat$out <- rbinom(n, 1, plogis(c0 + s * (0.6 * dat$xa + 0.5 * dat$db)))
  suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat))
}
CELLS <- data.frame(label = c("base n1000 G10", "base n1000 Gauto", "c0=-2 n500 Gauto", "s=2 n1000 Gauto"),
                    s = c(1, 1, 1, 2), c0 = c(0, 0, -2, 0), n = c(1000, 1000, 500, 1000), G = c("10", "auto", "auto", "auto"),
                    stringsAsFactors = FALSE)
for (k in seq_len(nrow(CELLS))) {
  ce <- CELLS[k, ]; G <- if (ce$G == "auto") "auto" else as.numeric(ce$G)
  set.seed(777000 + k)
  D <- matrix(NA_real_, 500, 10, dimnames = list(NULL, c("joint", "joint_df", "sym", "sym_df", "poly3", "poly3_df",
                                                         "poly2", "poly2_df", "stukel", "stukel_df")))
  for (i in 1:500) {
    fw <- warm(gen(ce$n, ce$c0, ce$s))
    j <- stk(fw, "joint"); a <- anova(fw, augment(fw, stukel_Z(fw)), test = "Rao")
    D[i, 1:2] <- c(abs(j$Statistic - a$Rao[2]) / max(a$Rao[2], 1e-8), j$df != a$Df[2])
    for (b in c("sym", "poly3", "poly2", "stukel")) {
      r <- quiet(def.gof(fw, G = G, basis = b, weights = "score"))
      a <- anova(fw, augment(fw, step_cov(fw, G, b)), test = "Rao")
      D[i, paste0(b, c("", "_df"))] <- c(abs(r$Test_Statistic - a$Rao[2]) / max(a$Rao[2], 1e-8), r$df != a$Df[2])
    }
  }
  cat(sprintf("  %-18s max rel |pkg - anova|: joint %.2e (df mismatches %d) | sym %.2e (%d) | poly3 %.2e (%d) | poly2 %.2e (%d) | stukel %.2e (%d)\n",
              ce$label, max(D[, "joint"]), sum(D[, "joint_df"]), max(D[, "sym"]), sum(D[, "sym_df"]),
              max(D[, "poly3"]), sum(D[, "poly3_df"]), max(D[, "poly2"]), sum(D[, "poly2_df"]),
              max(D[, "stukel"]), sum(D[, "stukel_df"])))
  for (v in c("joint", "sym", "poly3", "poly2", "stukel"))
    ROWS[[length(ROWS) + 1]] <- data.frame(part = "sim", design = ce$label, quantity = paste(v, "max rel diff over 500"),
                                           pkg = NA, ref = NA, abs_diff = NA, rel_diff = max(D[, v]),
                                           stringsAsFactors = FALSE)
}

W <- do.call(rbind, ROWS); W$pkg_version <- PKGV
write.csv(W, file.path(OUT, "pkg280_verify_num_warm_checks.csv"), row.names = FALSE)
cat(sprintf("\nfixture rows: max |diff| %.2e, max rel %.2e over %d comparisons\n",
            max(W$abs_diff, na.rm = TRUE), max(W$rel_diff[W$part != "sim"], na.rm = TRUE), sum(W$part != "sim")))
cat("done", format(Sys.time()), "\n")
