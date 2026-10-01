## pkg280_verify_corr_extra.R -- two spot checks of the correction.
## 1. NEWS says of the guard: "The information of such a column is rounding noise, and inverting it could give a p-value
##    of zero." Probe: two events at mirrored x (fitted slope 0), then one event x moved by delta, so the slope is small
##    but not zero; and the same with three events. For each: the ratio, the pre-guard result (2e0f6d1) and the installed
##    result, for the Stukel joint row and the stukel and sym score forms at G = 10.
## 2. The (y, predicted_probs, X) path of run.all.gof on a sample with no event: DEF rows NA with the no-fit Note.
## Run: Rscript pkg280_verify_corr_extra.R <clone dir> > ../paper_EDGE/theory/pkg280_verify_corr_extra.log 2>&1

args <- commandArgs(trailingOnly = TRUE)
PKG <- args[1]
options(width = 220, warn = 1)
suppressPackageStartupMessages(library(ebrahim.gof))
ns <- asNamespace("ebrahim.gof")
git_show <- function(sha, file) system2("git", c("-C", shQuote(PKG), "show", paste0(sha, ":", file)), stdout = TRUE)
PREG <- new.env(parent = ns)
for (ex in parse(text = git_show("2e0f6d1", "R/def_gof.R"), keep.source = FALSE))
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]])) eval(ex, PREG)
PRES <- new.env(parent = ns)
for (ex in parse(text = git_show("2e0f6d1", "R/run_all_gof.R"), keep.source = FALSE))
  if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) && as.character(ex[[2]]) == "gof_stukel") eval(ex, PRES)
cat("installed ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "| pre-guard = 2e0f6d1 |", format(Sys.time()), "\n")

q <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) NULL)
stk_ratio <- function(fit) {
  ph <- pmin(pmax(fitted(fit), 1e-6), 1 - 1e-6); pr <- as.numeric(fitted(fit)); W <- pr * (1 - pr); eta <- fit$linear.predictors
  Z <- cbind(0.5 * eta^2 * (ph >= 0.5), -0.5 * eta^2 * (ph < 0.5)); Z <- Z[, colSums(Z != 0) > 0, drop = FALSE]
  X <- model.matrix(fit)
  I <- crossprod(Z, W * Z) - crossprod(Z, W * X) %*% solve(crossprod(X, W * X), crossprod(X, W * Z))
  min(diag(I) / colSums(W * Z^2))
}
fmt <- function(r) if (is.null(r)) "error" else sprintf("S %-12.6g df %-3s p %-10.4g", as.numeric(r[[1]]), format(r[[2]]), as.numeric(r[[3]]))

cat("\n==== 1. near-flat fits: pre-guard (2e0f6d1) against installed ====\n")
x0 <- qnorm(ppoints(200))
for (ev in list(c(1, 200), c(1, 100, 200))) for (delta in c(0, 1e-12, 1e-9, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1)) {
  x <- x0; y <- rep(0, 200); y[ev] <- 1; x[ev[1]] <- x[ev[1]] + delta
  fit <- suppressWarnings(glm(y ~ x, family = binomial()))
  cat(sprintf("\nevents at rows %s, delta %g: slope %.3g, sd(eta) %.3g, Stukel ratio %.3g\n", paste(ev, collapse = ","), delta,
              coef(fit)[2], sd(fit$linear.predictors), q(stk_ratio(fit))))
  ctx <- ns$.gof_context(fit)
  s0 <- q(PRES$gof_stukel(ctx)); s1 <- q(ns$gof_stukel(ctx))
  cat(sprintf("  Stukel joint  pre %s | new %s | %s\n", fmt(s0[c("Statistic", "df", "p_value")]), fmt(s1[c("Statistic", "df", "p_value")]), s1$Note))
  for (b in c("stukel", "sym")) {
    d0 <- q(PREG$def.gof(fit, basis = b, weights = "score")); d1 <- q(def.gof(fit, basis = b, weights = "score"))
    cat(sprintf("  %-6s score  pre %s | new %s\n", b, fmt(if (is.null(d0)) NULL else d0[c("Test_Statistic", "df", "p_value")]),
                fmt(if (is.null(d1)) NULL else d1[c("Test_Statistic", "df", "p_value")])))
  }
}

cat("\n==== 2. run.all.gof(y, predicted_probs, X) on a sample with no event ====\n")
set.seed(9); x <- runif(200, -3, 3); X <- cbind(1, x); y <- rep(0, 200); ph <- plogis(-6 + 0.2 * x)
w <- character(0)
tb <- withCallingHandlers(run.all.gof(y, predicted_probs = ph, X = X, include_slow = FALSE),
                          warning = function(cn) { w <<- c(w, conditionMessage(cn)); invokeRestart("muffleWarning") })
tb <- as.data.frame(lapply(tb, function(v) v), stringsAsFactors = FALSE)   # drop the gof_battery print method
print(tb[tb$Test %in% c("DEF.poly2", "DEF.poly3", "DEF.stukel", "DEF.sym", "Stukel") | grepl("^Ensemble", tb$Test), c("Test", "p_value", "Note")], row.names = FALSE)
cat("warnings escaping:", if (length(w)) paste(unique(substr(w, 1, 100)), collapse = " / ") else "none", "\n")
cat("\ndone", format(Sys.time()), "\n")
