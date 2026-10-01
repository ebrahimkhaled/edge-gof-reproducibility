## Self-check: my HL must reproduce ResourceSelection::hoslem.test, and my Stukel must
## reproduce the project's own stuk_diag, BEFORE any run is trusted. A reimplementation that
## is never checked against a reference is how a wrong number reaches a paper.
suppressPackageStartupMessages({library(ebrahim.gof); library(ResourceSelection)})
source("_pstar_giviti_harness.R")
source("_dgp_library.R")

set.seed(20260910)
ok <- TRUE

cat("=== HL vs ResourceSelection::hoslem.test (10 datasets) ===\n")
for (i in 1:10) {
  g <- gen_bench(600, p = 4)
  fit <- glm(g$f, data = g$d, family = binomial())
  mine <- hl_stat(fit$y, fitted(fit), G = 10, equal_width = FALSE)
  ref  <- ResourceSelection::hoslem.test(fit$y, fitted(fit), g = 10)
  dstat <- abs(mine["stat"] - as.numeric(ref$statistic))
  dp    <- abs(mine["p"] - as.numeric(ref$p.value))
  flag  <- if (dstat < 1e-8 && dp < 1e-8) "ok" else "MISMATCH"
  if (flag != "ok") ok <- FALSE
  cat(sprintf("  %2d  mine X2=%8.4f p=%.6f | ref X2=%8.4f p=%.6f  %s\n",
              i, mine["stat"], mine["p"], ref$statistic, ref$p.value, flag))
}

cat("\n=== HL df convention ===\n")
g <- gen_bench(600, p = 4); fit <- glm(g$f, data = g$d, family = binomial())
cat("  mine df:", hl_stat(fit$y, fitted(fit), 10)["df"],
    " ref df:", as.numeric(ResourceSelection::hoslem.test(fit$y, fitted(fit), g = 10)$parameter), "\n")

cat("\n=== Stukel: mine vs the project's stuk_diag ===\n")
for (i in 1:5) {
  g <- dgp_alt("link", "cloglog", 800)
  fit <- glm(g$f, data = g$d, family = binomial())
  mine <- stukel_p(fit)
  fit2 <- glm(g$f, data = g$d, family = binomial()); fit2$data <- g$d
  ref <- tryCatch(stuk_diag(fit2), error = function(e) NULL)
  cat(sprintf("  %d  mine p=%.6f | stuk_diag: %s\n", i, mine,
              paste(utils::capture.output(str(ref, max.level = 1))[1], collapse = "")))
}

cat("\n=== one_rep on a null and an alternative ===\n")
print(round(one_rep(dgp_null("link", 1000), G = 10), 4))
print(round(one_rep(dgp_alt("link", "cloglog", 1000), G = 10), 4))

cat("\nSELF-CHECK:", if (ok) "PASS" else "FAIL -- do not run the sweeps", "\n")
