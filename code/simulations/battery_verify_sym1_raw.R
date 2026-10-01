## battery_verify_sym1_raw.R -- after E10.1: Stk.sym1 (raw fitted risks) against anova's Rao score test for adding
## eta|eta| to the working model, on sparse samples with and without a fitted risk outside [1e-6, 1 - 1e-6].
## Output: battery/_review/verify_sym1_raw.csv and .log
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

setwd(edge_path("code/simulations"))
suppressMessages({ library(statmod) })
source("_battery_tests.R")

tight <- glm.control(epsilon = 1e-14, maxit = 200)
rows <- list(); n_clamp <- 0; n_plain <- 0
for (seed in 1:3000) {
  if (n_clamp >= 25 && n_plain >= 25) break
  set.seed(seed)
  n <- sample(c(100, 150, 200, 300), 1)
  x <- as.numeric(scale(rchisq(n, 4)))
  y <- rbinom(n, 1, plogis(-4.9 + x))
  if (sum(y) < 1 || sum(y) == n) next
  f0 <- suppressWarnings(glm(y ~ x, family = binomial()))
  fit <- suppressWarnings(glm(y ~ x, family = binomial(), start = coef(f0), control = tight))
  pr <- as.numeric(fitted(fit))
  clamp <- any(pr < 1e-6 | pr > 1 - 1e-6)
  if (clamp && n_clamp >= 25) next
  if (!clamp && n_plain >= 25) next
  eta <- as.numeric(predict(fit, type = "link"))
  fq <- list(fit = fit, y = as.numeric(fit$y), eta = eta, ph = pmin(pmax(pr, 1e-6), 1 - 1e-6), p_raw = pr,
             X = model.matrix(fit), n = n, dmu = pr * (1 - pr))
  st <- tryCatch(suppressWarnings(bt_stukel(fq)), error = function(e) NULL)
  if (is.null(st)) next
  zs <- eta * abs(eta)
  f1 <- suppressWarnings(glm(y ~ x + zs, family = binomial(), start = c(coef(fit), 0), control = tight))
  rao <- tryCatch(anova(fit, f1, test = "Rao")[2, "Pr(>Chi)"], error = function(e) NA_real_)
  rows[[length(rows) + 1]] <- data.frame(seed = seed, n = n, events = sum(y), clamp = clamp,
                                         sym1 = st[["Stk.sym1"]], rao = rao, guard = st[["flag.info_guard"]])
  if (clamp) n_clamp <- n_clamp + 1 else n_plain <- n_plain + 1
}
res <- do.call(rbind, rows)
res$abs_diff <- abs(res$sym1 - res$rao)
dir.create("battery/_review", showWarnings = FALSE, recursive = TRUE)
write.csv(res, "battery/_review/verify_sym1_raw.csv", row.names = FALSE)
sink("battery/_review/verify_sym1_raw.log")
cat("samples:", nrow(res), " clamp-active:", sum(res$clamp), "\n")
for (cl in c(TRUE, FALSE)) {
  r <- res[res$clamp == cl & is.finite(res$sym1) & is.finite(res$rao), ]
  cat(sprintf("clamp = %s: compared %d, max |sym1 - rao| %.3g, median %.3g; sym1 NA %d, guard %d\n", cl, nrow(r),
              max(r$abs_diff), median(r$abs_diff), sum(res$clamp == cl & !is.finite(res$sym1)),
              sum(res$clamp == cl & res$guard == 1)))
}
sink()
cat(readLines("battery/_review/verify_sym1_raw.log"), sep = "\n")
