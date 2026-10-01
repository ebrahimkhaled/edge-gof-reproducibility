## battery_review_zero.R -- follow-up to battery_review_degenerate.R: samples with no event.
##  (a) constructed zero-event data sets (x = standardised chi-square(4) as gen_sparse_link, y = 0) for n in
##      {100, 150, 200, 300, 500}, 100 each, through battery_rep(): per test, how often a p-value is returned, how often
##      it is exactly 0, how often it is <= 0.05;
##  (b) the installed ebrahim.gof 2.8.0 on five of them (def.gof score and unit forms, gof_stukel joint), to see whether
##      the package behaves the same way;
##  (c) for every cell of the table that uses the sparse generator: event rate from a 400,000 draw, P(no event) =
##      (1 - rate)^n, and the size inflation this implies for the affected tests, against the Monte Carlo SE at the cell's B.
## Writes battery/_review/review_zero.log and review_zero.csv.
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
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
sink(file.path(OUT, "review_zero.log"), split = TRUE)
options(width = 220)
cat("battery_review_zero.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
RNGkind("L'Ecuyer-CMRG"); set.seed(20260914L)
Cells <- battery_cells()
ce0 <- as.list(Cells[Cells$block == "1b" & Cells$cell == "sparse49_n100", ])

## (a)
NS <- c(100L, 150L, 200L, 300L, 500L)
rows <- list(); keep <- list()
for (n in NS) for (k in 1:100) {
  x <- rchisq(n, 4); x <- as.numeric(scale(x))
  dat <- list(d = data.frame(x = x, y = rep(0, n)), f = y ~ x)
  ce <- ce0; ce$n <- n
  v <- battery_rep(dat, ce)
  rows[[length(rows) + 1]] <- as.data.frame(t(c(n_set = n, k = k, v)))
  if (n == 200L && k <= 5) keep[[k]] <- dat
}
R <- rbindlist(rows, fill = TRUE)
fwrite(R, file.path(OUT, "review_zero.csv"))
tc <- setdiff(names(R), c("n_set", "k", "n", "events", "fit_ok", "glm_conv")); tc <- tc[!grepl("^flag", tc)]
cat("(a) zero-event samples, 100 per n: share returning a p-value / share with p = 0 / share with p <= 0.05\n")
S <- rbindlist(lapply(tc, function(t) {
  z <- R[, .(ret = mean(is.finite(get(t))), zero = mean(is.finite(get(t)) & get(t) == 0), rej = mean(is.finite(get(t)) & get(t) <= 0.05)), by = n_set]
  data.table(test = t, txt = paste(sprintf("n=%d %.2f/%.2f/%.2f", z$n_set, z$ret, z$zero, z$rej), collapse = "  "), anyrej = max(z$rej))
}))
print(S[order(-anyrej), .(test, txt)], row.names = FALSE)
cat(sprintf("\nglm converged in %.0f%% of the zero-event fits; fitted probability range %s\n", 100 * mean(R$glm_conv == 1),
            "see (b)"))

## (b)
cat("\n(b) installed ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "on five zero-event data sets, n = 200\n")
ns <- asNamespace("ebrahim.gof")
for (k in seq_along(keep)) {
  dat <- keep[[k]]
  fit <- suppressWarnings(stats::glm(y ~ x, data = dat$d, family = stats::binomial()))
  pk <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) paste("error:", conditionMessage(e)))
  a <- pk(ebrahim.gof::def.gof(fit, G = 10, basis = "stukel", weights = "score")$p_value)
  b <- pk(ebrahim.gof::def.gof(fit, G = 10, basis = "sym", weights = "score")$p_value)
  c1 <- pk(ebrahim.gof::def.gof(fit, G = 10, basis = "stukel", weights = "unit")$p_value)
  j <- pk(get("gof_stukel", envir = ns)(get(".gof_context", envir = ns)(fit, G = 10), list(form = "joint"))$p_value)
  h <- battery_rep(dat, c(ce0[setdiff(names(ce0), "n")], list(n = 200L)))
  cat(sprintf("  set %d: glm coef %s, fitted p %.2e | package: stk score %s, sym score %s, stk unit %s, Stukel joint %s | harness: %s %s %s %s\n",
              k, paste(signif(coef(fit), 3), collapse = " "), fitted(fit)[1], format(a, digits = 3), format(b, digits = 3), format(c1, digits = 3),
              format(j, digits = 3), format(h[["EDGE.stk.sc.G10"]], digits = 3), format(h[["EDGE.sym.sc.G10"]], digits = 3),
              format(h[["EDGE.stk.u.G10"]], digits = 3), format(h[["Stk.joint"]], digits = 3)))
}

## (c)
cat("\n(c) cells using the sparse generator: P(no event) and implied size inflation at alpha = 0.05\n")
rej_by_n <- R[, lapply(.SD, function(p) mean(is.finite(p) & p <= 0.05)), by = n_set, .SDcols = tc]
aff <- tc[vapply(tc, function(t) max(rej_by_n[[t]]) > 0, logical(1))]
sp <- Cells[Cells$generator == "sparse", ]
out <- list()
for (i in seq_len(nrow(sp))) {
  ce <- as.list(sp[i, ]); ce$n <- 400000L
  set.seed(777L + i); dd <- bt_data(ce)$d
  rate <- mean(dd$y); n <- sp$n[i]
  p0 <- (1 - rate)^n
  nn <- NS[which.min(abs(NS - n))]
  q <- unlist(rej_by_n[n_set == nn, aff, with = FALSE])
  mcse <- sqrt(0.05 * 0.95 / sp$B[i])
  out[[i]] <- data.frame(block = sp$block[i], cell = sp$cell[i], n = n, B = sp$B[i], link = sp$link[i], event_rate = rate, p_no_event = p0,
                         max_inflation = p0 * max(q), worst_test = names(q)[which.max(q)], mcse_05 = mcse,
                         inflation_over_3mcse = p0 * max(q) > 3 * mcse, stringsAsFactors = FALSE)
}
O <- do.call(rbind, out)
print(O[O$p_no_event > 1e-6, ], row.names = FALSE, digits = 3)
cat(sprintf("\ncells with P(no event) > 1e-6: %d of %d sparse-generator cells; inflation above 3 MCSE in %d\n",
            sum(O$p_no_event > 1e-6), nrow(O), sum(O$inflation_over_3mcse)))
sink()
