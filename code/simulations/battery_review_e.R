## battery_review_e.R -- checks of the three B = 50 review cells run by run_M_battery.R into battery/_review/e
## (review item e): output columns against battery_names() and E4, replicate order and seeds, parallel output equal to
## serial battery_one() (worker independence, E0.4), and the block summary recomputed independently (rejection,
## size-adjusted power from the matched null, null size, declined, events < G subsets).
## Writes battery/_review/review_e.log.
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
OUT <- edge_battery("_review"); ROOT <- file.path(OUT, "e")
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_battery_cells.R"))
sink(file.path(OUT, "review_e.log"), split = TRUE)
cat("battery_review_e.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
RNGkind("L'Ecuyer-CMRG")
Cells <- battery_cells()
B <- 50L
ALL <- list(c("2", "probit_auc_n3500"), c("3", "stk_asym_n1000"), c("4", "sparse49_cloglog_n200"),
            c("1b", "null_auc_n3500"), c("1a", "null_link_n1000"), c("1b", "sparse49_n200"))
readP <- function(b, cn) fread(file.path(ROOT, b, paste0(cn, "_pvalues.csv.gz")))

cat("1. output files\n")
for (k in ALL) {
  ce <- as.list(Cells[Cells$block == k[1] & Cells$cell == k[2], ]); P <- readP(k[1], k[2])
  want <- c("rep", "seed", battery_names(ce))
  miss <- setdiff(want, names(P)); extra <- setdiff(names(P), want)
  flags <- grep("^flag\\.", names(P), value = TRUE)
  pcols <- setdiff(names(P), c("rep", "seed", "n", "events", "fit_ok", "glm_conv", flags))
  rng_ok <- all(sapply(pcols, function(cn) all(is.na(P[[cn]]) | (P[[cn]] >= 0 & P[[cn]] <= 1))))
  cat(sprintf("  %-3s %-22s rows %d, columns %d (missing %s, extra %s), order %s, seeds = seed_base + rep %s, n constant %s,\n      %d flag columns [%s], %d p-value columns all in [0, 1] or NA: %s, rep errors %d, fit_ok all %s\n",
      k[1], k[2], nrow(P), ncol(P), if (length(miss)) paste(miss, collapse = ",") else "none", if (length(extra)) paste(extra, collapse = ",") else "none",
      identical(as.integer(P$rep), 1:B), all(P$seed == ce$seed_base + P$rep), all(P$n == ce$n), length(flags), paste(flags, collapse = " "),
      length(pcols), rng_ok, sum(P$flag.rep_error %in% 1), all(P$fit_ok == 1)))
  na <- sapply(pcols, function(cn) sum(!is.finite(P[[cn]])))
  if (any(na > 0)) cat("      no p-value counts: ", paste(sprintf("%s %d", names(na)[na > 0], na[na > 0]), collapse = ", "), "\n")
  ## parallel output = serial battery_one
  md <- 0; nam <- 0
  for (r in c(1L, 23L, 50L)) {
    v <- battery_one(r, ce); pr <- unlist(P[r, names(v), with = FALSE])
    both <- is.finite(v) & is.finite(pr); md <- max(md, abs(v[both] - pr[both])); nam <- nam + sum(is.finite(v) != is.finite(pr))
  }
  cat(sprintf("      serial battery_one() for reps 1, 23, 50 vs the parallel file: max |diff| %.1e, NA mismatches %d\n", md, nam))
}

cat("\n2. summaries recomputed\n")
META <- c("rep", "seed", "n", "events", "fit_ok", "glm_conv")
for (k in ALL[1:3]) {
  ce <- as.list(Cells[Cells$block == k[1] & Cells$cell == k[2], ])
  S <- fread(file.path(ROOT, k[1], "_summary.csv"))[cell == k[2]]
  P <- readP(k[1], k[2]); N <- readP(ce$null_block, ce$null_cell)
  cat(sprintf("  %-3s %-22s matched null %s %s; summary rows %d; roles %s; null_cell column %s; columns: %s\n", k[1], k[2], ce$null_block, ce$null_cell,
              nrow(S), paste(unique(S$role), collapse = "/"), paste(unique(S$null_cell), collapse = "/"), paste(names(S), collapse = ",")))
  tests <- setdiff(names(P), META); tests <- tests[!grepl("^(flag\\.|lam[0-9]|chk\\.)", tests)]
  subsets <- list(all = rep(TRUE, nrow(P)))
  for (a in grep("^flag\\.evlt\\.", names(P), value = TRUE)) {
    fl <- P[[a]] %in% 1
    if (any(fl) && any(!fl)) { g <- sub("flag.evlt.", "", a, fixed = TRUE); subsets[[paste0("events<", g)]] <- fl; subsets[[paste0("events>=", g)]] <- !fl }
  }
  cat(sprintf("      subsets expected: %s; in summary: %s\n", paste(names(subsets), collapse = ", "), paste(unique(S$subset), collapse = ", ")))
  dmax <- c(rejection = 0, size_adj_power = 0, null_size = 0, declined = 0, mcse = 0); miss <- 0
  for (t in tests) for (sn in names(subsets)) for (a in c(0.01, 0.05, 0.10)) {
    ss <- subsets[[sn]]; p <- P[[t]][ss]; pn <- N[[t]]
    pn1 <- sort(ifelse(is.finite(pn), pn, 1)); crit <- pn1[ceiling(a * length(pn1) - 1e-9)]
    ref <- c(rejection = mean(is.finite(p) & p <= a), size_adj_power = mean(is.finite(p) & p <= crit),
             null_size = mean(is.finite(pn) & pn <= a), declined = mean(!is.finite(p)))
    ref["mcse"] <- sqrt(ref[["rejection"]] * (1 - ref[["rejection"]]) / length(p))
    row <- S[test == t & subset == sn & abs(alpha - a) < 1e-12]
    if (nrow(row) != 1) { miss <- miss + 1; next }
    for (nm in names(dmax)) dmax[nm] <- max(dmax[nm], abs(row[[nm]] - ref[[nm]]))
  }
  cat(sprintf("      max |summary - own| over %d tests x %d subsets x 3 alphas: %s; rows not found %d\n", length(tests), length(subsets),
              paste(sprintf("%s %.1e", names(dmax), dmax), collapse = ", "), miss))
  sh <- S[subset == "all" & abs(alpha - 0.05) < 1e-12 & test %in% c("EDGE.poly3.u.Grule", "EDGE.sym.sc.Grule", "Stk.joint", "Stk.LR", "GiViTI", "HL.G10", "PR")]
  print(sh[, .(test, B, rejection, size_adj_power, null_size, declined, status)], row.names = FALSE)
}

cat("\n3. progress logs\n")
for (b in c("1a", "1b", "2", "3", "4")) { cat("  --", b, "\n"); cat(paste0("  ", readLines(file.path(ROOT, b, "_progress.log"))), sep = "\n") }
sink()
