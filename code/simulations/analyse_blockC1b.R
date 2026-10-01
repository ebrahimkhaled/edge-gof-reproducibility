## analyse_blockC1b.R -- block Z Part A read against its declaration (sha256 f5e1a7f1...).
## Writes battery/C1b/_summary.csv and prints the comparison the paper will quote: the sign error
## against the exaggeration, test by test, at the same k.
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
suppressPackageStartupMessages(library(data.table))
setwd(SIMDIR)
BAT <- edge_battery()
source(file.path(SIMDIR, "_battery_cells.R"))

TESTS <- c(EDGE = "EDGE.poly3.u.Grule", `EDGE-sym` = "EDGE.sym.u.Grule", `HL (rule G)` = "HL.Grule",
           `HL (G=10)` = "HL.G10", `Stukel joint` = "Stk.joint", GiViTI = "GiViTI",
           `cubic LR` = "Cubic.LR")

## one read per file: block 9's deposits are large and every test asks the same file again
RD <- new.env(parent = emptyenv())
rate <- function(f, col) {
  if (!file.exists(f)) return(NA_real_)
  key <- basename(f)
  if (is.null(RD[[key]])) RD[[key]] <- fread(f)
  M <- RD[[key]]
  if (!col %in% names(M)) return(NA_real_)
  mean(M[[col]] <= 0.05, na.rm = TRUE)
}

rows <- list()
for (n in c(1000L, 5000L)) for (k in c(0L, as.integer(c(0.001, 0.002, 0.005, 0.01) * n))) {
  if (k == 0L) { c1b <- sprintf("c1b_clean_n%d", n); c1 <- sprintf("logit_clean_n%d", n) }
  else {
    ## both blocks name their cells from a VECTOR of rates, so 0.01 is written "010", not "01";
    ## rebuilding the tag one rate at a time would give the wrong name for the one-per-cent cells
    tag <- sprintf("%03d", as.integer(round(1000 * k / n)))
    c1b <- sprintf("c1b_r%s_n%d", tag, n); c1 <- sprintf("logit_C1_r%s_n%d", tag, n)
  }
  for (nm in names(TESTS)) {
    rows[[length(rows) + 1]] <- data.table(
      n = n, k = k, test = nm,
      sign_error  = rate(edge_battery("C1b", sprintf("%s_pvalues.csv.gz", c1b)), TESTS[[nm]]),
      exaggerated = rate(edge_battery("9", sprintf("%s_pvalues.csv.gz", c1)), TESTS[[nm]]))
  }
}
S <- rbindlist(rows)
S[, difference := sign_error - exaggerated]
fwrite(S, edge_battery("C1b", "_summary.csv"))

cat("\n== the sign error (C1b) against the exaggeration (C1), rejection at 0.05 ==\n")
for (nn in c(1000L, 5000L)) {
  cat(sprintf("\n-- n = %d --\n", nn))
  W <- dcast(S[n == nn], k ~ test, value.var = "sign_error")
  cat("sign error:\n"); print(W, digits = 3)
  W2 <- dcast(S[n == nn], k ~ test, value.var = "exaggerated")
  cat("exaggeration (block 9):\n"); print(W2, digits = 3)
}

cat("\n== the declared reading: does the directed test exceed 0.10 at ten corrupted records? ==\n")
x <- S[test == "EDGE" & ((n == 1000L & k == 10L) | (n == 5000L & k == 10L))]
print(x[, .(n, k, sign_error, exaggerated)])
cat(sprintf("\nA2 verdict: %s\n",
            if (any(x$sign_error > 0.10)) "GOES AGAINST the paper's scope claim, as declared"
            else "the scope claim survives the sign error"))
cat("\nlevel on clean data (both sample sizes, every test):\n")
print(dcast(S[k == 0L], test ~ n, value.var = "sign_error"), digits = 3)
