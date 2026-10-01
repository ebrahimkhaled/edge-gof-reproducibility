## analyse_blockTW.R -- the reading fixed in PREDECLARATION_blockTW_mechanism.md and its two addenda.
##   false alarms: rejection at 5% in the corrupted cells (holds at <= 0.10); size: clean logit cells;
##   power: raw, and size-adjusted at the clean logit cell of the same n and design.
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
suppressMessages(library(data.table))
rd <- function(dir) rbindlist(lapply(list.files(edge_battery(dir), "pvalues.csv.gz$", full.names = TRUE),
  function(f) { d <- fread(f); d[, cell := sub("_pvalues.csv.gz$", "", basename(f))]; d }), fill = TRUE)
TW <- rd("TW")
HH <- if (dir.exists(edge_battery("TW_HH"))) rd("TW_HH") else NULL
if (!is.null(HH) && nrow(HH)) TW <- merge(TW, HH[, !"EDGE.G10.check"], by = c("cell", "rep"), all.x = TRUE)
TW[, c("part", "truth", "corruption", "kk", "nn") := tstrsplit(cell, "_")]
TW[, n := as.integer(sub("n", "", nn))]; TW[, k := as.integer(sub("k", "", kk))]
ARMS <- intersect(c("EDGE.default", "EDGE.imhof", "EDGE.G10", "HL.G10", "HH.HLnp", "HH.X2np", "TWIN.score", "EDGE.Gn",
                    "SPZ", "Stk.joint", "Cubic.LR", "GiViTI", "PRG.Stk", "PRG.EDGE.G10", "PRG2.Stk", "PRG2.Cubic",
                    "PRG2.EDGE.default", "PRG2.EDGE.G10"), names(TW))
L <- melt(TW, id.vars = c("cell", "part", "truth", "corruption", "k", "n", "rep"), measure.vars = ARMS,
          variable.name = "test", value.name = "p")
L[, null_cell := sprintf("%s_logit_clean_k00_n%d", part, n)]
crit <- L[cell == null_cell, .(crit = as.numeric(quantile(p, 0.05, type = 1, na.rm = TRUE))), by = .(null_cell, test)]
L <- merge(L, crit, by = c("null_cell", "test"), all.x = TRUE)
S <- L[, .(raw = mean(!is.na(p) & p <= 0.05), adj = mean(!is.na(p) & p <= crit), missing = mean(is.na(p)), reps = .N),
       by = .(cell, part, truth, corruption, k, n, test)]
fwrite(S, edge_battery("TW", "_summary.csv"))
W <- dcast(S, part + n + truth + corruption + k ~ test, value.var = "raw")
setcolorder(W, c("part", "n", "truth", "corruption", "k", ARMS))
cat("\n--- rejection at 5% (raw) ---\n"); print(W, digits = 3)
cat("\n--- size-adjusted power, clean misfit cells ---\n")
print(dcast(S[corruption == "clean" & truth != "logit"], part + n + truth ~ test, value.var = "adj"), digits = 3)
cat("\n--- records dropped by the diagnostics rules (mean) ---\n")
print(TW[, .(PRG = mean(PRG.dropped, na.rm = TRUE), PRG_corrupt = mean(PRG.dropped_corrupt, na.rm = TRUE),
             PRG2 = mean(PRG2.dropped, na.rm = TRUE), PRG2_corrupt = mean(PRG2.dropped_corrupt, na.rm = TRUE)), by = cell])
