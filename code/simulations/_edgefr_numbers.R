## _edgefr_numbers.R -- the numbers Section 6.4 quotes, read from the deposits rather than typed.
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
BAT <- edge_battery()
R <- fread(edge_battery("analysis", "edgefr.csv"))
P <- fread(edge_battery("analysis", "edgefr_paired.csv"))

cat("== size-adjusted power, per scenario, default against EDGE-FR(1%) and ten groups ==\n")
W <- dcast(R[role %in% c("power (battery)", "power")], cell ~ variant, value.var = "adjusted")
W[, loss_FR := V0 - HYB05][, loss_G10 := V0 - G10]
print(W[order(-loss_FR), .(cell, V0 = round(V0, 3), HYB05 = round(HYB05, 3), G10 = round(G10, 3),
                           loss_FR = round(loss_FR, 3), loss_G10 = round(loss_G10, 3))])

cat("\n== paired, size-adjusted: EDGE-FR(1%) against the default ==\n")
print(P[a == "HYB05" & b == "V0", .(cell, FR = round(a_power, 3), default = round(b_power, 3),
                                    diff = round(diff, 3), FR_only = a_only, default_only = b_only,
                                    holm = signif(p_holm, 2))][order(diff)])
cat("\n== paired, size-adjusted: EDGE-FR(1%) against ten groups ==\n")
print(P[a == "HYB05" & b == "G10", .(cell, FR = round(a_power, 3), G10 = round(b_power, 3),
                                     diff = round(diff, 3), FR_only = a_only, G10_only = b_only,
                                     holm = signif(p_holm, 2))][order(-diff)])

cat("\n== block 9R, the four sharper remedies (nominal 0.05) ==\n")
V9 <- c("V0", "FR1", "FR2", "DROP1", "POOL1", "G10")
f9 <- Sys.glob(edge_battery("9R", "*_variants.csv.gz"))
B9 <- rbindlist(lapply(f9, function(f) {
  M <- fread(f); cl <- sub("_variants.csv.gz$", "", basename(f))
  data.table(cell = cl, variant = V9, rate = vapply(V9, function(v) mean(M[[v]] <= 0.05, na.rm = TRUE), 0))
}))
BATT <- c("cauchit_n1000", "t4_n1000", "loglog_n1000", "stk_short_n1000", "stk_long_n1000", "stk_asym_n1000")
cat("mean over the six battery scenarios (raw):\n")
print(round(dcast(B9[cell %in% BATT], . ~ variant, value.var = "rate", fun.aggregate = mean)[, -1], 4))
cat("\nprotection at n = 5000 (raw), k = 25 and k = 50:\n")
print(dcast(B9[cell %in% c("logit_C1_r005_n5000", "logit_C1_r010_n5000")], cell ~ variant, value.var = "rate"),
      digits = 3)
cat("\nlevel, clean logistic:\n")
print(dcast(B9[cell %in% c("logit_clean_n1000", "logit_clean_n5000")], cell ~ variant, value.var = "rate"),
      digits = 3)
cat("\ncomplementary log-log at n = 1000 (the cell the text quotes):\n")
print(dcast(B9[cell == "cloglog_clean_n1000"], cell ~ variant, value.var = "rate"), digits = 3)
