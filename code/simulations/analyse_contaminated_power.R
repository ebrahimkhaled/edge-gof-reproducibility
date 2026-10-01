## analyse_contaminated_power.R -- power under contamination, read the way Section 4.3 declares.
##
## Section 6.5 reported the cost of the protection as a drop in RAW rejection rate against a real
## misfit. At the cells it quotes the test's own false-alarm rate under the SAME contamination is
## 0.151 and 0.740, so a raw rate there is not power: part of it is the corruption. The declared rule
## is to read power at the critical value of the test's own matched null, and under contamination the
## matched null is the CONTAMINATED null -- the logistic truth at the same n and the same k, where the
## model is correct and the same records are corrupted.
##
## Writes battery/analysis/contaminated_power.csv.
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

TESTS <- c(EDGE = "EDGE.poly3.u.Grule", `EDGE-sym` = "EDGE.sym.u.Grule")
RD <- new.env(parent = emptyenv())
rd <- function(cell) {
  if (is.null(RD[[cell]])) {
    f <- edge_battery("9", sprintf("%s_pvalues.csv.gz", cell))
    RD[[cell]] <- if (file.exists(f)) fread(f) else NULL
  }
  RD[[cell]]
}
nm <- function(truth, k, n) {
  if (k == 0L) sprintf("%s_clean_n%d", truth, n)
  else sprintf("%s_C1_r%03d_n%d", truth, as.integer(round(1000 * k / n)), n)
}

pv <- function(truth, k, n, col) {
  M <- rd(nm(truth, k, n)); if (is.null(M) || !col %in% names(M)) return(NULL)
  p <- M[[col]]; p[!is.finite(p)] <- 1; p
}

rows <- list()
for (n in c(1000L, 5000L)) for (k in as.integer(c(0, 0.001, 0.002, 0.005, 0.01) * n)) {
  for (tn in names(TESTS)) {
    col <- TESTS[[tn]]
    pl <- pv("logit", k, n, col)                       # the matched null at the same contamination
    if (is.null(pl)) next
    cv_contam <- as.numeric(stats::quantile(pl, 0.05, type = 1))
    pl0 <- pv("logit", 0L, n, col)                     # the clean null, for the raw-versus-adjusted view
    cv_clean <- as.numeric(stats::quantile(pl0, 0.05, type = 1))
    for (truth in c("probit", "cloglog")) {
      pt <- pv(truth, k, n, col); if (is.null(pt)) next
      rows[[length(rows) + 1]] <- data.table(
        truth = truth, n = n, k = k, test = tn,
        false_alarm = mean(pl <= 0.05),
        raw = mean(pt <= 0.05),
        adj_clean_null = mean(pt <= cv_clean),
        adj_contaminated_null = mean(pt <= cv_contam),
        crit_contaminated = cv_contam)
    }
  }
}
P <- rbindlist(rows)
P[, loss_vs_clean := adj_clean_null[k == 0L] - adj_contaminated_null, by = .(truth, n, test)]
dir.create(edge_battery("analysis"), showWarnings = FALSE, recursive = TRUE)
fwrite(P, edge_battery("analysis", "contaminated_power.csv"))

cat("\n== power against a real misfit under contamination, three readings ==\n")
cat("   false_alarm = the SAME test on the logistic truth at the same k (its own contaminated null)\n")
cat("   raw = rejection at 0.05; adj = at the critical value of the clean / of the contaminated null\n\n")
for (tr in c("probit", "cloglog")) {
  cat("--", tr, "--\n")
  print(P[truth == tr & test == "EDGE",
          .(n, k, false_alarm = round(false_alarm, 3), raw = round(raw, 3),
            adj_clean = round(adj_clean_null, 3), adj_contam = round(adj_contaminated_null, 3))])
}
cat("\n== the cells Section 6.5 quotes ==\n")
print(P[test == "EDGE" & ((truth == "probit" & n == 5000L & k %in% c(0L, 25L, 50L)) |
                          (truth == "cloglog" & n == 1000L & k %in% c(0L, 10L))),
        .(truth, n, k, false_alarm = round(false_alarm, 3), raw = round(raw, 3),
          adj_contam = round(adj_contaminated_null, 3))])
