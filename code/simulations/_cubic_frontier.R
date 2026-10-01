## _cubic_frontier.R -- the cheap closed-form rivals on the scenarios Section 5.7 uses for its cost frontier,
## and on the asymmetric-link family. Size-adjusted at each test's own matched null, from the deposited
## per-replicate p-values, because the study's analysis files pair EDGE only with the partition tests.
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
suppressPackageStartupMessages(library(data.table))
SIM <- edge_path("code/simulations")
setwd(SIM)
BAT <- edge_battery()
source(file.path(SIM, "_battery_cells.R"))
C <- battery_cells()

TESTS <- c("EDGE.poly3.u.Grule", "Cubic.LR", "Stk.joint", "GiViTI", "HL.Grule")
RD <- new.env(parent = emptyenv())
pv <- function(cell) {
  if (!is.null(RD[[cell]])) return(RD[[cell]])
  row <- C[C$cell == cell, ]
  f <- edge_battery(row$block[1], sprintf("%s_pvalues.csv.gz", cell))
  if (!file.exists(f)) return(NULL)
  RD[[cell]] <- fread(f)
  RD[[cell]]
}
crit <- function(null_cell, test) {
  M <- pv(null_cell); if (is.null(M) || !test %in% names(M)) return(NA_real_)
  p <- M[[test]]; p[!is.finite(p)] <- 1
  as.numeric(stats::quantile(p, 0.05, type = 1))
}
adj_power <- function(cell, test) {
  M <- pv(cell); if (is.null(M) || !test %in% names(M)) return(NA_real_)
  nc <- C$null_cell[C$cell == cell][1]
  cv <- crit(nc, test); if (!is.finite(cv)) return(NA_real_)
  p <- M[[test]]; p[!is.finite(p)] <- 1
  mean(p <= cv)
}

## the twenty scenarios of the cost frontier: the alternatives block 8 ran the resampling rivals on
S8 <- fread(edge_battery("8", "_summary.csv"))
frontier <- S8[role == "alternative" & test == "proj", unique(cell)]
Fr <- rbindlist(lapply(frontier, function(cl)
  data.table(cell = cl, test = TESTS, power = vapply(TESTS, function(t) adj_power(cl, t), 0))))
cat("== cost frontier: mean size-adjusted power over the", length(frontier), "scenarios ==\n")
print(Fr[, .(scenarios = sum(is.finite(power)), mean = round(mean(power, na.rm = TRUE), 3)), by = test][order(-mean)])

## the asymmetric-link family (family 2 of the declaration)
Mem <- unique(fread(edge_battery("analysis", "rule_A_membership.csv"))[, .(cell, family)], by = "cell")
fam2 <- Mem$cell[Mem$family == 2]
A2 <- rbindlist(lapply(fam2, function(cl)
  data.table(cell = cl, test = TESTS, power = vapply(TESTS, function(t) adj_power(cl, t), 0))))
W <- dcast(A2, cell ~ test, value.var = "power")
W <- W[is.finite(EDGE.poly3.u.Grule) & is.finite(Cubic.LR)]
cat("\n== asymmetric links: the cubic calibration LR against the directed test ==\n")
cat("scenarios:", nrow(W), "\n")
cat("mean power  EDGE", round(mean(W$EDGE.poly3.u.Grule), 3), " cubic", round(mean(W$Cubic.LR), 3), "\n")
cat("cubic ahead in", sum(W$Cubic.LR > W$EDGE.poly3.u.Grule), "scenarios; largest margin",
    round(max(W$Cubic.LR - W$EDGE.poly3.u.Grule), 3), "\n")
cat("EDGE leads or ties (within 0.01) in", sum(W$EDGE.poly3.u.Grule - W$Cubic.LR > -0.01), "\n")
print(W[order(Cubic.LR - EDGE.poly3.u.Grule)][, .(cell, EDGE = round(EDGE.poly3.u.Grule, 3),
                                                  cubic = round(Cubic.LR, 3),
                                                  margin = round(Cubic.LR - EDGE.poly3.u.Grule, 3))])
