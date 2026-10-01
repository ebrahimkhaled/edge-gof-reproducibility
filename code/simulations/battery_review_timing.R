## battery_review_timing.R -- serial replicate timing (review item f). BLAS on one thread, one R process, nothing else
## running. Times whole replicates at n = 1000, 5000, 16,000 (plus 20,000 and 50,000 once), a component breakdown for one
## replicate at each n, compares with battery/_estimate_cells.csv, and rescales the 42.4 CPU-hour projection.
## Writes battery/_review/review_timing.log and review_timing.csv.
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
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))
sink(file.path(OUT, "review_timing.log"), split = TRUE)
cat("battery_review_timing.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
RNGkind("L'Ecuyer-CMRG")
Cells <- battery_cells()
EST <- fread(edge_battery("_estimate_cells.csv"))

TGT <- list(c("1b", "null_base_n1000", 20), c("3", "stk_short_n1000", 20), c("1a", "null_quad_n1000", 20),
            c("3", "cauchit_n5000", 8), c("1a", "null_link_n5000", 8), c("3", "quad_0.05_n5000", 8),
            c("2", "probit_base_n16000", 4), c("2", "t4_base_n20000", 3), c("7", "runD_probit_n50000", 2))
invisible(battery_one(1L, as.list(Cells[Cells$block == "1b" & Cells$cell == "null_base_n610", ])))   # load namespaces
rows <- list()
for (tg in TGT) {
  ce <- as.list(Cells[Cells$block == tg[1] & Cells$cell == tg[2], ]); k <- as.integer(tg[3])
  if (!length(ce$cell)) { cat("missing cell", tg[2], "\n"); next }
  invisible(battery_one(1000L, ce))                                              # warm-up at this n
  el <- system.time(for (r in seq_len(k)) battery_one(r, ce))[["elapsed"]] / k
  est <- EST[block == tg[1] & cell == tg[2]]$sec_per_rep
  rows[[length(rows) + 1]] <- data.frame(block = tg[1], cell = tg[2], n = ce$n, reps = k, sec_per_rep = el,
                                         builder_sec_per_rep = if (length(est)) est else NA, ratio = if (length(est)) el / est else NA)
  cat(sprintf("%-3s %-22s n = %6d  %2d reps  %.3f s/rep  (builder estimate %s)\n", tg[1], tg[2], ce$n, k, el,
              if (length(est)) sprintf("%.3f", est) else "-"))
}
TT <- do.call(rbind, rows)

## component breakdown, one replicate per n
cat("\ncomponents of one replicate (seconds)\n")
comp <- function(ce) {
  set.seed(ce$seed_base + 7L); dat <- bt_data(ce)
  tm <- c()
  tm["data"] <- system.time(dat2 <- { set.seed(ce$seed_base + 7L); bt_data(ce) })[["elapsed"]]
  tm["fit"] <- system.time(fq <- bt_fit(dat))[["elapsed"]]
  arms <- bt_arms(ce)
  tm["edge_all_arms"] <- system.time(for (a in unique(arms)) bt_edge_arm(fq, a, ce))[["elapsed"]]
  tm["stukel_4forms"] <- system.time(bt_stukel(fq))[["elapsed"]]
  tm["giviti_95"] <- system.time(bt_giviti(fq$y, fq$p_raw, 0.95))[["elapsed"]]
  tm["giviti_50"] <- system.time(bt_giviti(fq$y, fq$p_raw, 0.50))[["elapsed"]]
  tm["cubic"] <- system.time(bt_cubic(fq$y, fq$eta))[["elapsed"]]
  if (fq$n < BT_RIVAL_NMAX) {
    tm["hlw"] <- system.time(bt_hlw(fq$y, fq$ph, 10))[["elapsed"]]
    tm["ph"] <- system.time(bt_ph(fq$y, fq$ph, 10))[["elapsed"]]
    tm["tsiatis"] <- system.time(suppressWarnings(bt_tsiatis(fq, 10)))[["elapsed"]]
    tm["xie"] <- system.time(suppressWarnings(bt_xie(fq)))[["elapsed"]]
    tm["pr"] <- system.time(bt_pr(fq))[["elapsed"]]
  }
  tm
}
for (nm in c("null_base_n1000", "cauchit_n5000", "probit_base_n16000", "runD_probit_n50000")) {
  ce <- as.list(Cells[Cells$cell == nm, ][1, ])
  tm <- comp(ce)
  cat(sprintf("  %-20s %s\n", nm, paste(sprintf("%s %.3f", names(tm), tm), collapse = "  ")))
}

## builder's per-cell estimate: how many cells rest on each timed key, and the heaviest cells
cat("\nheaviest 15 cells in _estimate_cells.csv\n")
print(EST[order(-cpu_s)][1:15], row.names = FALSE)
cat("\nCPU-hours by n band in _estimate_cells.csv\n")
EST[, band := cut(n, c(0, 700, 1500, 3000, 7000, 12000, 30000, Inf), right = FALSE)]
print(EST[, .(cells = .N, reps = sum(B), cpu_h = sum(cpu_s, na.rm = TRUE) / 3600, s_per_rep_mean = mean(sec_per_rep, na.rm = TRUE)), by = band][order(band)])

## rescaled projection: per-n ratio of observed to builder seconds, interpolated on log n
TT$logn <- log(TT$n)
rat <- TT[is.finite(TT$ratio), ]
rr <- aggregate(ratio ~ n, data = rat, FUN = mean)
f_ratio <- function(n) stats::approx(log(rr$n), rr$ratio, xout = log(pmax(pmin(n, max(rr$n)), min(rr$n))), rule = 2)$y
EST[, ratio := ifelse(is.na(n), 1, f_ratio(ifelse(is.na(n), 1000, n)))]
cpu_builder <- sum(EST$cpu_s, na.rm = TRUE) / 3600
cpu_rescaled <- sum(EST$cpu_s * EST$ratio, na.rm = TRUE) / 3600
cat(sprintf("\nbuilder projection %.1f CPU-h; rescaled by the observed/builder ratio at the nearest timed n: %.1f CPU-h\n", cpu_builder, cpu_rescaled))
for (w in c(12, 16, 20)) for (sp in c(7, 8, 0.8 * min(w, 12) + 0.25 * max(w - 12, 0))) {
  cat(sprintf("  %2d workers, speed-up %4.1f: builder %.1f wall-h, rescaled %.1f wall-h\n", w, sp, cpu_builder / sp, cpu_rescaled / sp))
}
fwrite(TT, file.path(OUT, "review_timing.csv"))
sink()
