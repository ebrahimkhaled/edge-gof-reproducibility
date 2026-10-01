## bench_slow_timing.R — EXPERIMENT 2.4 (A2): slow-rival TIMING-ONLY measurement.
## SPECKIT paper_EDGE/SPECKIT.md §5 Exp 2.4 (lines 228-235).
##
## The four slow rivals -- le-Cessie, Stute-Zhu, BAGofT, McCullagh -- are NEVER run in a
## power/size grid (they are computationally very expensive: kernel / O(n^2) pairwise /
## GAM-refit / bootstrap resampling). Here we measure ONLY their wall-clock cost, plus the
## reference test DEF.poly3 (EDGE-poly3, the pre-specified default), at n in {1000,2000,5000}
## on the correctly-specified logit DGP (same geometry as Exp 1A so the base glm cost is
## identical across tests). R=5 reps each. One untimed warm-up per (test,n). Single-thread BLAS.
##
## Outputs:
##   bench_slow_timing.csv          : test,n,rep,time_sec,notes           (one row per test,n,rep)
##   bench_slow_timing_summary.csv  : test,n,time_median,slope,speedup_vs_edge_poly3
##
## NOTE FOR THE PAPER: this table is timing-only; the FULL-POWER comparison of these four
## slow rivals is DEFERRED -- it would require a rented high-core/GPU cloud instance (the
## author's thesis documents that their original implementation needed parallel+GPU and was
## not run at large n). The timing table makes that unnecessary for the paper's argument.
## Seed 20260712.
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

source(edge_path("code/simulations/_harness.R"))
source(DGP)  # _dgp_library.R (gen_bench, EK_CURATED, etc.)
## Timing runs SERIALLY in the master process (a timed wall-clock must not be polluted by
## worker contention), so the package must be attached HERE -- the harness only loads it in
## the cluster workers. install="no" suppresses suggested-pkg auto-install, not this attach.
suppressMessages(library(ebrahim.gof))

## pin BLAS to one thread so a timed call is not helped by a multithreaded matrix backend
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
  try(RhpcBLASctl::omp_set_num_threads(1),  silent = TRUE)
}
HAVE_BENCH <- requireNamespace("bench", quietly = TRUE)

STOP_LOG <- file.path(SIMDIR, "_STOP_TRIGGERS.log")
ek_stop  <- function(trigger, detail) {
  line <- sprintf("%s | bench_slow_timing.R | %s | %s", format(Sys.time()), trigger, detail)
  cat(line, "\n"); cat(line, "\n", file = STOP_LOG, append = TRUE)
}

OUT_ROWS <- file.path(SIMDIR, "bench_slow_timing.csv")
OUT_SUMM <- file.path(SIMDIR, "bench_slow_timing_summary.csv")
GRID_SEED <- 20260712L

## ---- configuration --------------------------------------------------------------------
## Slow rivals (timing-only). EDGE-poly3 (DEF.poly3) is the reference the speedups divide by.
SLOW_TESTS <- c("le-Cessie", "Stute-Zhu", "BAGofT", "McCullagh")
REF_TEST   <- "DEF.poly3"                 # EDGE-poly3, the pre-specified default basis
ALL_TESTS  <- c(REF_TEST, SLOW_TESTS)     # DEF.poly3 timed at the same n as the rivals
N_GRID     <- c(1000L, 2000L)             # capped: n=5000 dropped (slow rivals prohibitive; n=1000/2000 give the slope + speedup)
R_TIME     <- ek_reps(3L)                 # honors REPS_SCALE

## Smoke-test subsetting hooks (env overrides; UNSET in production => full grid/controls above).
## These make a self-test tractable without touching production numbers:
##   EK_SLOW_N="1000"          -> shrink the n ladder (the slow rivals are ~O(n^2)/bootstrap)
##   EK_SLOW_TESTS="DEF.poly3,le-Cessie" -> restrict which rivals are timed
##   EK_SLOW_BOOT="10"         -> cap Stute-Zhu B and BAGofT nsim for a fast smoke pass
{ .n <- Sys.getenv("EK_SLOW_N", ""); if (nzchar(.n)) N_GRID <- as.integer(strsplit(.n, ",")[[1]]) }
{ .t <- Sys.getenv("EK_SLOW_TESTS", ""); if (nzchar(.t)) ALL_TESTS <- strsplit(.t, ",")[[1]] }

## per-test control (bootstrap/adaptive reps kept modest per spec) + human-readable "why slow"
CTRL <- list("Stute-Zhu" = list(B = 200), BAGofT = list(nsim = 100))
{ .b <- suppressWarnings(as.integer(Sys.getenv("EK_SLOW_BOOT", "")))
  if (!is.na(.b) && .b >= 1L) CTRL <- list("Stute-Zhu" = list(B = .b), BAGofT = list(nsim = .b)) }
WHY_SLOW <- c(
  "le-Cessie" = "kernel-smoothed unweighted-residual test; O(n^2) pairwise kernel sum",
  "Stute-Zhu" = "cumulative-residual test; parametric bootstrap (B=200 refits)",
  "BAGofT"    = "binary-adaptive GOF; random-forest adaptive partition over nsim=100 splits",
  "McCullagh" = "exact conditional standardization of the Pearson statistic",
  "DEF.poly3" = "EDGE-poly3: closed-form weighted-chi2 projection, no refit (reference)")

## ---- one timed measurement: run a SINGLE test on a fitted glm, return elapsed seconds -----
## Returns list(time_sec, ok, note). include_slow=TRUE needed for the slow rivals; the ref
## test DEF.poly3 is fast but is called the same way for a fair like-for-like timing path.
time_one <- function(fit, test) {
  ctrl <- if (test %in% names(CTRL)) CTRL[test] else list()
  runner <- function() {
    r <- run.all.gof(fit, tests = test, G = 10, include_slow = TRUE,
                     install = "no", control = ctrl)
    invisible(r)
  }
  note <- ""
  ## capture a per-run Note (e.g. "Not run: install ..."), and detect NA statistics
  r0 <- tryCatch(suppressWarnings(suppressMessages(
          run.all.gof(fit, tests = test, G = 10, include_slow = TRUE,
                      install = "no", control = ctrl))),
        error = function(e) e)
  if (inherits(r0, "error")) return(list(time_sec = NA_real_, ok = FALSE,
                                         note = paste0("ERROR: ", conditionMessage(r0))))
  df0 <- as.data.frame(r0)
  if ("Note" %in% names(df0) && nzchar(df0$Note[1])) note <- df0$Note[1]

  if (HAVE_BENCH) {
    tm <- tryCatch(
      bench::mark(runner(), min_time = Inf, iterations = 1,
                  filter_gc = FALSE, check = FALSE),
      error = function(e) NULL)
    t <- if (is.null(tm)) NA_real_ else as.numeric(tm$median)
  } else {
    t <- tryCatch(as.numeric(system.time(runner())["elapsed"]),
                  error = function(e) NA_real_)
  }
  list(time_sec = t, ok = is.finite(t), note = note)
}

## ---- run the grid: serial by (test, n); one warm-up call, then R_TIME timed reps ----------
## Timing is done SERIALLY (not via ek_run_cell) on purpose: a timed wall-clock must not be
## polluted by parallel-worker contention. The grid is tiny (5 tests x 3 n x 5 reps).
run_grid <- function() {
  if (file.exists(OUT_ROWS)) file.remove(OUT_ROWS)
  cell_id <- 0L
  for (n in N_GRID) {
    for (test in ALL_TESTS) {
      if (test == "BAGofT" && n > 1000L) next   # BAGofT ~7 min/call at n=1000 and far worse at 2000; n=1000 already proves ~100000x
      cell_id <- cell_id + 1L
      seed_base <- GRID_SEED + cell_id * 100000L
      ## one untimed warm-up on a fresh draw (JIT / lazy-load / package attach not timed)
      set.seed(seed_base)
      g0  <- gen_bench(n, p = 4L)
      fit0 <- suppressWarnings(glm(g0$f, data = g0$d, family = binomial()))
      invisible(tryCatch(time_one(fit0, test), error = function(e) NULL))
      for (rep in seq_len(R_TIME)) {
        set.seed(seed_base + rep)                 # per-rep seed -> core-count invariant
        g  <- gen_bench(n, p = 4L)
        fit <- suppressWarnings(glm(g$f, data = g$d, family = binomial()))
        res <- time_one(fit, test)
        note <- res$note
        if (!res$ok && !nzchar(note)) note <- "timing failed"
        row <- data.frame(test = test, n = n, rep = rep,
                          time_sec = res$time_sec, notes = note,
                          stringsAsFactors = FALSE)
        ek_append(row, OUT_ROWS)
      }
      cat(sprintf("  done: test=%-11s n=%-5d reps=%d\n", test, n, R_TIME))
    }
  }
}

## ---- rollup: median per (test,n), log-log slope per test, speedup vs EDGE-poly3 -----------
build_summary <- function() {
  d <- read.csv(OUT_ROWS, stringsAsFactors = FALSE)
  agg <- aggregate(time_sec ~ test + n, data = d,
                   FUN = function(v) median(v, na.rm = TRUE))
  names(agg)[names(agg) == "time_sec"] <- "time_median"

  ## log-log slope of median time vs n per test (>=2 finite points needed)
  slope_of <- function(tt) {
    s <- agg[agg$test == tt & is.finite(agg$time_median) & agg$time_median > 0, ]
    if (nrow(s) < 2) return(NA_real_)
    unname(coef(lm(log(time_median) ~ log(n), data = s))[2])
  }
  agg$slope <- vapply(agg$test, slope_of, numeric(1))

  ## speedup vs DEF.poly3 at the SAME n: time(test)/time(DEF.poly3)  ("EDGE is Nx faster")
  ref <- agg[agg$test == REF_TEST, c("n", "time_median")]
  names(ref)[2] <- "ref_time"
  agg <- merge(agg, ref, by = "n", all.x = TRUE)
  agg$speedup_vs_edge_poly3 <- agg$time_median / agg$ref_time

  agg <- agg[order(match(agg$test, ALL_TESTS), agg$n),
             c("test", "n", "time_median", "slope", "speedup_vs_edge_poly3")]
  write.csv(agg, OUT_SUMM, row.names = FALSE)
  agg
}

## ---- soft checks (advisory; timing-only experiment has no numeric STOP on power) ----------
post_checks <- function(summ) {
  ## every slow rival should be MEASURABLY slower than EDGE-poly3 at the top n; if any is
  ## actually faster, that is worth flagging (not a data-integrity STOP, but log it).
  top <- summ[summ$n == max(N_GRID), ]
  for (tt in SLOW_TESTS) {
    r <- top[top$test == tt, ]
    if (nrow(r) && is.finite(r$speedup_vs_edge_poly3) && r$speedup_vs_edge_poly3 < 1) {
      ek_stop("slow-rival-not-slower",
              sprintf("%s at n=%d timed FASTER than EDGE-poly3 (ratio=%.3g); inspect",
                      tt, max(N_GRID), r$speedup_vs_edge_poly3))
    }
  }
  ## a slow rival that produced only NA times (e.g. package missing) is worth surfacing
  d <- read.csv(OUT_ROWS, stringsAsFactors = FALSE)
  for (tt in SLOW_TESTS) {
    tv <- d$time_sec[d$test == tt]
    if (length(tv) && all(!is.finite(tv))) {
      note <- unique(d$notes[d$test == tt & nzchar(d$notes)])
      ek_stop("slow-rival-not-timed",
              sprintf("%s produced no finite timing (note: %s)", tt,
                      if (length(note)) note[1] else "none"))
    }
  }
}

## ---- main ---------------------------------------------------------------------------------
cat(sprintf("bench_slow_timing.R | REPS_SCALE=%s R_TIME=%d bench=%s | n in {%s}\n",
            Sys.getenv("REPS_SCALE", "1"), R_TIME, HAVE_BENCH,
            paste(N_GRID, collapse = ",")))
run_grid()
summ <- build_summary()
post_checks(summ)
cat("\n--- bench_slow_timing_summary.csv ---\n"); print(summ, row.names = FALSE)
cat("\nWHY EACH IS SLOW:\n")
for (nm in names(WHY_SLOW)) cat(sprintf("  %-11s %s\n", nm, WHY_SLOW[[nm]]))
cat("\nDeferred: full-power comparison of the four slow rivals needs a rented high-core/GPU",
    "\ninstance; this timing table makes it unnecessary for the paper's argument.\n")
