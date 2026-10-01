## bench_compute.R — EDGE paper, EXPERIMENT 1 (compute/scalability benchmark) + 1C (Stukel
## separation vs n). SPECKIT §5 Exp 1 / 1A / 1B / 1C. Writes:
##   bench_compute_time.csv          (per n,p,test,rep)
##   bench_compute_time_summary.csv  (median/p10/p90 per n,p,test)
##   bench_compute_scaling_fits.csv  (per-test log-log slope of time vs n + Stukel/EF ratio)
##   bench_stukel_failure.csv        (1C: stukel_fail_rate / ef_na_rate / event_rate vs n)
## STOP triggers (spec): 1B EDGE log-log slope > 1.25 ; 1C EDGE ef_fail_rate > 0.
## Foundation contract: _harness.R (ek_cluster/ek_run_cell/ek_append/ek_reps) + _dgp_library.R
## (gen_bench / gen_sparse_link / stuk_diag / run_curated_full). Package = ebrahim.gof v2.1.1.
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

suppressMessages(library(parallel))
SIMDIR <- edge_path("code/simulations")
source(file.path(SIMDIR, "_harness.R"))   # defines SIMDIR, DGP, EK_NCORES, ek_* helpers
source(DGP)                               # _dgp_library.R (also loaded inside each worker)

## Self-test override: the orchestrator can cap cores to avoid oversubscription.
if (nzchar(Sys.getenv("EK_NCORES_OVERRIDE")))
  EK_NCORES <- max(1L, as.integer(Sys.getenv("EK_NCORES_OVERRIDE")))

## Timing backend availability (probed on the master; the workers probe again inside one_rep_time).
## `bench` gives high-resolution median wall-time + memory allocation; `RhpcBLASctl` pins BLAS to a
## single thread so a timed matrix op is not silently helped by a multithreaded backend. Both are
## OPTIONAL: guarded with requireNamespace so the script still runs (falling back to system.time and
## NA memory) if either is absent.
HAVE_BENCH <- requireNamespace("bench", quietly = TRUE)
HAVE_BLASCTL <- requireNamespace("RhpcBLASctl", quietly = TRUE)

STOPLOG   <- file.path(SIMDIR, "_STOP_TRIGGERS.log")
GRID_SEED <- 20260706L                    # Exp 1 grid seed (spec 1A/1B)
SEED_1C   <- 20260707L                    # Exp 1C grid seed (spec 1C)

ek_stop <- function(script, trigger, detail) {
  line <- sprintf("%s\t%s\t%s\t%s", format(Sys.time()), script, trigger, detail)
  cat("STOP-TRIGGER:", line, "\n")
  cat(line, "\n", file = STOPLOG, append = TRUE)
}

## ---- test sets --------------------------------------------------------------
## Timed in isolation (spec 1A). "glm_fit" = the base glm as reference; le-Cessie is the
## O(n^2) foil, kept ONLY as the quadratic-cost reference and CAPPED at n<=2000 (it is genuinely
## slow -- ~16s at n=2000 -- and its full cost curve is already measured in bench_slow_timing.R, so
## here we just need one or two small-n points to anchor the O(n^2) shape). All others are the fast
## curated/partition rivals.
BENCH_TESTS  <- c("DEF.poly3","DEF.stukel","EF","HL","Pigeon-Heyse","Osius-Rojek",
                  "Stukel","Tsiatis","Xie","Pulkstenis-Robinson","le-Cessie","glm_fit")
SLOW_TESTS   <- c("le-Cessie")            # need include_slow=TRUE in run.all.gof
LECESSIE_CAP <- 2000L                     # le-Cessie only timed at n <= this (was 5000)
LECESSIE_RTIME <- 5L                      # le-Cessie is timed with only this many reps (it is 16s/call)
EDGE_BASES   <- c("DEF.poly3","DEF.stukel")  # STOP-trigger tests for the slope check

## Per-call timing-batch floor (seconds). Production default 0.20; smoke runs may lower it via
## EK_TIME_FLOOR to trade timing precision for speed. Exported to every worker in main().
TIME_FLOOR_VAL <- as.numeric(Sys.getenv("EK_TIME_FLOOR", "0.20"))

## ---- grids ------------------------------------------------------------------
N_LADDER <- c(200,500,1000,2000,5000,10000,20000,50000)   # TOP=50000 (spec A4, do NOT exceed)
P_GRID   <- c(4,8,16,32)                                   # k = p-4 nuisance z1..zk
Rtime_of <- function(n) if (n <= 20000L) 30L else 10L      # spec 1A reps

## Smoke-test subsetting hooks (env-var overrides; unset in production => full grids above).
## e.g. EK_SMOKE_N="200,500" EK_SMOKE_P="4" EK_SMOKE_N1C="100,200" for a fast self-test.
.env_nums <- function(v) { s <- Sys.getenv(v, ""); if (!nzchar(s)) return(NULL)
  as.numeric(strsplit(s, ",")[[1]]) }
if (!is.null(.env_nums("EK_SMOKE_N"))) N_LADDER <- .env_nums("EK_SMOKE_N")
if (!is.null(.env_nums("EK_SMOKE_P"))) P_GRID   <- .env_nums("EK_SMOKE_P")

## =============================================================================
## 1A/1B — WALL-CLOCK SCALING
## =============================================================================
## one_rep: for a given cell=(n,p) it generates ONE bench dataset, fits the base glm once,
## does a warm-up battery call, then times (and, where possible, measures allocated memory for)
## EACH test in isolation. Returns a NAMED numeric vector: time.<test>, mem.<test>, glm_iter,
## refit_iter, refit_converged, n_events, event_rate. Each of the B=R_time reps is one such
## timing pass (parSapply seeds set.seed(seed_base+rep) => reproducible, core-count invariant).
##
## Timing backend: if `bench` is installed we use bench::mark (high-resolution median + allocated
## memory); otherwise we fall back to the system.time micro-batch (per-call estimate resolution-
## limited from below) with NA memory. BLAS is pinned to one thread per worker (via RhpcBLASctl,
## if present) so a timed matrix op is not silently helped by a multithreaded backend. Both
## backends are OPTIONAL and guarded with requireNamespace so this still runs if either is absent.
one_rep_time <- function(rep, cell) {
  n <- cell$n; p <- cell$p

  ## pin BLAS/OMP to a single thread on this worker (guarded; no-op if RhpcBLASctl absent)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1),  silent = TRUE)
  }
  have_bench <- requireNamespace("bench", quietly = TRUE)

  G   <- gen_bench(n, p)
  dat <- G$d; f <- G$f
  ## base fit (timed separately as glm_fit); model=FALSE/y=FALSE to keep memory small
  fit <- suppressWarnings(glm(f, data = dat, family = binomial(), model = TRUE))
  n_events   <- sum(dat$y)
  event_rate <- n_events / n
  glm_iter   <- as.numeric(fit$iter)

  ## Stukel refit diagnostics (augmented-model iter + convergence) via the shared diag helper's
  ## refit; we build the same augmented fit once to record refit_iter/refit_converged.
  e  <- predict(fit); da <- fit$data
  da$za <- 0.5*e^2*(e>=0); da$zb <- -0.5*e^2*(e<0)
  fa <- tryCatch(suppressWarnings(glm(update(formula(fit), . ~ . + za + zb),
                                      data = da, family = binomial())),
                 error = function(er) NULL)
  refit_iter      <- if (is.null(fa)) NA_real_ else as.numeric(fa$iter)
  refit_converged <- if (is.null(fa)) 0 else as.numeric(isTRUE(fa$converged))

  ## ---- per-test timing+memory: returns c(time_sec, mem_bytes) ----
  ## bench::mark path: min_time=Inf + iterations=1 => exactly one measured call, median time +
  ## allocated memory. system.time fallback micro-batches to clear TIME_FLOOR (mem_bytes = NA).
  TIME_FLOOR <- TIME_FLOOR_VAL   # exported to workers from main() (see clusterExport below)
  BATCH_CAP  <- 4096L  # never call a single test more than this many times per timing (fallback only)
  time_call_sys <- function(expr_fun) {
    ## one measured call first
    e1 <- as.numeric(system.time(expr_fun())["elapsed"])
    if (e1 >= TIME_FLOOR || e1 * 2 >= TIME_FLOOR) return(e1)   # already resolvable
    K <- 4L
    repeat {
      el <- as.numeric(system.time(for (i in seq_len(K)) expr_fun())["elapsed"])
      if (el >= TIME_FLOOR || K >= BATCH_CAP) return(el / K)
      ## grow K to project onto TIME_FLOOR (guard against el==0)
      K <- min(BATCH_CAP, as.integer(max(K * 2, ceiling(K * TIME_FLOOR / max(el, 1e-4)))))
    }
  }
  measure_call <- function(expr_fun) {
    if (have_bench) {
      tm <- tryCatch(bench::mark(expr_fun(), min_time = Inf, iterations = 1,
                                 filter_gc = FALSE, check = FALSE),
                     error = function(e) NULL)
      if (is.null(tm)) return(c(NA_real_, NA_real_))
      c(as.numeric(tm$median), as.numeric(tm$mem_alloc))
    } else {
      c(time_call_sys(expr_fun), NA_real_)   # system.time fallback: no memory measurement
    }
  }
  measure_one <- function(tst) {
    if (identical(tst, "glm_fit"))
      return(measure_call(function() suppressWarnings(glm(f, data = dat, family = binomial()))))
    slow <- tst %in% SLOW_TESTS
    measure_call(function() suppressWarnings(
      run.all.gof(fit, tests = tst, include_slow = slow, install = "no")))
  }

  ## which tests are eligible at this n (le-Cessie capped at n<=LECESSIE_CAP, now 2000) AND at
  ## this rep: le-Cessie is ~16s/call, so we only measure it on the first LECESSIE_RTIME reps
  ## (R_time=5, not 30). Reps beyond that leave le-Cessie NA and are dropped in the reshape.
  tests_here <- BENCH_TESTS
  if (n > LECESSIE_CAP)   tests_here <- setdiff(tests_here, "le-Cessie")
  if (rep > LECESSIE_RTIME) tests_here <- setdiff(tests_here, "le-Cessie")

  ## warm-up: one untimed pass of every eligible test (JIT/first-call costs excluded)
  for (tst in tests_here) invisible(try(measure_one(tst), silent = TRUE))

  ## timed+measured pass
  tvec <- setNames(rep(NA_real_, length(BENCH_TESTS)), paste0("time.", BENCH_TESTS))
  mvec <- setNames(rep(NA_real_, length(BENCH_TESTS)), paste0("mem.",  BENCH_TESTS))
  for (tst in tests_here) {
    tm <- measure_one(tst)
    tvec[paste0("time.", tst)] <- tm[1]
    mvec[paste0("mem.",  tst)] <- tm[2]
  }

  c(tvec, mvec,
    glm_iter = glm_iter, refit_iter = refit_iter, refit_converged = refit_converged,
    n_events = n_events, event_rate = event_rate)
}

run_1A <- function(cl) {
  path_long <- file.path(SIMDIR, "bench_compute_time.csv")
  if (file.exists(path_long)) file.remove(path_long)
  cell_id <- 0L
  for (n in N_LADDER) for (p in P_GRID) {
    cell_id <- cell_id + 1L
    B <- ek_reps(Rtime_of(n))
    seed_base <- GRID_SEED + cell_id * 100000L
    cell <- list(n = n, p = p)
    df <- ek_run_cell(cl, B, seed_base, one_rep_time, cell)
    ## reshape wide->long over tests, one row per (n,p,test,rep)
    df$rep  <- seq_len(nrow(df))
    df$seed <- seed_base + df$rep
    rows <- do.call(rbind, lapply(BENCH_TESTS, function(tst) {
      data.frame(n = n, p = p, test = tst, rep = df$rep, seed = df$seed,
                 time_sec = df[[paste0("time.", tst)]],
                 mem_bytes = df[[paste0("mem.", tst)]],
                 glm_iter = df$glm_iter,
                 refit_iter = ifelse(tst == "Stukel", df$refit_iter, NA_real_),
                 refit_converged = ifelse(tst == "Stukel", df$refit_converged, NA_real_),
                 n_events = df$n_events, event_rate = df$event_rate,
                 stringsAsFactors = FALSE)
    }))
    rows <- rows[!is.na(rows$time_sec), , drop = FALSE]   # drop capped le-Cessie cells
    ek_append(rows, path_long)
    cat(sprintf("[1A] n=%d p=%d B=%d done\n", n, p, B))
  }
  path_long
}

## ---- 1A rollup + 1B scaling fits -------------------------------------------
summarise_1A <- function(path_long) {
  d <- read.csv(path_long, stringsAsFactors = FALSE)
  key <- interaction(d$n, d$p, d$test, drop = TRUE)
  agg <- lapply(split(d, key), function(g) {
    data.frame(n = g$n[1], p = g$p[1], test = g$test[1],
               R_time = nrow(g),
               time_median = median(g$time_sec, na.rm = TRUE),
               time_p10 = as.numeric(quantile(g$time_sec, 0.10, na.rm = TRUE)),
               time_p90 = as.numeric(quantile(g$time_sec, 0.90, na.rm = TRUE)),
               mem_median = median(g$mem_bytes, na.rm = TRUE),
               mcse_time = sd(g$time_sec, na.rm = TRUE) / sqrt(nrow(g)),
               stringsAsFactors = FALSE)
  })
  summ <- do.call(rbind, agg)
  summ <- summ[order(summ$test, summ$p, summ$n), ]
  path_summ <- file.path(SIMDIR, "bench_compute_time_summary.csv")
  if (file.exists(path_summ)) file.remove(path_summ)
  ek_append(summ, path_summ)
  cat("[1A] wrote", path_summ, "\n")
  summ
}

scaling_fits_1B <- function(summ, script) {
  ## per-test log-log slope of time_median vs n, pooled over p (spec: identify the slope).
  fits <- lapply(split(summ, summ$test), function(g) {
    g <- g[is.finite(g$time_median) & g$time_median > 0 & g$n > 0, , drop = FALSE]
    if (nrow(g) < 3 || length(unique(g$n)) < 3)
      return(data.frame(test = g$test[1], loglog_slope = NA_real_,
                        n_points = nrow(g), stringsAsFactors = FALSE))
    fit <- lm(log(time_median) ~ log(n), data = g)
    data.frame(test = g$test[1], loglog_slope = unname(coef(fit)[2]),
               n_points = nrow(g), stringsAsFactors = FALSE)
  })
  fitdf <- do.call(rbind, fits)

  ## money column: time(Stukel)/time(EF) at each n (median over p)
  med_np <- function(tst) {
    s <- summ[summ$test == tst, ]
    tapply(s$time_median, s$n, median, na.rm = TRUE)
  }
  st <- med_np("Stukel"); ef <- med_np("EF")
  ns <- sort(unique(summ$n))
  ratio <- data.frame(test = "ratio_Stukel_over_EF", n = ns,
                      time_Stukel = as.numeric(st[as.character(ns)]),
                      time_EF = as.numeric(ef[as.character(ns)]),
                      stringsAsFactors = FALSE)
  ratio$stukel_over_ef <- ratio$time_Stukel / ratio$time_EF

  ## write both blocks (slope table then the ratio-vs-n table) with a common schema-ish shape
  path_fit <- file.path(SIMDIR, "bench_compute_scaling_fits.csv")
  if (file.exists(path_fit)) file.remove(path_fit)
  blockA <- data.frame(kind = "loglog_slope", test = fitdf$test, n = NA_integer_,
                       value = fitdf$loglog_slope, extra = fitdf$n_points,
                       stringsAsFactors = FALSE)
  blockB <- data.frame(kind = "stukel_over_ef", test = "Stukel/EF", n = ratio$n,
                       value = ratio$stukel_over_ef, extra = ratio$time_Stukel,
                       stringsAsFactors = FALSE)
  ek_append(rbind(blockA, blockB), path_fit)
  cat("[1B] wrote", path_fit, "\n")

  ## STOP: any EDGE basis log-log slope > 1.25
  edge <- fitdf[fitdf$test %in% EDGE_BASES & is.finite(fitdf$loglog_slope), ]
  bad  <- edge[edge$loglog_slope > 1.25, ]
  if (nrow(bad) > 0)
    ek_stop(script, "1B_EDGE_loglog_slope_gt_1.25",
            paste(sprintf("%s slope=%.3f", bad$test, bad$loglog_slope), collapse = "; "))
  invisible(fitdf)
}

## =============================================================================
## 1C — STUKEL SEPARATION vs n
## =============================================================================
## Per rep: gen_sparse_link(n,intercept,cov) -> fit glm(...,data=) -> stuk_diag(fit) fail flag,
## and confirm the EDGE tests are computable (run_curated_full non-NA for the EF omnibus and
## the DEF bases). Returns named vector: stukel_fail, stukel_na, ef_na, refit_iter, event_rate.
## `degenerate` flags a draw with FEWER EVENTS THAN GROUPS (n_events < G): no grouped GOF test
## of any kind can bin such a sample, so it is uninformative rather than an EDGE bug. It is
## recorded and reported separately (edge_degenerate_rate) and EXCLUDED from ef_fail_rate, which
## the STOP trigger reserves for EDGE failing on a *computable* (events >= G) sparse sample.
one_rep_1C <- function(rep, cell) {
  n <- cell$n; intercept <- cell$intercept; cov <- cell$cov; Gk <- 10L
  G   <- gen_sparse_link(n, intercept = intercept, beta = 1.0, cov = cov)
  dat <- G$d; f <- G$f
  ne  <- sum(dat$y)
  out <- c(stukel_fail = NA_real_, stukel_na = 1, ef_na = NA_real_, degenerate = 0,
           refit_iter = NA_real_, event_rate = mean(dat$y), n_events = ne)
  ## sub-G events (or all-events): no grouped test can compute -> degenerate, not an EDGE failure
  if (ne < Gk || (n - ne) < Gk || length(unique(dat$y)) < 2L) { out["degenerate"] <- 1; return(out) }
  fit <- tryCatch(suppressWarnings(glm(f, data = dat, family = binomial())),
                  error = function(e) NULL)
  if (is.null(fit)) { out["degenerate"] <- 1; return(out) }  # glm itself couldn't fit
  ## Stukel refit-based separation diagnostic (fit MUST carry data=, which glm(...,data=) does)
  sd <- tryCatch(stuk_diag(fit), error = function(e) c(p = NA, fail = 1))
  out["stukel_fail"] <- as.numeric(sd["fail"])
  out["stukel_na"]   <- as.numeric(is.na(sd["p"]))
  ## augmented refit iterations (for refit_iter_median reporting)
  e <- predict(fit); da <- fit$data; da$za <- 0.5*e^2*(e>=0); da$zb <- -0.5*e^2*(e<0)
  fa <- tryCatch(suppressWarnings(glm(update(formula(fit), . ~ . + za + zb),
                                      data = da, family = binomial())), error = function(er) NULL)
  out["refit_iter"] <- if (is.null(fa)) NA_real_ else as.numeric(fa$iter)
  ## EDGE must never fail on a computable sample: run the curated battery, check EF + DEF non-NA
  cur <- tryCatch(run_curated_full(fit, G = Gk), error = function(e) NULL)
  if (is.null(cur)) {
    out["ef_na"] <- 1
  } else {
    edge_p <- cur[c("p.EF","p.DEF.poly2","p.DEF.poly3","p.DEF.stukel")]
    out["ef_na"] <- as.numeric(any(is.na(edge_p)))
  }
  out
}

## Summarise a single 1C cell into one row of bench_stukel_failure.csv. Rates for Stukel and
## EDGE are over the COMPUTABLE (non-degenerate, events>=G) subset; the degenerate subset (fewer
## events than groups) is reported separately as edge_degenerate_rate and never counted as an
## EDGE failure (no grouped test can bin such a sample).
summ_1C_cell <- function(df, design, n, B, intercept, covariate_df) {
  comp <- df$degenerate == 0                       # computable draws
  ncomp <- sum(comp)
  fail_rate <- if (ncomp) mean(df$stukel_fail[comp], na.rm = TRUE) else NA_real_
  ef_fail   <- if (ncomp) mean(df$ef_na[comp],       na.rm = TRUE) else 0   # EDGE fail on computable
  data.frame(
    design = design, n = n, B = B,
    stukel_fail_rate = fail_rate,
    stukel_na_rate   = if (ncomp) mean(df$stukel_na[comp], na.rm = TRUE) else NA_real_,
    ef_na_rate       = ef_fail,
    ef_fail_rate     = ef_fail,                     # STOP trigger watches this (computable only)
    edge_degenerate_rate = mean(df$degenerate),     # draws with fewer events than groups
    n_computable     = ncomp,
    refit_iter_median = if (ncomp) median(df$refit_iter[comp], na.rm = TRUE) else NA_real_,
    event_rate_mean  = mean(df$event_rate, na.rm = TRUE),
    mcse_fail        = sqrt(pmax(0, fail_rate*(1-fail_rate)) / max(1, ncomp)),
    intercept        = intercept,
    covariate_df     = covariate_df,
    stringsAsFactors = FALSE)
}

run_1C <- function(cl, script) {
  path_1c <- file.path(SIMDIR, "bench_stukel_failure.csv")
  if (file.exists(path_1c)) file.remove(path_1c)
  B_fail <- ek_reps(5000L)

  ## primary sweep: n x {sparse(intercept -4.9), moderate(intercept -1)}, chisq4 covariate
  N_1C <- c(100,150,200,300,500,1000,2000,5000)
  if (!is.null(.env_nums("EK_SMOKE_N1C"))) N_1C <- .env_nums("EK_SMOKE_N1C")
  designs <- list(
    list(design = "sparse",   intercept = -4.9, cov = "chisq4", cov_df = 4),
    list(design = "moderate", intercept = -1.0, cov = "chisq4", cov_df = 4)
  )
  cell_id <- 0L
  sparse_fail_at_small <- c()   # track sparse fail rate at n<=300 for the fallback trigger
  for (dg in designs) for (n in N_1C) {
    cell_id <- cell_id + 1L
    seed_base <- SEED_1C + cell_id * 100000L
    cell <- list(n = n, intercept = dg$intercept, cov = dg$cov)
    df <- ek_run_cell(cl, B_fail, seed_base, one_rep_1C, cell)
    row <- summ_1C_cell(df, dg$design, n, B_fail, dg$intercept, dg$cov_df)
    ek_append(row, path_1c)
    cat(sprintf("[1C] %s n=%d fail=%.3f ef_na=%.3f ev=%.3f\n",
                dg$design, n, row$stukel_fail_rate, row$ef_na_rate, row$event_rate_mean))
    ## STOP: EDGE ef_fail_rate > 0 (EDGE must never fail to compute)
    if (isTRUE(row$ef_fail_rate > 0))
      ek_stop(script, "1C_EDGE_ef_fail_rate_gt_0",
              sprintf("%s n=%d ef_fail_rate=%.4f", dg$design, n, row$ef_fail_rate))
    if (dg$design == "sparse" && n <= 300) sparse_fail_at_small <- c(sparse_fail_at_small, row$stukel_fail_rate)
  }

  ## FALLBACK LADDER (spec 1C): if sparse stukel_fail_rate ~0 at n<=300, search harsher
  ## intercepts / smaller n / chisq1 covariate; report the first cell with fail_rate>=0.05.
  do_fallback <- is.null(.env_nums("EK_SMOKE_N1C"))   # skip the ladder in smoke subsets
  if (do_fallback && length(sparse_fail_at_small) && max(sparse_fail_at_small, na.rm = TRUE) < 0.05) {
    cat("[1C] fallback ladder engaged (sparse fail ~0 at n<=300)\n")
    ladder <- list(
      list(design = "sparse_i5.5",  intercept = -5.5, cov = "chisq4", cov_df = 4, ns = c(100,150,200,300)),
      list(design = "sparse_i6.2",  intercept = -6.2, cov = "chisq4", cov_df = 4, ns = c(100,150,200,300)),
      list(design = "sparse_i6.2_x1", intercept = -6.2, cov = "chisq1", cov_df = 1, ns = c(60,80,100,150,200,300)),
      list(design = "sparse_i6.2_small", intercept = -6.2, cov = "chisq4", cov_df = 4, ns = c(60,80,100))
    )
    found <- FALSE
    for (lg in ladder) for (n in lg$ns) {
      cell_id <- cell_id + 1L
      seed_base <- SEED_1C + cell_id * 100000L
      cell <- list(n = n, intercept = lg$intercept, cov = lg$cov)
      df <- ek_run_cell(cl, B_fail, seed_base, one_rep_1C, cell)
      row <- summ_1C_cell(df, lg$design, n, B_fail, lg$intercept, lg$cov_df)
      ek_append(row, path_1c)
      cat(sprintf("[1C-fallback] %s n=%d fail=%.3f ef_na=%.3f ev=%.3f\n",
                  lg$design, n, row$stukel_fail_rate, row$ef_na_rate, row$event_rate_mean))
      if (isTRUE(row$ef_fail_rate > 0))
        ek_stop(script, "1C_EDGE_ef_fail_rate_gt_0",
                sprintf("%s n=%d ef_fail_rate=%.4f", lg$design, n, row$ef_fail_rate))
      if (!found && isTRUE(row$stukel_fail_rate >= 0.05)) {
        found <- TRUE
        cat(sprintf("[1C] FIRST cell with stukel_fail_rate>=0.05: %s n=%d fail=%.3f\n",
                    lg$design, n, row$stukel_fail_rate))
      }
    }
    if (!found)
      cat("[1C] NOTE: no ladder cell reached stukel_fail_rate>=0.05; ",
          "claim must be downgraded per spec CLAIM-CALIBRATION RULE.\n", sep = "")
  }
  path_1c
}

## =============================================================================
## DRIVER
## =============================================================================
main <- function() {
  script <- "bench_compute.R"
  cat(sprintf("bench_compute.R start | EK_NCORES=%d | REPS_SCALE=%s\n",
              EK_NCORES, Sys.getenv("REPS_SCALE", "1")))
  do_1A <- Sys.getenv("EK_DO_1A", "1") == "1"
  do_1C <- Sys.getenv("EK_DO_1C", "1") == "1"

  cl <- ek_cluster(GRID_SEED)             # Pattern A: one long-lived cluster (spec 24-core plan)
  on.exit(stopCluster(cl), add = TRUE)
  ## push the script-level constants the one_rep_* closures reference onto every worker
  clusterExport(cl, c("BENCH_TESTS","SLOW_TESTS","LECESSIE_CAP","LECESSIE_RTIME","TIME_FLOOR_VAL"),
                envir = environment())

  if (do_1A) {
    path_long <- run_1A(cl)
    summ <- summarise_1A(path_long)
    scaling_fits_1B(summ, script)
  }
  if (do_1C) run_1C(cl, script)

  cat("bench_compute.R done\n")
}

if (identical(environment(), globalenv()) && !interactive()) main()
