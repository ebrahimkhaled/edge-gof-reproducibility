## run_M_blockEXT_SP.R -- block EXT_SP: the cumulative calibration test of Sadatsafavi and Petkau (2024) on block EXT.
## Contract: declarations/PREDECLARATION_blockEXT_SP.md (sha256 7455de2d...), hashed before any replicate ran.
##
## Every one of block EXT's 21 cells is regenerated with its own seeds, checked against EXT's stored events and
## Spiegelhalter p-values, and given cumulcalib(y, p, method = "BB"). Each chunk of replicates goes to disk as soon
## as it finishes (parts/), so a restart picks up where the last run stopped; a cell's final file is assembled from
## its parts once all replicates are in.
##
##   Rscript run_M_blockEXT_SP.R [--workers 14] [--reps 1000] [--chunk 100] [--cells 1,4,14] [--smoke] [--timing]
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
suppressPackageStartupMessages({ library(parallel); library(data.table) })
args <- commandArgs(TRUE)
opt <- function(k, d) { i <- match(k, args); if (is.na(i)) d else args[i + 1] }
SMOKE <- "--smoke" %in% args
W <- as.integer(opt("--workers", 14))
REPS <- if (SMOKE) 3L else as.integer(opt("--reps", 1000))
CHUNK <- as.integer(opt("--chunk", 100))
OUT <- edge_battery(if (SMOKE) "EXT_SP_smoke" else "EXT_SP")
PARTS <- file.path(OUT, "parts")
dir.create(PARTS, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "_progress.log")
lg <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
  cat(s, "\n"); cat(s, "\n", file = LOG, append = TRUE) }

## block EXT's cells, as in run_M_blockEXT.R
CELLS <- rbind(
  data.frame(fam = "null", C = 0, n = c(250L, 500L, 1000L), k = 0L),
  data.frame(fam = rep(c("large", "slope", "ushape", "thresh", "inter"), each = 2),
             C = c(.2, .4, .8, .6, .4, .8, 2, 4, .6, 1), n = 500L, k = 0L),
  data.frame(fam = rep(c("exag", "sign"), each = 4), C = rep(c(4, -4), each = 4), n = 1000L,
             k = rep(c(1L, 2L, 5L, 10L), 2)))
CELLS$id <- seq_len(nrow(CELLS))
pick <- opt("--cells", NA)
if (!is.na(pick)) CELLS <- CELLS[CELLS$id %in% as.integer(strsplit(pick, ",")[[1]]), ]

## block EXT's data, drawn from the replicate's own seed with the RNG kind fixed every time
ext_data <- function(cc, r) {
  n <- cc$n; b0 <- c(-0.5, 0.5, 0.8, 0.6, 0.4, 0.3)
  set.seed(if (cc$id <= 13) 7000000L + cc$id * 10000L + r else 7600000L + cc$id * 10000L + r,
           kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  X <- matrix(stats::rnorm(n * 5), n, 5)
  eta_true <- as.numeric(cbind(1, X) %*% b0)
  if (cc$id <= 13) {
    eta0 <- eta_true
    d <- switch(cc$fam, null = 0, large = cc$C, slope = (cc$C - 1) * eta0,
                ushape = cc$C * (X[, 1]^2 - 1), thresh = cc$C * pmax(X[, 2] - 1, 0),
                inter = cc$C * X[, 1] * X[, 3])
    y <- stats::rbinom(n, 1L, stats::plogis(eta0 + d))
  } else {
    y <- stats::rbinom(n, 1L, stats::plogis(eta_true))
    bad <- sample.int(n, cc$k)
    Xr <- X; Xr[bad, 2] <- cc$C * X[bad, 2]
    eta0 <- as.numeric(cbind(1, Xr) %*% b0)
  }
  list(y = y, p = stats::plogis(eta0))
}

one_rep <- function(r, cc) {
  D <- ext_data(cc, r); y <- D$y; p <- D$p; v <- p * (1 - p)
  zsp <- sum((y - p) * (1 - 2 * p)) / sqrt(sum((1 - 2 * p)^2 * v))     # the check against block EXT
  msg <- character(0)
  sp <- tryCatch(withCallingHandlers(cumulcalib::cumulcalib(y, p, method = c("BB", "BM")),
                                     message = function(m) { msg <<- c(msg, conditionMessage(m)); invokeRestart("muffleMessage") },
                                     warning = function(w) { msg <<- c(msg, conditionMessage(w)); invokeRestart("muffleWarning") }),
                 error = function(e) NULL)
  bb <- if (is.null(sp)) rep(NA_real_, 3) else unname(c(sp$by_method$BB$pval, sp$by_method$BB$pval_by_component))
  c(id = cc$id, rep = r, events = sum(y), spz_check = 2 * stats::pnorm(-abs(zsp)),
    SP_BB = bb[1], SP_BB_mean = bb[2], SP_BB_bridge = bb[3],
    SP_BM = if (is.null(sp)) NA_real_ else sp$by_method$BM$pval, notes = length(msg))
}

cell_name <- function(cc) sprintf("c%02d_%s_C%s_n%04d_k%02d", cc$id, cc$fam, cc$C, cc$n, cc$k)
ext_ref <- function(cc) fread(edge_battery("EXT", paste0(cell_name(cc), "_pvalues.csv.gz")))

## ---------------------------------------------------------------- timing on one core
if ("--timing" %in% args) {
  source(file.path(SIMDIR, "edge_external.R"))
  ## the Windows clock ticks in about 15 ms, so each timed run repeats the call m times (inner_loop) and divides
  tm <- function(f, m, k = 5) { f(); s <- vapply(seq_len(k), function(i) system.time(for (j in seq_len(m)) f())[["elapsed"]] / m, 0)
    c(stats::median(s), min(s)) }
  TT <- rbindlist(lapply(c(1000L, 5000L, 20000L, 100000L), function(n) {
    D <- ext_data(list(id = 1L, fam = "null", C = 0, n = n, k = 0L), 1L)
    m <- max(5L, as.integer(200000L / n))
    a <- tm(function() cumulcalib::cumulcalib(D$y, D$p, method = "BB"), m)
    b <- tm(function() edge_external(D$y, D$p, G = 10), m)
    data.table(test = c("Sadatsafavi-Petkau BB (cumulcalib)", "EDGE external, ten groups"), n = n,
               median_sec = c(a[1], b[1]), min_sec = c(a[2], b[2]), reps = 5L, inner_loop = m)
  }))
  fwrite(TT, edge_battery("T", "timing_sp.csv")); print(TT)
  quit(save = "no")
}

## ---------------------------------------------------------------- the replicates
lg("block EXT_SP: %d cells x %d replicates, chunks of %d, %d workers%s", nrow(CELLS), REPS, CHUNK, W,
   if (SMOKE) " [SMOKE]" else "")
cl <- makePSOCKcluster(W)
clusterExport(cl, c("ext_data", "one_rep"))
invisible(clusterEvalQ(cl, { loadNamespace("cumulcalib"); NULL }))
for (i in seq_len(nrow(CELLS))) {
  cc <- CELLS[i, ]; nm <- cell_name(cc)
  f <- file.path(OUT, paste0(nm, "_pvalues.csv.gz"))
  if (file.exists(f)) { lg("cell %d present, skipped", cc$id); next }
  ref <- ext_ref(cc)
  t0 <- proc.time()[["elapsed"]]
  starts <- seq(1L, REPS, by = CHUNK)
  for (s in starts) {
    e <- min(s + CHUNK - 1L, REPS)
    pf <- file.path(PARTS, sprintf("%s_r%04d-%04d.csv.gz", nm, s, e))
    if (file.exists(pf)) next                                        # resume: this chunk is already on disk
    R <- rbindlist(lapply(parLapplyLB(cl, s:e, one_rep, cc = cc), function(v) as.data.table(as.list(v))))
    m <- ref[match(R$rep, ref$rep)]
    if (anyNA(m$rep) || !all(R$events == m$events) || !isTRUE(all.equal(R$spz_check, m$spiegelhalter)))
      stop(sprintf("block EXT_SP: cell %d replicates %d-%d differ from block EXT's data sets", cc$id, s, e))
    tmp <- sub("\\.csv\\.gz$", ".part.csv.gz", pf); fwrite(R, tmp); file.rename(tmp, pf)
  }
  P <- rbindlist(lapply(file.path(PARTS, sprintf("%s_r%04d-%04d.csv.gz", nm, starts, pmin(starts + CHUNK - 1L, REPS))), fread))
  setorder(P, rep)
  stopifnot(nrow(P) == REPS, identical(as.integer(P$rep), seq_len(REPS)))
  tmp <- sub("\\.csv\\.gz$", ".part.csv.gz", f); fwrite(P, tmp); file.rename(tmp, f)
  lg("cell %2d %-6s C=%-4s n=%4d k=%2d | SP_BB %.3f  SP_BM %.3f  (EXT EDGE_G10 %.3f) | notes %d | data checked | %.1f min",
     cc$id, cc$fam, cc$C, cc$n, cc$k, mean(!is.na(P$SP_BB) & P$SP_BB < .05), mean(!is.na(P$SP_BM) & P$SP_BM < .05),
     mean(ref$EDGE_G10[seq_len(REPS)] < .05), sum(P$notes), (proc.time()[["elapsed"]] - t0) / 60)
}
stopCluster(cl)
lg("block EXT_SP: run finished")
