## run_M_battery.R -- driver of the EDGE restructure battery (paper_EDGE/theory/PREDECLARATION_restructure_battery.md, E3-E4).
##
## Usage
##   Rscript run_M_battery.R --block 0 [--B 200] [--workers 20] [--worker-grid 12,16,20]
##   Rscript run_M_battery.R --block 3 [--cells a,b,...] [--workers 16]
##   Rscript run_M_battery.R --block 3 --B 10 --workers 8 --root battery/dryrun      (a test run)
##   Rscript run_M_battery.R --estimate [--block all] [--workers 16] [--root battery/_test2]   (CPU-hours; efficiency from <root>/0)
##   Rscript run_M_battery.R --summary --block 3                                       (rebuild the block summary)
##   Rscript run_M_battery.R --list                                                    (write battery/cells.csv)
##
## Output (E4): battery/<block>/<cell>_pvalues.csv.gz, one row per replicate (seed, n, events, flags, every p-value);
## battery/<block>/_progress.log; battery/<block>/_summary.csv. A cell whose output exists is skipped, so a block can be
## restarted after an interruption. --B marks a test run: its outputs go to battery/_test unless --root is given, and a
## missing block 0 pass is then only a warning.
##
## Block 0 (E3 row 0): (i) old-seed identity with the stored July p-values (EDGE poly2/poly3/stk unit at G = 10, HL, HL_F,
## to 1e-6); (ii) harness against the installed ebrahim.gof >= 2.8.0 (def.gof unit and score, gof_stukel joint, and the
## ported rivals, to 1e-8) on every replicate except samples with no event or no non-event, plus six such constructed samples
## and a behaviour probe of the package (R2 construction, a no-event sample, a flat fit); (iii) worker timing. It writes
## battery/RULE_weighting.txt (Section A verbatim) and exits 1 unless (i) and (ii) pass, 2 if the package is older than
## 2.8.0 or fails the probe. Blocks 1a-7 refuse to start before block 0 has passed; block 2 also refuses without
## RULE_weighting.txt.
##
## Block 6 real one-off fits (E8.5): the output has one row for the stored row order (rep 0) and one per random row order
## (rep k, seed 20260914 + k); battery/6/_row_orders.csv gives per test the stored value, the median and the 5%-95% range.
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
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))

PKG_MIN <- "2.8.0"
BAT     <- edge_battery()
RULE    <- edge_battery("RULE_weighting.txt")
PRE     <- file.path(SIMDIR, "..", "paper_EDGE", "theory", "PREDECLARATION_restructure_battery.md")
SHA_AD  <- "aad03b19c20d393739191f47b587b072708e8e2fef5d497cbad914a9e594f2b7"   # sections A-D, frozen
SHA_AE4 <- "554e43b66bb4faac95195db7491303365f8445c02968edbaad0ee2523edc2e06"   # sections A-E4, written before any cell ran
ALPHAS  <- c(0.01, 0.05, 0.10)
META    <- c("rep", "seed", "n", "events", "fit_ok", "glm_conv")
RIVALS  <- c("HL_w", "PH", "Tsiatis", "Xie", "PR")

## ---- options ----------------------------------------------------------------------------------------------------------
parse_opts <- function(a) {
  opt <- list(block = NULL, cells = NULL, B = NULL, workers = NULL, root = NULL, "worker-grid" = "12,16,20",
              estimate = FALSE, summary = FALSE, list = FALSE)
  flags <- c("estimate", "summary", "list")
  i <- 1L
  while (i <= length(a)) {
    key <- sub("^--", "", a[i]); val <- NULL
    if (grepl("=", key)) { val <- sub("^[^=]*=", "", key); key <- sub("=.*$", "", key) }
    if (!key %in% names(opt)) stop("unknown option --", key)
    if (key %in% flags) { opt[[key]] <- TRUE; i <- i + 1L; next }
    if (is.null(val)) { i <- i + 1L; val <- a[i] }
    opt[[key]] <- val; i <- i + 1L
  }
  opt
}
opt <- parse_opts(commandArgs(trailingOnly = TRUE))
TEST_MODE <- !is.null(opt$B)
ROOT <- if (!is.null(opt$root)) opt$root else if (TEST_MODE) edge_battery("_test") else BAT
if (!grepl("^([A-Za-z]:|/)", ROOT)) ROOT <- file.path(SIMDIR, ROOT)
ROOT <- normalizePath(ROOT, winslash = "/", mustWork = FALSE)
if (!startsWith(ROOT, normalizePath(BAT, winslash = "/", mustWork = FALSE)))
  stop("--root must lie inside ", BAT)
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)

pkg_version <- function() tryCatch(suppressWarnings(utils::packageDescription("ebrahim.gof")$Version), error = function(e) NA_character_)
pkg_ok <- function() { v <- pkg_version(); !is.na(v) && utils::compareVersion(v, PKG_MIN) >= 0 }
default_workers <- function() {
  f <- file.path(ROOT, "0", "workers_choice.txt")
  if (file.exists(f)) return(as.integer(readLines(f, warn = FALSE)[1]))
  max(1L, min(20L, detectCores() - 2L))
}
W <- if (!is.null(opt$workers)) as.integer(opt$workers) else default_workers()

## ---- the weighting rule (Section A) ---------------------------------------------------------------------------------------
sha256_bytes <- function(r) digest::digest(r, algo = "sha256", serialize = FALSE)
rule_bytes <- function() {
  b <- readBin(PRE, "raw", file.info(PRE)$size)
  iE <- grepRaw("## E. Cell list", b, fixed = TRUE)
  if (!length(iE) || sha256_bytes(b[seq_len(iE - 1)]) != SHA_AD) stop("sections A-D of the pre-declaration no longer match their sha256")
  iE7 <- grepRaw("### E7. Sample sizes", b, fixed = TRUE)
  if (length(iE7) && sha256_bytes(b[seq_len(iE7 - 1)]) != SHA_AE4) warning("sections A-E4 differ from the sha256 recorded before launch")
  iA <- grepRaw("## A. Choosing", b, fixed = TRUE); iB <- grepRaw("## B. Tests", b, fixed = TRUE)
  b[iA:(iB - 1)]
}
write_rule <- function() {
  r <- rule_bytes()
  if (file.exists(RULE)) {
    if (!identical(readBin(RULE, "raw", file.info(RULE)$size), r)) stop("battery/RULE_weighting.txt is not Section A verbatim")
  } else {
    dir.create(BAT, showWarnings = FALSE); writeBin(r, RULE)
  }
  cat("RULE_weighting.txt sha256", digest::digest(file = RULE, algo = "sha256"), "\n")
}
check_rule <- function() {
  if (!file.exists(RULE)) stop("refusing to start block 2: battery/RULE_weighting.txt is missing (run block 0 first)")
  if (!identical(readBin(RULE, "raw", file.info(RULE)$size), rule_bytes())) stop("refusing to start block 2: RULE_weighting.txt is not Section A verbatim")
  cat("RULE_weighting.txt sha256", digest::digest(file = RULE, algo = "sha256"), "\n")
}

## ---- workers --------------------------------------------------------------------------------------------------------------
pin_blas <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
}
bt_cluster <- function(w, load_pkg = FALSE) {
  pin_blas()
  ## one retry on another port: a socket clash with a cluster that another R session is opening at the same moment
  ## stopped one dry-run block; the port is taken from the process id so no random number is drawn in the master
  cl <- tryCatch(makeCluster(w), error = function(e) {
    message("makeCluster failed (", conditionMessage(e), "); retrying once on another port")
    Sys.sleep(3); makeCluster(w, port = 11000L + (Sys.getpid() %% 997L))
  })
  parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
  clusterCall(cl, function(f, pkg) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
    source(f); RNGkind("L'Ecuyer-CMRG")
    if (pkg) suppressPackageStartupMessages(library(ebrahim.gof))
    NULL
  }, file.path(SIMDIR, "_battery_tests.R"), load_pkg)
  cl
}
run_cell <- function(cl, cell, B) {
  chunks <- splitIndices(B, min(B, length(cl) * 4L))
  if (cell$generator == "real") chunks <- lapply(chunks, function(i) i - 1L)       # rep 0 = stored row order, 1 to B - 1 = random orders
  M <- do.call(rbind, parLapplyLB(cl, chunks, bt_chunk, cell = cell))
  as.data.frame(M[order(M[, "rep"]), , drop = FALSE])
}
write_cell <- function(M, path) {
  tmp <- paste0(path, ".part")
  fwrite(M, tmp, compress = "gzip")
  if (file.exists(path)) file.remove(path)
  file.rename(tmp, path)
}
plog <- function(block, msg) {
  line <- sprintf("%s  %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", sep = "")
  cat(line, "\n", file = file.path(ROOT, block, "_progress.log"), append = TRUE, sep = "")
}
## rough seconds per replicate on one core, for the ETA only (battery/_estimate_cells.csv replaces it when present)
cost_guess <- function(cells) {
  f <- edge_battery("_estimate_cells.csv")
  n <- ifelse(is.na(cells$n), 1000, cells$n)
  g <- (0.02 + 7e-5 * n) * ifelse(cells$type == "external", 0.3, ifelse(cells$type == "edgeonly", 0.25, 1))
  if (file.exists(f)) {
    e <- fread(f); m <- match(paste(cells$block, cells$cell), paste(e$block, e$cell))
    g[!is.na(m)] <- e$sec_per_rep[m[!is.na(m)]]
  }
  g
}

## ---- stored July p-values (block 0 (i)) and the projection pairing (block 4) ----------------------------------------------------
stored_pvalues <- function(cell, B) {
  cache <- edge_battery("_cache"); dir.create(cache, showWarnings = FALSE)
  f <- file.path(cache, sprintf("stored_%s.csv", cell$cell))
  if (file.exists(f)) { S <- fread(f); if (max(S$rep) >= B) return(S[rep <= B]) }
  Bc <- max(B, 200L)
  if (cell$stored == "sim_null") {
    S <- fread(file.path(SIMDIR, "sim_null_pvalues.csv"))[family == cell$family & n == cell$n & G == 10L & rep <= Bc]
    S <- S[, c("rep", "seed", grep("^p\\.", names(S), value = TRUE)), with = FALSE]
    setnames(S, sub("^p\\.", "", names(S)))
  } else {
    S <- if (cell$stored == "sim_power_broad")
      fread(file.path(SIMDIR, "sim_power_broad_pvalues.csv"), select = c("family", "param", "n", "G", "rep", "seed", "test", "p_value"),
            colClasses = list(character = "param"))[family == cell$family & param == cell$param & n == cell$n & G == 10L & rep <= Bc]
    else fread(file.path(SIMDIR, "sim_edge_loses_pvalues.csv"))[scenario == "crossover" & n == cell$n & G == 10L & rep <= Bc]
    S <- dcast(S, rep + seed ~ test, value.var = "p_value")
  }
  setorder(S, rep); fwrite(S, f)
  S[rep <= B]
}
ID_MAP <- data.frame(harness = c("EDGE.poly2.u.G10", "EDGE.poly3.u.G10", "EDGE.stk.u.G10", "HL.G10", "HLF.G10",
                                 "HL_w", "PH", "Tsiatis", "Xie", "Stk.marg"),
                     stored = c("DEF.poly2", "DEF.poly3", "DEF.stukel", "HL", "EF", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "Stukel"),
                     required = c(rep(TRUE, 5), rep(FALSE, 5)), stringsAsFactors = FALSE)
pair_proj <- function(M, cell) {
  sk <- strsplit(cell$stored_key, "|", fixed = TRUE)[[1]]            # scenario|param of grid_proj_power.R
  S <- fread(file.path(SIMDIR, "proj_power_grid_pvalues.csv"), colClasses = list(character = "param"))[scenario == sk[1] & param == sk[2] & n == cell$n]
  m <- match(M$seed, S$seed)
  M$proj <- S$p_proj[m]                      # stored projection p-value, same data through the reused seed
  M$chk.p_edge3_stored <- S$p_edge3[m]       # identity check: must equal EDGE.poly3.u.G10
  M
}

## ---- summary ------------------------------------------------------------------------------------------------------------------
## always every finished cell of the block, whatever --cells was, so a split run never loses part of the summary (review F3)
block_summary <- function(block) {
  sel <- Cells[Cells$block == block, ]
  cache <- new.env()
  readp <- function(b, cn) {
    key <- paste(b, cn)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      f <- file.path(ROOT, b, paste0(cn, "_pvalues.csv.gz"))
      assign(key, if (file.exists(f)) fread(f) else NULL, envir = cache)
    }
    get(key, envir = cache)
  }
  mkrow <- function(ce, role, sn, Bs, t, a, rate, mcse, sap, nsz, dec, rgp, status)
    data.frame(block = block, cell = ce$cell, role = role, null_cell = ce$null_cell, n = ce$n, subset = sn, B = Bs, test = t,
               alpha = a, rejection = rate, mcse = mcse, size_adj_power = sap, null_size = nsz, declined = dec,
               rejection_given_p = rgp, status = status, stringsAsFactors = FALSE)
  rows <- list()
  for (i in seq_len(nrow(sel))) {
    ce <- sel[i, ]; P <- readp(block, ce$cell)
    if (is.null(P)) next
    if (ce$generator == "real") P <- P[P$rep == 0, ]               # the stored row order; the random orders go to _row_orders.csv
    role <- if (ce$role == "alternative" && is.na(ce$null_cell)) "alternative, no matched null declared" else ce$role
    N <- if (!is.na(ce$null_cell)) readp(ce$null_block, ce$null_cell) else NULL
    tests <- setdiff(names(P), META); tests <- tests[!grepl("^(flag\\.|lam[0-9]|chk\\.|stat\\.|info\\.)", tests)]
    ## flag rates: the Stukel one-df fallback (B), stuk_diag (Figure 8), events < G, samples with no event or no non-event,
    ## replicate errors (review F5)
    for (fl in grep("^flag\\.", names(P), value = TRUE))
      rows[[length(rows) + 1]] <- mkrow(ce, role, "all", nrow(P), fl, NA_real_, mean(P[[fl]] %in% 1), NA_real_, NA_real_, NA_real_,
                                        mean(is.na(P[[fl]])), NA_real_, "flag rate: share of replicates with the flag = 1")
    rivals_off <- "flag.rivals_run" %in% names(P) && all(P$flag.rivals_run %in% 0)
    subsets <- list(all = rep(TRUE, nrow(P)))
    for (a in grep("^flag\\.evlt\\.", names(P), value = TRUE)) {        # events < G reported separately when both occur
      fl <- P[[a]] %in% 1
      if (any(fl) && any(!fl)) { g <- sub("flag.evlt.", "", a, fixed = TRUE)
        subsets[[paste0("events<", g)]] <- fl; subsets[[paste0("events>=", g)]] <- !fl }
    }
    for (t in tests) {
      p <- P[[t]]
      pn <- if (!is.null(N) && t %in% names(N)) N[[t]] else NULL
      status <- if (rivals_off && t %in% RIVALS) "not run (n >= 10,000)"
                else if (all(!is.finite(p))) "no p-value in any replicate"
                else if (role == "alternative" && is.null(N)) "matched null not run yet" else ""
      for (sn in names(subsets)) {
        ss <- subsets[[sn]]; ps <- p[ss]; Bs <- length(ps)
        for (a in ALPHAS) {
          rate <- mean(is.finite(ps) & ps <= a)
          crit <- if (!is.null(pn)) as.numeric(stats::quantile(ifelse(is.finite(pn), pn, 1), a, type = 1)) else NA_real_
          ## C4 "both ways": unconditional (rejection, declined counted as no rejection) and among samples with a p-value
          rgp <- if (any(is.finite(ps))) sum(is.finite(ps) & ps <= a) / sum(is.finite(ps)) else NA_real_
          rows[[length(rows) + 1]] <- mkrow(ce, role, sn, Bs, t, a, rate, sqrt(rate * (1 - rate) / Bs),
            if (is.finite(crit)) mean(is.finite(ps) & ps <= crit) else NA_real_,
            if (!is.null(pn)) mean(is.finite(pn) & pn <= a) else NA_real_, mean(!is.finite(ps)), rgp, status)
        }
      }
    }
  }
  if (!length(rows)) return(invisible(NULL))
  if (any(sel$generator == "real")) row_order_summary(block, readp)
  S <- rbindlist(rows)
  fwrite(S, file.path(ROOT, block, "_summary.csv"))
  cat(sprintf("summary: %s (%d rows)\n", file.path(ROOT, block, "_summary.csv"), nrow(S)))
  invisible(S)
}

## E8.5: each real one-off fit on its stored row order (rep 0) and on the random orders (rep >= 1). Per test: the stored
## value, the median and the 5% and 95% quantiles over the random orders (type 7), the orders with a p-value, the share of
## orders with p <= 0.05, the distinct fitted risks and, for a grouped test, the group boundaries of its G that split tied
## risks in the stored order and the share of orders with such a boundary. Per-order p-values stay in the cell's file.
row_order_summary <- function(block, readp) {
  sel <- Cells[Cells$block == block & Cells$generator == "real", ]
  rows <- list()
  for (i in seq_len(nrow(sel))) {
    P <- readp(block, sel$cell[i])
    if (is.null(P) || !any(P$rep == 0)) next
    P <- as.data.frame(P); s0 <- P[P$rep == 0, , drop = FALSE][1, ]; Q <- P[P$rep > 0, , drop = FALSE]
    tests <- setdiff(names(P), META); tests <- tests[!grepl("^(flag\\.|lam[0-9]|chk\\.|stat\\.|info\\.)", tests)]
    for (t in tests) {
      arm <- regmatches(t, regexpr("(G10|Grule|G[0-9]+)$", t))
      if (!length(arm) && sel$type[i] == "external" && t %in% c(names(BT_EXT_GROUPED), "HL.ext")) arm <- "Gdef"
      tb <- if (length(arm)) paste0(if (grepl("^HLF\\.", t)) "info.tied_boundaries_hlf." else "info.tied_boundaries.", arm) else ""
      tq <- if (tb %in% names(P)) Q[[tb]] else NULL
      q <- Q[[t]]; qf <- q[is.finite(q)]
      rows[[length(rows) + 1]] <- data.frame(block = block, cell = sel$cell[i], n = s0$n, orders = nrow(Q), test = t,
        p_stored = s0[[t]], p_median = if (length(qf)) stats::median(qf) else NA_real_,
        p_q05 = if (length(qf)) as.numeric(stats::quantile(qf, 0.05)) else NA_real_,
        p_q95 = if (length(qf)) as.numeric(stats::quantile(qf, 0.95)) else NA_real_,
        orders_with_p = length(qf), share_p_le_05 = if (nrow(Q)) mean(is.finite(q) & q <= 0.05) else NA_real_,
        distinct_risks = s0$info.distinct_risks, tied_boundaries_stored = if (tb %in% names(P)) s0[[tb]] else NA_real_,
        share_orders_tied_boundary = if (!is.null(tq) && nrow(Q)) mean(is.finite(tq) & tq > 0) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) return(invisible(NULL))
  S <- rbindlist(rows)
  fwrite(S, file.path(ROOT, block, "_row_orders.csv"))
  cat(sprintf("row orders: %s (%d rows)\n", file.path(ROOT, block, "_row_orders.csv"), nrow(S)))
  invisible(S)
}

## ---- block 0 ------------------------------------------------------------------------------------------------------------------
pkg_check_rep <- function(rep, cell) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  set.seed(cell$seed_base + rep)
  dat <- bt_data(cell)
  h <- battery_rep(dat, cell)                    # no random numbers are drawn after the data
  fit <- suppressWarnings(stats::glm(dat$f, data = dat$d, family = stats::binomial()))
  G2 <- c(G10 = 10, Grule = bt_rule_G(nrow(dat$d)))
  pk <- c()
  for (a in names(G2)) for (b in BT_BASES) for (w in c("unit", "score"))
    pk[paste0("EDGE.", b, ".", if (w == "unit") "u" else "sc", ".", a)] <- tryCatch(
      ebrahim.gof::def.gof(fit, G = G2[[a]], basis = if (b == "stk") "stukel" else b, weights = w)$p_value, error = function(e) NA_real_)
  ns <- asNamespace("ebrahim.gof")
  g <- function(fun, G, o = list()) tryCatch(as.numeric(get(fun, envir = ns)(get(".gof_context", envir = ns)(fit, G = G), o)$p_value),
                                             error = function(e) NA_real_)
  pk["Stk.joint"] <- g("gof_stukel", 10, list(form = "joint"))
  pk["Stk.LR"]    <- g("gof_stukel", 10, list(form = "lr"))
  pk["Stk.marg"]  <- g("gof_stukel", 10, list(form = "marginal"))
  for (a in names(G2)) { pk[paste0("HL.", a)] <- g("gof_hl", G2[[a]]); pk[paste0("HLF.", a)] <- g("gof_ef", G2[[a]]) }
  if (nrow(dat$d) < BT_RIVAL_NMAX) {
    pk["HL_w"] <- g("gof_hlw", 10); pk["PH"] <- g("gof_ph_test", 10)
    pk["Tsiatis"] <- suppressWarnings(g("gof_tsiatis", 10)); pk["Xie"] <- suppressWarnings(g("gof_xie", 10)); pk["PR"] <- g("gof_pr", 10)
  }
  c(rep = rep, degenerate = h[["flag.degenerate"]], guard = h[["flag.info_guard"]], setNames(h[names(pk)], paste0("h.", names(pk))),
    setNames(pk, paste0("p.", names(pk))))
}
## the installed package on a sample with no event or no non-event (E8.1): def.gof gives NA with one def_degenerate warning
## in every basis and form, and gof_stukel gives NA in every form
pkg_declines <- function(fit, G = 10) {
  ns <- asNamespace("ebrahim.gof"); cls <- character(0); pd <- c()
  for (b in c("poly3", "poly2", "stukel", "sym")) for (w in c("unit", "score"))
    pd[paste(b, w)] <- tryCatch(withCallingHandlers(as.numeric(ebrahim.gof::def.gof(fit, G = G, basis = b, weights = w)$p_value),
      warning = function(x) { cls <<- c(cls, class(x)[1]); invokeRestart("muffleWarning") }), error = function(e) -1)
  ctx <- get(".gof_context", envir = ns)(fit, G = G)
  ps <- vapply(c("joint", "lr", "marginal"), function(f) tryCatch(as.numeric(suppressWarnings(
    get("gof_stukel", envir = ns)(ctx, list(form = f))$p_value)), error = function(e) -1), numeric(1))
  other <- unique(cls[cls != "def_degenerate"])
  list(ok = all(is.na(pd)) && sum(cls == "def_degenerate") == 8 && !length(other) && all(is.na(ps)),
       detail = sprintf("def.gof %d of 8 NA with %d def_degenerate warnings%s; gof_stukel %d of 3 NA", sum(is.na(pd)),
                        sum(cls == "def_degenerate"), if (length(other)) paste0(" (also ", paste(other, collapse = ", "), ")") else "",
                        sum(is.na(ps))))
}

## behaviour probe of block 0 (ii): the installed ebrahim.gof must act as the harness on three constructed data sets
## (E8.1-E8.3), or a real block 0 exits 2.
##   R2, the package's own test (Mersenne-Twister, set.seed(4): 475 risks U(0.02, 0.45) and 25 at 0.5001, G = 20, stukel):
##     the half-column of length 1.6e-7 is kept, so def.gof's unit p-value is finite and equals the harness's three-column
##     value, and the score form has df 3 and the harness's statistic (1e-10);
##   no event, n = 200: def.gof gives an NA p-value with a def_degenerate warning in every basis and form, gof_stukel NA;
##   a flat fit (two events at mirrored x, so the fitted logit is constant): the stukel and sym score forms give NA with a
##     def_no_information warning and the joint Stukel form NA, where the harness's information guard fires.
pkg_probe <- function() {
  kinds <- RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  on.exit(RNGkind(kinds[1], kinds[2], kinds[3]))
  ns <- asNamespace("ebrahim.gof")
  row <- function(cell, test, ok, diff, note) data.frame(cell = cell, test = test, required = TRUE, reps = 1L, degenerate_reps = 0L,
    max_abs_diff = diff, na_mismatch = NA_integer_, tol = "1e-10", ok = isTRUE(ok), note = note, stringsAsFactors = FALSE)
  err <- function(cell) function(e) row(cell, "error", FALSE, NA_real_, conditionMessage(e))
  r1 <- "probe: R2 construction (package test), G = 20, stukel"
  R1 <- tryCatch({
    set.seed(4)
    ph <- c(sort(stats::runif(475, 0.02, 0.45)), rep(0.5001, 25)); X <- cbind(1, stats::qlogis(ph)); y <- stats::rbinom(500, 1, ph)
    fk <- list(n = 500L, ph = ph, y = y, dmu = ph * (1 - ph), X = X, A = crossprod(X, ph * (1 - ph) * X))
    gk <- bt_groups(fk, 20L); Zk <- bt_basis(gk$pbar, "stk")
    hu <- bt_edge_unit(gk, fk$A, Zk); hs <- bt_edge_score(gk, fk$A, Zk)
    pu <- ebrahim.gof::def.gof(y, ph, X = X, G = 20, basis = "stukel")
    ps <- ebrahim.gof::def.gof(y, ph, X = X, G = 20, basis = "stukel", weights = "score")
    d <- max(abs(pu$p_value - hu$p), abs(pu$Test_Statistic - hu$S) / max(1, hu$S), abs(ps$Test_Statistic - hs$S) / max(1, hs$S))
    row(r1, "def.gof unit p finite = harness three-column value; score df 3 = harness",
        ncol(Zk) == 3L && is.finite(pu$p_value) && isTRUE(d <= 1e-10) && isTRUE(ps$df == 3) && hs$k == 3L, d,
        sprintf("package unit S %.10f p %.10f, harness S %.10f p %.10f (3 columns); score S %.6f df %s, harness df %d",
                pu$Test_Statistic, pu$p_value, hu$S, hu$p, ps$Test_Statistic, format(ps$df), hs$k))
  }, error = err(r1))
  r2 <- "probe: no event, n = 200"
  R2 <- tryCatch({
    set.seed(20260921L); x0 <- as.numeric(scale(stats::rchisq(200, 4)))
    f0 <- suppressWarnings(stats::glm(y ~ x, data = data.frame(x = x0, y = rep(0, 200)), family = stats::binomial()))
    dg <- pkg_declines(f0)
    row(r2, "def.gof NA with def_degenerate (4 bases x 2 forms); gof_stukel NA (3 forms)", dg$ok, NA_real_, dg$detail)
  }, error = err(r2))
  r3 <- "probe: flat fit (two events at mirrored x), n = 200, G = 10"
  R3 <- tryCatch({
    xg <- stats::qnorm(stats::ppoints(200)); yg <- rep(0, 200); yg[c(1, 200)] <- 1
    dg <- data.frame(x = xg, y = yg)
    fg <- suppressWarnings(stats::glm(y ~ x, data = dg, family = stats::binomial()))
    h <- battery_rep(list(d = dg, f = y ~ x), list(generator = "sparse", n = 200L, type = "full", ao = FALSE, G_extra = ""))
    cls <- character(0)
    pk <- vapply(c("stukel", "sym"), function(b) withCallingHandlers(as.numeric(ebrahim.gof::def.gof(fg, G = 10, basis = b, weights = "score")$p_value),
      warning = function(w) { cls <<- c(cls, class(w)[1]); invokeRestart("muffleWarning") }), numeric(1))
    pj <- as.numeric(suppressWarnings(get("gof_stukel", envir = ns)(get(".gof_context", envir = ns)(fg, G = 10), list(form = "joint"))$p_value))
    hh <- h[c("EDGE.stk.sc.G10", "EDGE.sym.sc.G10", "Stk.joint")]
    row(r3, "stukel and sym score NA with def_no_information; joint Stukel NA; harness NA with flag.info_guard",
        all(is.na(pk)) && sum(cls == "def_no_information") == 2 && is.na(pj) && all(!is.finite(hh)) && h[["flag.info_guard"]] == 1, NA_real_,
        sprintf("package score %s, joint %s, warnings %s; harness %s, flag.info_guard %d", paste(format(pk), collapse = "/"), format(pj),
                paste(unique(cls), collapse = ", "), paste(format(hh), collapse = "/"), as.integer(h[["flag.info_guard"]])))
  }, error = err(r3))
  rows <- rbind(R1, R2, R3)
  list(ok = all(rows$ok), rows = rows,
       detail = paste(sprintf("%s %s", c("R2", "no event", "flat fit"), ifelse(rows$ok, "ok", "FAIL")), collapse = ", "))
}
PKG_REQUIRED <- c(paste0("EDGE.", rep(BT_BASES, each = 4), ".", c("u", "sc"), ".", rep(c("G10", "Grule"), each = 2)),
                  "Stk.joint", "HL.G10", "HL.Grule", "HLF.G10", "HLF.Grule", RIVALS)

run_block0 <- function(B0) {
  dir.create(file.path(ROOT, "0"), recursive = TRUE, showWarnings = FALSE)
  write_rule()
  sel <- Cells[Cells$block == "0", ]
  if (TEST_MODE) sel$B <- B0 else sel$B <- pmax(sel$B, B0)
  B0 <- sel$B[1]
  plog("0", sprintf("block 0: B = %d, workers %d, ebrahim.gof installed %s", B0, W, pkg_version()))

  ## (i) old-seed identity
  cl <- bt_cluster(W)
  ID <- list()
  for (i in which(!is.na(sel$stored))) {
    ce <- as.list(sel[i, ])
    t1 <- Sys.time(); M <- run_cell(cl, ce, B0)
    write_cell(M, file.path(ROOT, "0", paste0(ce$cell, "_pvalues.csv.gz")))
    S <- stored_pvalues(ce, B0)
    seeds_ok <- nrow(S) == nrow(M) && all(S$seed == M$seed)
    for (k in seq_len(nrow(ID_MAP))) {
      h <- M[[ID_MAP$harness[k]]]; s <- S[[ID_MAP$stored[k]]]
      both <- is.finite(h) & is.finite(s)
      md <- if (any(both)) max(abs(h[both] - s[both])) else NA_real_
      mism <- sum(is.finite(h) != is.finite(s))
      ID[[length(ID) + 1]] <- data.frame(cell = ce$cell, harness = ID_MAP$harness[k], stored = ID_MAP$stored[k],
        required = ID_MAP$required[k], reps = nrow(M), seeds_match = seeds_ok, max_abs_diff = md, na_mismatch = mism,
        ok = seeds_ok && mism == 0 && isTRUE(md <= 1e-6), stringsAsFactors = FALSE)
    }
    plog("0", sprintf("(i) %-24s %d replicates in %.1f s", ce$cell, B0, as.numeric(difftime(Sys.time(), t1, units = "secs"))))
  }
  ID <- do.call(rbind, ID)
  fwrite(ID, file.path(ROOT, "0", "identity.csv"))
  pass_i <- all(ID$ok[ID$required])
  print(ID, row.names = FALSE, digits = 3)
  plog("0", sprintf("(i) identity with the stored July p-values (required tests): %s", if (pass_i) "PASS" else "FAIL"))

  ## (ii) harness against the package
  pass_ii <- NA
  if (!pkg_ok()) {
    plog("0", sprintf("(ii) PENDING: ebrahim.gof %s installed, %s or later needed; the package was not loaded", pkg_version(), PKG_MIN))
  } else {
    stopCluster(cl); cl <- bt_cluster(W, load_pkg = TRUE)
    clusterExport(cl, "pkg_check_rep")
    pr <- pkg_probe()
    plog("0", sprintf("(ii) behaviour probe of the installed ebrahim.gof: %s (%s)", if (pr$ok) "PASS" else "FAIL", pr$detail))
    PK <- list()
    for (i in which(!grepl("^w_", sel$cell))) {
      ce <- as.list(sel[i, ])
      chunks <- splitIndices(B0, min(B0, W * 4L))
      M <- as.data.frame(do.call(rbind, parLapplyLB(cl, chunks, function(r, cell) do.call(rbind, lapply(r, pkg_check_rep, cell = cell)), cell = ce)))
      ## a sample with no event or no non-event has no fitted model: the harness returns no p-value for any test there
      ## (E8.1), the package only for its def.gof and Stukel rows (its HL, HL_F, Tsiatis and Xie rows still return one).
      ## Such replicates are the only ones left out of the comparison, and must be NA in the harness; replicates where the
      ## information guard fired are compared like any other (E8.2).
      use <- M$degenerate %in% 0
      grd <- M$guard %in% 1                   # a score statistic left out a column with information below 1e-10
      for (t in sub("^p\\.", "", grep("^p\\.", names(M), value = TRUE))) {
        guarded <- grepl("\\.sc\\.", t) || t == "Stk.joint"
        keep <- use
        h <- M[[paste0("h.", t)]][keep]; p <- M[[paste0("p.", t)]][keep]
        both <- is.finite(h) & is.finite(p)
        md <- if (any(both)) max(abs(h[both] - p[both])) else NA_real_
        tol <- if (t %in% c("Stk.LR", "Stk.marg")) 1e-6 else 1e-8
        hdeg <- all(!is.finite(M[[paste0("h.", t)]][!use]))
        PK[[length(PK) + 1]] <- data.frame(cell = ce$cell, test = t, required = t %in% PKG_REQUIRED, reps = sum(keep),
          degenerate_reps = sum(!use), max_abs_diff = md, na_mismatch = sum(is.finite(h) != is.finite(p)), tol = format(tol),
          ok = sum(is.finite(h) != is.finite(p)) == 0 && (is.na(md) || md <= tol) && hdeg,
          note = paste(c(if (sum(!use)) sprintf("harness NA in the %d degenerate replicates: %s", sum(!use), hdeg),
                         if (guarded && sum(use & grd)) sprintf("%d replicates with the information guard, compared", sum(use & grd))),
                       collapse = "; "),
          stringsAsFactors = FALSE)
      }
      plog("0", sprintf("(ii) %-24s checked against the package (%d degenerate replicates left out)", ce$cell, sum(!use)))
    }

    ## constructed samples with no event (y = 0) and no non-event (y = 1), n = 200 (E8.1, E8.6): every harness test is NA, and
    ## the package gives no p-value in def.gof (4 bases x 2 forms, each with a def_degenerate warning) and gof_stukel (joint, lr,
    ## marginal). Its HL, HL_F, Tsiatis and Xie rows, which still return a p-value there, are listed for information.
    zc <- as.list(Cells[Cells$block == "0" & Cells$cell == "v_sparse49_n300", ]); zc$n <- 200L
    nsp <- asNamespace("ebrahim.gof")
    pv <- function(expr) tryCatch(format(suppressWarnings(expr), digits = 3), error = function(e) "error")
    for (k in 1:6) {
      set.seed(20260914L + k); x <- as.numeric(scale(rchisq(200, 4))); yv <- if (k <= 3) 0 else 1
      dat <- list(d = data.frame(x = x, y = rep(yv, 200)), f = y ~ x)
      h <- battery_rep(dat, zc)
      tc <- names(h)[!grepl("^(flag\\.|n$|events$|fit_ok$|glm_conv$)", names(h))]
      fit <- suppressWarnings(stats::glm(y ~ x, data = dat$d, family = stats::binomial()))
      dg <- pkg_declines(fit)
      ctx <- get(".gof_context", envir = nsp)(fit, G = 10)
      PK[[length(PK) + 1]] <- data.frame(cell = sprintf("constructed, all y = %d, n = 200, set %d", yv, k),
        test = "every harness test NA; package def.gof (8) and gof_stukel (3) NA", required = TRUE, reps = 1L, degenerate_reps = 1L,
        max_abs_diff = NA_real_, na_mismatch = NA_integer_, tol = "",
        ok = all(!is.finite(h[tc])) && h[["flag.degenerate"]] == 1 && dg$ok,
        note = sprintf("package: %s; for information HL %s, HL_F %s, Tsiatis %s, Xie %s", dg$detail,
          pv(get("gof_hl", envir = nsp)(ctx)$p_value), pv(get("gof_ef", envir = nsp)(ctx)$p_value),
          pv(get("gof_tsiatis", envir = nsp)(ctx)$p_value), pv(get("gof_xie", envir = nsp)(ctx)$p_value)),
        stringsAsFactors = FALSE)
    }

    ## the column rule (E8.3) on a Stukel half-column that is tiny but not zero (one group at mean risk 0.5001), G = 20
    set.seed(20260915L); ph <- c(sort(runif(475, 0.02, 0.45)), rep(0.5001, 25)); yk <- rbinom(500, 1, ph); Xk <- cbind(1, stats::qlogis(ph))
    fk <- list(n = 500L, ph = ph, y = yk, dmu = ph * (1 - ph), X = Xk, A = crossprod(Xk, ph * (1 - ph) * Xk))
    gk <- bt_groups(fk, 20L); Zk <- bt_basis(gk$pbar, "stk")
    for (w in c("unit", "score")) {
      hp <- if (w == "unit") bt_edge_unit(gk, fk$A, Zk)$p else bt_edge_score(gk, fk$A, Zk)$p
      pp <- tryCatch(suppressWarnings(ebrahim.gof::def.gof(yk, predicted_probs = ph, X = Xk, G = 20, basis = "stukel", weights = w)$p_value),
                     error = function(e) NA_real_)
      PK[[length(PK) + 1]] <- data.frame(cell = "constructed, Stukel half-column 1.6e-7 long, G = 20", test = paste0("EDGE.stk.", w),
        required = TRUE, reps = 1L, degenerate_reps = 0L, max_abs_diff = abs(hp - pp), na_mismatch = as.integer(is.finite(hp) != is.finite(pp)),
        tol = "1e-8", ok = isTRUE(abs(hp - pp) <= 1e-8),
        note = "the tiny half-column is kept and every column scaled to unit length, in harness and package (E8.3)",
        stringsAsFactors = FALSE)
    }
    PK <- do.call(rbind, c(PK, list(pr$rows)))
    fwrite(PK, file.path(ROOT, "0", "package_check.csv"))
    pass_ii <- all(PK$ok[PK$required & !grepl("^probe: ", PK$cell)])
    if (!pr$ok) {                           # the harness follows the corrected 2.8.0 build, so the package must behave as it does
      if (TEST_MODE) plog("0", "(ii) test run: the installed package fails the behaviour probe; accepted for a test run only")
      else { plog("0", "(ii) PENDING: the installed ebrahim.gof fails the behaviour probe (R2 construction, no-event sample, flat fit); install the corrected 2.8.0 build, then run block 0 again"); pass_ii <- NA }
    }
    print(PK[!PK$ok | PK$required & PK$test %in% c("Stk.joint", "EDGE.sym.sc.Grule") | grepl("^probe: ", PK$cell), ], row.names = FALSE, digits = 3)
    plog("0", sprintf("(ii) harness against ebrahim.gof %s (required tests): %s", pkg_version(), if (pass_ii) "PASS" else "FAIL"))
  }
  stopCluster(cl)

  ## (iii) workers (E0.8)
  grid <- as.integer(strsplit(opt[["worker-grid"]], ",")[[1]])
  if (TEST_MODE) { grid <- grid[grid <= 8L]; if (!length(grid)) grid <- c(4L, 8L) }
  else if (!all(c(12L, 16L, 20L) %in% grid)) stop("E0.8: the real block 0 times 12, 16 and 20 workers (--worker-grid ", opt[["worker-grid"]], ")")
  wc <- as.list(sel[grepl("^w_", sel$cell), ][1, ])
  pin_blas()
  invisible(battery_one(4L, wc))                                      # warm-up: namespaces and first-call costs (review f1)
  ts <- system.time(for (r in 1:3) battery_one(r, wc))[["elapsed"]] / 3
  WT <- list()
  for (w in grid) {
    clw <- bt_cluster(w)
    invisible(clusterCall(clw, function(ce) { battery_one(1L, ce); NULL }, wc))     # warm-up, not timed
    t1 <- Sys.time(); run_cell(clw, wc, B0); el <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    stopCluster(clw)
    WT[[length(WT) + 1]] <- data.frame(workers = w, B = B0, n = wc$n, wall_s = el, reps_per_s = B0 / el,
                                       serial_s_per_rep = ts, speedup = B0 * ts / el, efficiency = B0 * ts / el / w)
    plog("0", sprintf("(iii) %2d workers: %d replicates at n = %d in %.1f s (speed-up %.1f)", w, B0, wc$n, el, B0 * ts / el))
  }
  WT <- do.call(rbind, WT)
  fwrite(WT, file.path(ROOT, "0", "workers.csv"))
  best <- WT$workers[which.max(WT$reps_per_s)]
  writeLines(as.character(best), file.path(ROOT, "0", "workers_choice.txt"))
  plog("0", sprintf("(iii) fastest of %s workers: %d (written to %s)", paste(grid, collapse = ", "), best, file.path(ROOT, "0", "workers_choice.txt")))

  ok <- pass_i && isTRUE(pass_ii)
  mk <- file.path(ROOT, "0", "_PASSED")
  if (ok) writeLines(c(sprintf("B=%d", B0), format(Sys.time()), paste("ebrahim.gof", pkg_version())), mk)
  else if (file.exists(mk)) file.remove(mk)
  plog("0", sprintf("block 0: %s", if (ok) "PASSED" else if (pass_i && is.na(pass_ii)) "NOT PASSED (package check pending)" else "FAILED"))
  quit(save = "no", status = if (ok) 0L else if (pass_i && is.na(pass_ii)) 2L else 1L)
}

## ---- the timing sample ---------------------------------------------------------------------------------------------------------
run_estimate <- function(blocks) {
  pin_blas()
  z <- Cells[Cells$block %in% blocks, ]
  z$key <- paste(z$block, z$type, z$generator, z$family, z$xdist, z$ao, z$n, z$G_extra, z$dataset, z$pstar)
  reps <- z[!duplicated(z$key), ]
  cat(sprintf("timing %d representative cells serially (BLAS on one thread)\n", nrow(reps)))
  reps$sec_per_rep <- NA_real_
  for (i in seq_len(nrow(reps))) {
    ce <- as.list(reps[i, ])
    k <- if (ce$generator == "real") 1L else if (ce$n <= 5000) 3L else if (ce$n <= 20000) 2L else 1L
    invisible(battery_one(k + 1L, ce))                                   # first call loads namespaces
    reps$sec_per_rep[i] <- system.time(for (r in seq_len(k)) battery_one(r, ce))[["elapsed"]] / k
  }
  z$sec_per_rep <- reps$sec_per_rep[match(z$key, reps$key)]
  z$cpu_s <- z$B * z$sec_per_rep
  fwrite(z[, c("block", "cell", "n", "B", "type", "sec_per_rep", "cpu_s")], edge_battery("_estimate_cells.csv"))
  wf <- file.path(ROOT, "0", "workers.csv")                             # --root battery/_test2 reads a test block 0's timing
  eff <- if (file.exists(wf)) { wt <- fread(wf); wt$efficiency[which.min(abs(wt$workers - W))] } else 7.5 / 20
  src <- if (file.exists(wf)) "block 0 worker timing" else "MAP_battery.md 4.1 (7.5x at 20 workers)"
  E <- as.data.table(z)[, list(cells = .N, replicates = sum(B), cpu_h = sum(cpu_s) / 3600), by = block]
  E$wall_h_ideal <- E$cpu_h / W
  E$wall_h_observed <- E$cpu_h / (W * eff)
  E <- rbind(E, data.frame(block = "total", cells = sum(E$cells), replicates = sum(E$replicates), cpu_h = sum(E$cpu_h),
                           wall_h_ideal = sum(E$wall_h_ideal), wall_h_observed = sum(E$wall_h_observed)))
  fwrite(E, edge_battery("_estimate_blocks.csv"))
  cat(sprintf("\nprojection at %d workers; observed efficiency %.2f from %s\n", W, eff, src))
  print(E, digits = 3)
}

## ---- main ---------------------------------------------------------------------------------------------------------------------
Cells <- battery_cells()
if (isTRUE(opt$list)) {
  fwrite(Cells, edge_battery("cells.csv"))
  print(table(Cells$block)); cat("replicates:", sum(Cells$B), "\n")
  quit(save = "no", status = 0L)
}
if (isTRUE(opt$estimate)) {
  blocks <- if (is.null(opt$block) || opt$block == "all") setdiff(BC_BLOCKS, "0") else strsplit(opt$block, ",")[[1]]
  run_estimate(blocks); quit(save = "no", status = 0L)
}
if (is.null(opt$block) || !opt$block %in% BC_BLOCKS) stop("--block must be one of ", paste(BC_BLOCKS, collapse = ", "))
block <- opt$block
dir.create(file.path(ROOT, block), recursive = TRUE, showWarnings = FALSE)
sel <- Cells[Cells$block == block, ]
if (!is.null(opt$cells)) {
  want <- strsplit(opt$cells, ",")[[1]]
  if (any(!want %in% sel$cell)) stop("unknown cells in block ", block, ": ", paste(setdiff(want, sel$cell), collapse = ", "))
  sel <- sel[sel$cell %in% want, ]
}
if (isTRUE(opt$summary)) { block_summary(block); quit(save = "no", status = 0L) }
if (block == "0") run_block0(if (TEST_MODE) as.integer(opt$B) else 200L)

p0 <- file.path(ROOT, "0", "_PASSED")
if (!file.exists(p0)) {
  if (TEST_MODE) warning("block 0 has not passed in ", ROOT, "; continuing because this is a test run (--B)")
  else stop("refusing to start block ", block, ": block 0 has not passed (", p0, ")")
} else if (!TEST_MODE) {
  b0 <- as.integer(sub("B=", "", readLines(p0, warn = FALSE)[1]))
  if (b0 < 200L) stop("refusing to start block ", block, ": block 0 passed only at B = ", b0)
}
if (block == "2") check_rule()

sel$B_run <- if (TEST_MODE) ifelse(sel$generator == "real", pmin(sel$B, 1L + as.integer(opt$B)), as.integer(opt$B)) else sel$B   # real: stored order + B orders
sel$out <- file.path(ROOT, block, paste0(sel$cell, "_pvalues.csv.gz"))
done <- file.exists(sel$out)
wcost <- sel$B_run * cost_guess(sel)
plog(block, sprintf("block %s: %d cells (%d done), %d replicates to run, %d workers, root %s%s", block, nrow(sel), sum(done),
                    sum(sel$B_run[!done]), W, ROOT, if (TEST_MODE) " [test run]" else ""))
for (i in which(done)) {                                            # a finished cell is skipped only if it has B rows (review F11)
  nr <- nrow(fread(sel$out[i], select = "rep"))
  if (nr != sel$B_run[i]) stop(sprintf("%s exists with %d rows but B is %d: move the file away and start again", sel$out[i], nr, sel$B_run[i]))
}
if (any(!done)) {
  cl <- bt_cluster(W)
  t0 <- Sys.time(); wtot <- sum(wcost[!done]); wdone <- 0; k <- 0L
  for (i in which(!done)) {
    k <- k + 1L; ce <- as.list(sel[i, ]); t1 <- Sys.time()
    M <- run_cell(cl, ce, sel$B_run[i])
    extra <- ""
    if (isTRUE(ce$pair_proj)) {
      M <- pair_proj(M, ce)
      both <- is.finite(M$chk.p_edge3_stored) & is.finite(M$EDGE.poly3.u.G10)
      extra <- sprintf("  stored p_edge3 max|diff| %.1e", if (any(both)) max(abs(M$chk.p_edge3_stored[both] - M$EDGE.poly3.u.G10[both])) else NA)
    }
    write_cell(M, sel$out[i])
    wdone <- wdone + wcost[i]
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    eta <- el / wdone * (wtot - wdone)
    plog(block, sprintf("[%d/%d] %-30s B=%-5d n=%-6s %7.1f s  rep errors %d  %4.1f%% of cost  ETA %s (%.2f h)%s",
                        k, sum(!done), ce$cell, sel$B_run[i], ce$n, as.numeric(difftime(Sys.time(), t1, units = "secs")),
                        sum(M$flag.rep_error %in% 1), 100 * wdone / wtot, format(Sys.time() + eta, "%Y-%m-%d %H:%M"), eta / 3600, extra))
  }
  stopCluster(cl)
}
block_summary(block)
