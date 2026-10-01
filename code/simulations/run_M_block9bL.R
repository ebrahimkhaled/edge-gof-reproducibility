## run_M_block9bL.R -- block 9bL: le Cessie on the 300 data sets of block 9b.
## Contract: paper_EDGE/theory/PREDECLARATION_block9bL_lecessie.md (sha256 e3c89be7...), frozen before any data set
## was tested. No claim: every rate is reported beside EDGE-poly3 (unit, rule G) on the same replicates.
##
##   Rscript run_M_block9bL.R
##
## Serial by design (section 3.3). Each data set is regenerated from the seed block 9b stored (b9b_gen after
## set.seed), kept only if EDGE-poly3 unit at G = 10 recomputed on it matches the stored value to 1e-8 and seed,
## n and events agree, and then tested with block 8L's le Cessie call (l8_lecessie), which is checked against its
## definition before anything runs.
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
options(block8L.source_only = TRUE)
suppressMessages(source(file.path(SIMDIR, "run_M_block8L.R")))      # l8_lecessie, l8_transpose_check, rv_identity
suppressMessages(source(file.path(SIMDIR, "_block9b_rivals.R")))    # b9b_cells, b9b_gen
OUT9 <- edge_battery("9bL")
dir.create(OUT9, showWarnings = FALSE, recursive = TRUE)
lg <- function(msg) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "", msg)
  cat(line, "\n", sep = ""); cat(line, "\n", sep = "", file = file.path(OUT9, "_progress.log"), append = TRUE)
}

## the implementation is the declared one (section 1)
if (utils::packageVersion("ebrahim.gof") < L8_MINVER) stop("ebrahim.gof >= ", L8_MINVER, " is required")
tc <- l8_transpose_check()
if (!isTRUE(abs(tc[["package"]] - tc[["declared"]]) < 1e-10) || !isTRUE(abs(tc[["package"]] - tc[["vendor"]]) > 1e-6))
  stop("the installed le Cessie is not the corrected (I-H)'R(I-H)")
lg(sprintf("block 9bL: le Cessie check passed (package %.12f = declared %.12f; vendor %.12f)", tc[["package"]], tc[["declared"]], tc[["vendor"]]))

C <- b9b_cells()
RNGkind("L'Ecuyer-CMRG")
for (i in seq_len(nrow(C))) {
  cell <- as.list(C[i, ])
  out <- file.path(OUT9, sprintf("%s_lecessie_pvalues.csv.gz", cell$cell))
  if (file.exists(out)) { lg(sprintf("%s: finished file present, skipped", cell$cell)); next }
  S <- data.table::fread(edge_battery("9b", sprintf("%s_pvalues.csv.gz", cell$cell)),
                         select = c("rep", "seed", "n", "events", "EDGE.poly3.u.G10"))
  if (nrow(S) != 100L || any(S$rep != 1:100)) stop(cell$cell, ": block 9b's file does not hold replicates 1-100")
  rows <- lapply(seq_len(nrow(S)), function(j) {
    r <- S$rep[j]; seed <- cell$seed_base + r
    set.seed(seed)
    g <- b9b_gen(cell$mult)
    dat <- list(d = g$d, f = g$f); n <- nrow(g$d); ev <- sum(g$d$y)
    fq <- bt_fit(dat)
    p10 <- if (is.null(fq)) NA_real_ else tryCatch(bt_edge_arm(fq, 10L, cell)[["EDGE.poly3.u"]], error = function(e) NA_real_)
    st <- list(p = S$EDGE.poly3.u.G10[j], seed = S$seed[j], n = S$n[j], events = S$events[j])
    ok <- rv_identity(list(p = p10, n = n, events = ev), st, seed)
    res <- c(rep = r, seed = seed, n = n, events = ev, id.absdiff = abs(p10 - st$p), id_ok = as.numeric(ok),
             lecessie = NA_real_, lc.stat = NA_real_, lc.df = NA_real_, sec = NA_real_)
    if (ok && !is.null(fq)) {
      t0 <- proc.time()[["elapsed"]]
      lc <- l8_lecessie(fq$fit)
      res[c("lecessie", "lc.stat", "lc.df")] <- lc[c("p", "stat", "df")]
      res["sec"] <- proc.time()[["elapsed"]] - t0
    }
    res
  })
  M <- as.data.frame(do.call(rbind, rows))
  data.table::fwrite(M, out)
  lg(sprintf("%-6s kept %d of %d, identity failures %d, no p-value %d, le Cessie rejection %.3f, median %.2f s",
             cell$cell, sum(M$id_ok == 1), nrow(M), sum(M$id_ok == 0), sum(M$id_ok == 1 & !is.finite(M$lecessie)),
             mean(M$id_ok == 1 & is.finite(M$lecessie) & M$lecessie <= 0.05), stats::median(M$sec, na.rm = TRUE)))
}

## section 3.2: the rates beside EDGE-poly3 (unit, rule G) on the same replicates, paired, Holm across the three cells
Sm <- data.table::rbindlist(lapply(C$cell, function(ce) {
  L <- data.table::fread(file.path(OUT9, sprintf("%s_lecessie_pvalues.csv.gz", ce)))
  B <- data.table::fread(edge_battery("9b", sprintf("%s_pvalues.csv.gz", ce)), select = c("rep", "EDGE.poly3.u.Grule"))
  K <- L[id_ok == 1]
  pe <- B$EDGE.poly3.u.Grule[match(K$rep, B$rep)]
  r1 <- is.finite(K$lecessie) & K$lecessie <= 0.05; r2 <- is.finite(pe) & pe <= 0.05
  nb <- sum(r1 & !r2); nc <- sum(!r1 & r2)
  data.table::data.table(cell = ce, reps_kept = nrow(K), identity_failures = sum(L$id_ok == 0),
                         max_identity_absdiff = max(K$id.absdiff, na.rm = TRUE), no_pvalue = sum(!is.finite(K$lecessie)),
                         lecessie_rejection = mean(r1), edge_rejection = mean(r2), lecessie_only = nb, edge_only = nc,
                         mcnemar_p = if (nb + nc == 0) 1 else stats::binom.test(nb, nb + nc, 0.5)$p.value,
                         median_sec = stats::median(K$sec, na.rm = TRUE))
}))
Sm[, holm_p := stats::p.adjust(mcnemar_p, method = "holm")]
data.table::fwrite(Sm, file.path(OUT9, "_summary.csv"))
print(Sm)
lg("block 9bL: run finished")
