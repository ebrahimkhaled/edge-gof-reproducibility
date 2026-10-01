## run_M_block9.R -- block 9 of the EDGE restructure: a few corrupted records.
## Contract: paper_EDGE/theory/PREDECLARATION_block9_contamination.md (sha256 99eb4c27...), sections 1 and 2.
##
##   Rscript run_M_block9.R [--cells a,b,...] [--workers 20] [--reps N] [--out battery/9] [--force]
##
## 84 cells, B = 1000 replicates each, seeds set.seed(300000000 + cell_id * 10000 + rep). Each replicate draws the
## base design (x ~ U(-3,3), d ~ Bernoulli(0.5), eta = 0.6x + 0.5d), applies the cell's corruption, fits y ~ x + d
## and runs the battery's own test code on it, so block 9 and blocks 0-8 compute the same statistics.
##
## Nothing in blocks 0-8 is read or written. A run refuses to write anywhere inside those folders, and it takes a
## lock so two runs cannot interleave on the same output. A cell's replicates accumulate in <cell>_pvalues.csv.gz.part,
## written whole and swapped in, so an interrupted run resumes without recomputing and never leaves a truncated file.
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
suppressMessages({ library(parallel); library(data.table) })

SIMDIR <- edge_path("code/simulations")
ROOT   <- edge_battery()
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_block9_contam.R"))

B9_BATCH <- 4L                                        # replicates per worker per batch (one .part write per batch)

## ---- options ---------------------------------------------------------------------------------------------------
b9_opts <- function(a = commandArgs(TRUE)) {
  opt <- list()
  i <- 1L
  while (i <= length(a)) {
    if (!startsWith(a[i], "--")) stop("unexpected argument ", a[i])
    key <- sub("^--", "", a[i])
    if (key == "force") { opt[[key]] <- TRUE; i <- i + 1L }
    else { if (i == length(a)) stop("--", key, " needs a value"); opt[[key]] <- a[i + 1L]; i <- i + 2L }
  }
  opt
}
OPT <- b9_opts()

b9_abs <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
OUT <- b9_abs(if (is.null(OPT$out)) edge_battery("9") else
              if (grepl("^([A-Za-z]:|/)", OPT$out)) OPT$out else file.path(SIMDIR, OPT$out))

## the folders this run must never write into: every finished block and the dry run
busy <- b9_abs(edge_battery(c(as.character(0:8), "dryrun", "analysis")))
if (any(startsWith(paste0(OUT, "/"), paste0(busy, "/")))) stop("--out cannot be a folder of blocks 0-8: ", OUT)

W <- if (is.null(OPT$workers)) 20L else as.integer(OPT$workers)
if (is.na(W) || W < 1L) stop("--workers needs a positive integer")
TEST_RUN <- !is.null(OPT$reps)
REPS_OVERRIDE <- if (TEST_RUN) as.integer(OPT$reps) else NA_integer_
if (TEST_RUN && (is.na(REPS_OVERRIDE) || REPS_OVERRIDE < 1L)) stop("--reps needs a positive integer")
if (TEST_RUN && identical(OUT, b9_abs(edge_battery("9"))))
  stop("a test run (--reps) cannot write to the real block 9 folder; pass --out")

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "_progress.log")
b9_log <- function(msg) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "", msg)
  cat(line, "\n", sep = ""); cat(line, "\n", sep = "", file = LOG, append = TRUE)
}

## ---- cells ----------------------------------------------------------------------------------------------------
CELLS <- b9_cells()
if (nrow(CELLS) != 84L) stop("block 9 must have 84 cells, not ", nrow(CELLS))
if (!is.null(OPT$cells)) {
  want <- trimws(strsplit(OPT$cells, ",")[[1]])
  if (any(!want %in% CELLS$cell)) stop("not a block 9 cell: ", paste(setdiff(want, CELLS$cell), collapse = ", "))
  CELLS <- CELLS[CELLS$cell %in% want, , drop = FALSE]
}
b9_out_path <- function(ce) file.path(OUT, paste0(ce$cell, "_pvalues.csv.gz"))
b9_reps <- function(ce) if (TEST_RUN) REPS_OVERRIDE else as.integer(ce$B)

## ---- refusals, all found before any computation ----------------------------------------------------------------
b9_check <- function(C) {
  bad <- character(0)
  for (i in seq_len(nrow(C))) {
    ce <- as.list(C[i, ]); R <- b9_reps(ce)
    o <- b9_out_path(ce); pt <- paste0(o, ".part")
    g <- tryCatch(b9_gen(ce$truth, ce$n, ce$corruption, ce$rate), error = function(e) NULL)
    if (is.null(g)) { bad <- c(bad, sprintf("%s: the generator fails", ce$cell)); next }
    if (nrow(g$d) != ce$n || !setequal(names(g$d), c("y", "x", "d")))
      bad <- c(bad, sprintf("%s: the generator returns the wrong data frame", ce$cell))
    if (length(g$corrupt) != ce$k)
      bad <- c(bad, sprintf("%s: the generator corrupted %d records, the cell table says %d",
                            ce$cell, length(g$corrupt), ce$k))
    if (file.exists(o)) {
      nr <- nrow(fread(o, select = "rep"))
      if (nr != R) bad <- c(bad, sprintf("%s exists with %d rows but %d replicates are planned: move it away", o, nr, R))
    } else if (file.exists(pt) || file.exists(paste0(pt, ".tmp"))) {
      P <- b9_read_part(pt)
      if (is.null(P)) bad <- c(bad, sprintf("%s (and its .tmp) cannot be read: move them away", pt))
      else if (anyDuplicated(P$rep) || any(!P$rep %in% seq_len(R)))
        bad <- c(bad, sprintf("%s does not belong to this run: move it away", pt))
    }
  }
  bad
}

b9_read_part <- function(part) {
  for (f in c(part, paste0(part, ".tmp"))) {
    if (!file.exists(f)) next
    P <- tryCatch(as.data.frame(fread(f)), error = function(e) NULL)
    if (!is.null(P)) return(P)
  }
  NULL
}

b9_pid_alive <- function(pid) {
  if (length(pid) != 1L || is.na(pid)) return(FALSE)
  if (.Platform$OS.type == "windows") {
    tl <- tryCatch(system2("tasklist", c("/FI", shQuote(paste("PID eq", pid), type = "cmd"), "/NH"),
                           stdout = TRUE, stderr = FALSE), error = function(e) character(0))
    return(any(grepl(paste0(" ", pid, " "), paste0(" ", tl, " "), fixed = TRUE)))
  }
  isTRUE(tools::pskill(pid, 0L))
}

## ---- progress --------------------------------------------------------------------------------------------------
b9_cost <- function(n) n / 1000                        # a replicate costs about linearly in n
b9_prog <- function(total) {
  t0 <- Sys.time(); done <- 0
  list(add = function(x) done <<- done + x,
       eta = function() {
         if (done <= 0) return("-")
         el <- as.numeric(difftime(Sys.time(), t0, units = "secs")); rem <- el / done * max(total - done, 0)
         sprintf("%s (%.2f h)", format(Sys.time() + rem, "%Y-%m-%d %H:%M"), rem / 3600)
       })
}

## ---- workers ---------------------------------------------------------------------------------------------------
b9_cluster <- function(w) {
  cl <- makePSOCKcluster(w)
  parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
  clusterEvalQ(cl, RNGkind("L'Ecuyer-CMRG"))
  clusterExport(cl, "SIMDIR", envir = environment())
  clusterEvalQ(cl, {
    source(file.path(SIMDIR, "_battery_tests.R"))
    source(file.path(SIMDIR, "_block9_contam.R"))
    TRUE
  })
  cl
}
b9_task <- function(rep, cell) b9_one(rep, cell)

## ---- one cell --------------------------------------------------------------------------------------------------
b9_run_cell <- function(cl, ce, R, prog) {
  out <- b9_out_path(ce); part <- paste0(out, ".part")
  M <- b9_read_part(part)
  todo <- setdiff(seq_len(R), M$rep)
  b9_log(sprintf("%s: %d of %d replicates to run%s", ce$cell, length(todo), R,
                 if (is.null(M)) "" else sprintf(" (%d resumed from .part)", nrow(M))))
  bs <- length(cl) * B9_BATCH
  for (b in split(todo, ceiling(seq_along(todo) / bs))) {
    rows <- clusterApplyLB(cl, b, b9_task, cell = ce)
    M <- rbind(M, as.data.frame(do.call(rbind, rows)))
    M <- M[order(M$rep), , drop = FALSE]; rownames(M) <- NULL
    tmp <- paste0(part, ".tmp")                        # write whole, then swap in: a crash never truncates .part
    fwrite(M, tmp, compress = "gzip")
    if (file.exists(part)) file.remove(part)
    if (!file.rename(tmp, part)) stop("could not rename ", tmp)
    prog$add(length(b) * b9_cost(ce$n))
    b9_log(sprintf("  %-26s %d/%d replicates written to .part  ETA %s", ce$cell, nrow(M), R, prog$eta()))
  }
  if (is.null(M) || nrow(M) != R || any(M$rep != seq_len(R))) stop("internal: ", part, " is incomplete")
  if (file.exists(out)) file.remove(out)
  if (!file.rename(part, out)) stop("could not rename ", part)
  if (file.exists(paste0(part, ".tmp"))) file.remove(paste0(part, ".tmp"))
  M
}

## ---- main ------------------------------------------------------------------------------------------------------
b9_log(sprintf("block 9: %d cells, B = %s, %d workers, out %s", nrow(CELLS),
               if (TEST_RUN) paste0(REPS_OVERRIDE, " (TEST RUN)") else "1000", W, OUT))

bad <- b9_check(CELLS)
if (length(bad)) { for (b in bad) b9_log(paste("refused:", b)); stop("block 9 refused to start; see above") }

lock <- file.path(OUT, "_lock")
if (file.exists(lock) && !isTRUE(OPT$force)) {
  pid <- suppressWarnings(as.integer(readLines(lock, n = 1L)))
  if (b9_pid_alive(pid)) stop("another block 9 run (process ", pid, ") holds ", lock)
  b9_log(sprintf("a stale lock from process %s is being replaced", format(pid)))
}
writeLines(c(as.character(Sys.getpid()), format(Sys.time())), lock)
on.exit(if (file.exists(lock)) file.remove(lock), add = TRUE)

plan <- CELLS
plan$done <- vapply(seq_len(nrow(plan)), function(i) file.exists(b9_out_path(as.list(plan[i, ]))), logical(1))
for (i in which(plan$done)) {
  ce <- as.list(plan[i, ])
  b9_log(sprintf("%s: finished file present (%d rows), skipped", ce$cell, nrow(fread(b9_out_path(ce), select = "rep"))))
}
todo <- which(!plan$done)
if (!length(todo)) { b9_log("every selected cell is finished"); quit(save = "no") }

total <- sum(vapply(todo, function(i) b9_reps(as.list(plan[i, ])) * b9_cost(plan$n[i]), numeric(1)))
prog <- b9_prog(total)
cl <- b9_cluster(W)
on.exit(stopCluster(cl), add = TRUE)

k <- 0L
for (i in todo) {
  ce <- as.list(plan[i, ]); k <- k + 1L
  M <- b9_run_cell(cl, ce, b9_reps(ce), prog)
  err <- sum(M[["flag.b9_error"]] %in% 1, na.rm = TRUE)
  deg <- sum(M[["flag.degenerate"]] %in% 1, na.rm = TRUE)
  b9_log(sprintf("[%d/%d] %-26s n=%-5d k=%-3d kept %d, errors %d, degenerate %d  ETA %s",
                 k, length(todo), ce$cell, ce$n, ce$k, nrow(M), err, deg, prog$eta()))
}
b9_log("block 9 finished")
