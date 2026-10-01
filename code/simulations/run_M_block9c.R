## run_M_block9c.R -- block 9c: what governs EDGE's tolerance of corrupted predictions.
## Contract: paper_EDGE/theory/PREDECLARATION_block9c_groupsize.md (sha256 6320ff35...), sections 2 and 3.
##
##   Rscript run_M_block9c.R [--workers 20] [--reps N] [--out battery/9c]
##
## 10 cells: n in {1000, 5000} x k in {0, 5, 10, 25, 50}, B = 1000, C1 corruption (x times 4, outcome drawn
## at the original x), logistic truth so every rejection is a false alarm. Each replicate is tested at every
## declared G in the SAME fit, so the comparison across G is paired on the replicate -- that is the whole
## point of the block: hold the data fixed and move only the partition.
##
## Seeds: set.seed(500000000 + cell_id * 10000 + rep), disjoint from blocks 0-8, 9 (3e8) and 9b (4e8).
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

B9C_SEED0 <- 500000000
B9C_N     <- c(1000L, 5000L)
B9C_K     <- c(0L, 5L, 10L, 25L, 50L)
## section 2: the extra partitions, beyond the G10 and rule-G arms bt_arms() always provides
B9C_EXTRA <- list("1000" = "20,25", "5000" = "25,50,100")

b9c_cells <- function() {
  C <- CJ(n = B9C_N, k = B9C_K, sorted = FALSE)
  C <- C[order(n, k)]
  C[, cell_id := seq_len(.N)]
  C[, seed_base := B9C_SEED0 + cell_id * 10000]
  C[, cell := sprintf("g_n%d_k%02d", n, k)]
  C[, `:=`(block = "9c", B = 1000L, type = "full", ao = FALSE, generator = "b9",
           truth = "logit", corruption = ifelse(k == 0L, "clean", "C1"),
           formula = "y ~ x + d")]
  C[, rate := k / n]
  C[, G_extra := vapply(as.character(n), function(s) B9C_EXTRA[[s]], character(1))]
  as.data.frame(C)
}

## ---- options -----------------------------------------------------------------------------------------
opt <- list(); a <- commandArgs(TRUE); i <- 1L
while (i <= length(a)) { k <- sub("^--", "", a[i]); opt[[k]] <- a[i + 1L]; i <- i + 2L }
W   <- if (is.null(opt$workers)) 20L else as.integer(opt$workers)
OUT <- normalizePath(if (is.null(opt$out)) edge_battery("9c") else opt$out, winslash = "/", mustWork = FALSE)
REPS <- if (is.null(opt$reps)) NA_integer_ else as.integer(opt$reps)
busy <- normalizePath(edge_battery(c(as.character(0:9), "dryrun", "analysis")), winslash = "/", mustWork = FALSE)
if (any(startsWith(paste0(OUT, "/"), paste0(busy, "/")))) stop("--out cannot be a folder of blocks 0-9: ", OUT)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT, "_progress.log")
say <- function(m) { l <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "", m); cat(l, "\n", sep = "")
                     cat(l, "\n", sep = "", file = LOG, append = TRUE) }

C <- b9c_cells()
say(sprintf("block 9c: %d cells, B = %s, %d workers, out %s", nrow(C),
            if (is.na(REPS)) "1000" else paste0(REPS, " (TEST RUN)"), W, OUT))

## refuse before computing: the generator must produce exactly the k the table names
for (i in seq_len(nrow(C))) {
  ce <- as.list(C[i, ])
  g <- b9_gen(ce$truth, ce$n, ce$corruption, ce$rate)
  if (length(g$corrupt) != ce$k) stop(sprintf("%s: generator corrupted %d records, the table says %d",
                                              ce$cell, length(g$corrupt), ce$k))
  if (nrow(g$d) != ce$n) stop(ce$cell, ": wrong number of rows")
}
say("generator checks passed for all 10 cells")

cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, RNGkind("L'Ecuyer-CMRG"))
clusterExport(cl, "SIMDIR", envir = environment())
clusterEvalQ(cl, { source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_block9_contam.R")); TRUE })

t0 <- Sys.time()
for (i in seq_len(nrow(C))) {
  ce <- as.list(C[i, ])
  out <- file.path(OUT, paste0(ce$cell, "_pvalues.csv.gz"))
  if (file.exists(out)) { say(sprintf("%s: finished file present, skipped", ce$cell)); next }
  R <- if (is.na(REPS)) as.integer(ce$B) else REPS
  rows <- clusterApplyLB(cl, seq_len(R), function(r, cell) b9_one(r, cell), cell = ce)
  M <- as.data.frame(do.call(rbind, rows))
  M <- M[order(M$rep), , drop = FALSE]
  tmp <- paste0(out, ".tmp"); fwrite(M, tmp, compress = "gzip")
  if (file.exists(out)) file.remove(out)
  if (!file.rename(tmp, out)) stop("could not rename ", tmp)
  arms <- paste(names(bt_arms(ce)), collapse = ", ")
  say(sprintf("[%d/%d] %-14s n=%-5d k=%-3d arms {%s}  kept %d, errors %d  elapsed %.1f min",
              i, nrow(C), ce$cell, ce$n, ce$k, arms, nrow(M), sum(M[["flag.b9_error"]] %in% 1),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
say("block 9c finished")
