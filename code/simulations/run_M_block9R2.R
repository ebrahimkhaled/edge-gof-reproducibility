## run_M_block9R2.R -- block 9R2: the hybrid partition (coarse tails, fine middle) against the default and the coarse one.
## Contract: paper_EDGE/theory/PREDECLARATION_block9R2_hybrid.md (sha256 9943e551...), frozen before any scenario ran.
##
##   Rscript run_M_block9R2.R [--workers 20] [--reps N] [--out battery/9R2]
##
## The scenarios, the seeds and the regeneration are block 9R's, so that the two tables are read side by side. V0 must
## reproduce the stored EDGE-poly3 unit p-value at the rule G to 1e-8 (the identity gate) and G10 must reproduce block
## 9R's own G10 column exactly; both are counted and reported.
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
source(file.path(SIMDIR, "_block9_contam.R"))
source(file.path(SIMDIR, "_block9R2_hybrid.R"))
BAT <- edge_battery()

r92_cells <- function() {
  B9 <- b9_cells()
  b9keep <- c("logit_clean_n1000", "logit_clean_n5000",
              sprintf("logit_C1_r%s_n1000", c("001", "002", "005", "010")),
              sprintf("logit_C1_r%s_n5000", c("001", "002", "005", "010")),
              "probit_clean_n1000", "probit_clean_n5000", "cloglog_clean_n1000", "cloglog_clean_n5000",
              "probit_C1_r005_n1000", "cloglog_C1_r005_n1000")
  a <- B9[B9$cell %in% b9keep, ]
  A <- data.frame(source = "b9", block = a$block, cell = a$cell, n = a$n, seed_base = a$seed_base,
                  reps = 1000L, stringsAsFactors = FALSE)
  C <- battery_cells()
  batkeep <- c("cauchit_n1000", "t4_n1000", "loglog_n1000", "stk_short_n1000", "stk_long_n1000", "stk_asym_n1000")
  b <- C[C$cell %in% batkeep, ]
  B <- data.frame(source = "bat", block = b$block, cell = b$cell, n = b$n, seed_base = b$seed_base,
                  reps = 1000L, stringsAsFactors = FALSE)
  P <- rbind(A, B)
  if (nrow(P) != 22L) stop("block 9R2: expected 22 scenarios, found ", nrow(P))
  P
}

r92_one <- function(rep, ce, cell) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  set.seed(ce$seed_base + rep)
  dat <- if (ce$source == "b9") { g <- b9_gen(cell$truth, cell$n, cell$corruption, cell$rate); list(d = g$d, f = g$f) }
         else bt_data(cell)
  n <- nrow(dat$d); ev <- sum(dat$d$y)
  fq <- if (min(ev, n - ev) == 0) NULL else bt_fit(dat)
  c(rep = rep, n = n, events = ev, fr2_pvalues(fq))
}
r92_task <- function(rep, ce, cell) r92_one(rep, ce, cell)

OPT <- list(); a <- commandArgs(TRUE); i <- 1L
while (i <= length(a)) { OPT[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
W <- if (is.null(OPT$workers)) 20L else as.integer(OPT$workers)
REPS <- if (is.null(OPT$reps)) NA_integer_ else as.integer(OPT$reps)
OUT <- normalizePath(if (!is.null(OPT$out)) OPT$out else edge_battery(if (is.na(REPS)) "9R2" else "9R2_test"),
                     winslash = "/", mustWork = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lg <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
                      cat(s, "\n"); cat(s, "\n", file = file.path(OUT, "_progress.log"), append = TRUE) }

P <- r92_cells()
lg("block 9R2: %d scenarios, %d workers, out %s%s", nrow(P), W, OUT, if (is.na(REPS)) "" else sprintf(" [%d replicates]", REPS))
cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl))
clusterCall(cl, function(simdir) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  SIMDIR <<- simdir
  source(file.path(simdir, "_battery_tests.R")); source(file.path(simdir, "_battery_cells.R"))
  source(file.path(simdir, "_block9_contam.R")); source(file.path(simdir, "_block9R2_hybrid.R"))
  RNGkind("L'Ecuyer-CMRG"); NULL
}, SIMDIR)
clusterExport(cl, c("r92_one", "r92_task"))

for (i in seq_len(nrow(P))) {
  ce <- as.list(P[i, ]); R <- if (is.na(REPS)) ce$reps else min(REPS, ce$reps)
  f <- file.path(OUT, sprintf("%s_variants.csv.gz", ce$cell))
  if (file.exists(f)) { lg("%s: present, skipped", ce$cell); next }
  t0 <- Sys.time()
  cell <- if (ce$source == "b9") as.list(b9_cells()[b9_cells()$cell == ce$cell, ])
          else as.list(battery_cells()[battery_cells()$cell == ce$cell, ])
  rows <- clusterApplyLB(cl, seq_len(R), r92_task, ce = ce, cell = cell)
  M <- as.data.table(do.call(rbind, rows))
  ## gate 1: V0 is the stored EDGE-poly3 unit p-value at the rule G
  st <- fread(edge_battery(ce$block, sprintf("%s_pvalues.csv.gz", ce$cell)), select = c("rep", "EDGE.poly3.u.Grule"))
  m <- st$EDGE.poly3.u.Grule[match(M$rep, st$rep)]
  bad1 <- sum(!(abs(M$V0 - m) <= 1e-8 | (!is.finite(M$V0) & !is.finite(m))), na.rm = TRUE)
  ## gate 2: G10 is block 9R's own G10 column on the same regenerated data
  f9 <- edge_battery("9R", sprintf("%s_variants.csv.gz", ce$cell))
  bad2 <- NA_integer_
  if (file.exists(f9)) {
    s9 <- fread(f9, select = c("rep", "G10"))
    g <- s9$G10[match(M$rep, s9$rep)]
    bad2 <- sum(!(abs(M$G10 - g) <= 1e-8 | (!is.finite(M$G10) & !is.finite(g))), na.rm = TRUE)
  }
  M[, id_ok := as.integer(abs(V0 - m) <= 1e-8 | (!is.finite(V0) & !is.finite(m)))]
  fwrite(M, f, compress = "gzip")
  lg("[%2d/%d] %-22s n=%-5d R=%d  gates %d/%s  rejection V0 %.3f HYB05 %.3f HYB10 %.3f HYB20 %.3f G10 %.3f  %.1f min",
     i, nrow(P), ce$cell, ce$n, R, bad1, if (is.na(bad2)) "-" else as.character(bad2),
     mean(M$V0 <= 0.05, na.rm = TRUE), mean(M$HYB05 <= 0.05, na.rm = TRUE), mean(M$HYB10 <= 0.05, na.rm = TRUE),
     mean(M$HYB20 <= 0.05, na.rm = TRUE), mean(M$G10 <= 0.05, na.rm = TRUE),
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
lg("block 9R2: run finished")
