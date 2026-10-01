## run_M_blockSTREAM.R -- block STREAM: edge.stream() under drift, repeated looks, power and update cost.
## Contract: paper_EDGE/theory/PREDECLARATION_blockSTREAM.md (hashed before any replicate ran).
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
PKG <- Sys.getenv("EBRAHIM_GOF_SRC")
OUT <- edge_battery("STREAM"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
REPS <- 2000L
suppressMessages(library(ebrahim.gof))  # archive: was load_all() of the author's development tree; install ebrahim.gof (>= 2.9.0) from CRAN instead
set.seed(9000001)
p_ref <- plogis(rnorm(10000, -1.5, 1))
BR <- as.numeric(quantile(p_ref, (1:9) / 10, type = 7, names = FALSE))

one <- function(cell, r) {
  set.seed(9100000L + 10000L * cell + r)
  s <- edge.stream(breaks = BR, G = 10)
  mu <- if (cell == 2) -0.5 else -1.5; sg <- if (cell == 2) 1.3 else 1
  slope <- if (cell == 4) 0.8 else 1
  reads <- switch(as.character(cell), "1" = c(10000, 100000), "2" = c(10000, 100000),
                  "3" = seq(1000, 10000, 1000), "4" = c(2000, 10000))
  out <- numeric(0); seen <- 0
  while (seen < max(reads)) {
    p <- plogis(rnorm(100, mu, sg))
    y <- rbinom(100, 1, plogis(slope * qlogis(p)))
    s <- update(s, y, p); seen <- seen + 100
    if (seen %in% reads) out <- c(out, summary(s)$p_value)
  }
  c(cell = cell, rep = r, out, smallest = summary(s)$smallest_group)
}

cl <- parallel::makeCluster(16)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
parallel::clusterExport(cl, c("one", "BR", "PKG"))
invisible(parallel::clusterEvalQ(cl, suppressMessages(library(ebrahim.gof))))  # archive: was load_all() of the author's development tree; install ebrahim.gof (>= 2.9.0) from CRAN instead
res <- list()
for (cell in 1:4) {
  M <- do.call(rbind, parallel::parLapplyLB(cl, seq_len(REPS), function(r, cell) one(cell, r), cell = cell))
  fwrite_csv <- file.path(OUT, sprintf("cell%d_pvalues.csv", cell)); utils::write.csv(M, fwrite_csv, row.names = FALSE)
  pv <- M[, 3:(ncol(M) - 1), drop = FALSE]
  res[[cell]] <- if (cell == 3)
    data.frame(cell = 3, read = c("any of ten looks at 0.05", "any of ten looks at 0.005", "final look only"),
               rate = c(mean(apply(pv < 0.05, 1, any)), mean(apply(pv < 0.005, 1, any)), mean(pv[, 10] < 0.05)),
               smallest_group = median(M[, ncol(M)]))
  else data.frame(cell = cell, read = paste0("n=", if (cell == 4) c(2000, 10000) else c(10000, 100000)),
                  rate = colMeans(pv < 0.05), smallest_group = median(M[, ncol(M)]))
  print(res[[cell]])
}
parallel::stopCluster(cl)

## update cost with a million records already streamed
s <- edge.stream(breaks = BR, G = 10)
for (i in 1:100) { p <- plogis(rnorm(10000, -1.5, 1)); s <- update(s, rbinom(10000, 1, p), p) }
p <- plogis(rnorm(100, -1.5, 1)); y <- rbinom(100, 1, p)
tu <- median(vapply(1:200, function(i) system.time(for (j in 1:100) update(s, y, p))[["elapsed"]] / 100, 0))
ts <- median(vapply(1:200, function(i) system.time(for (j in 1:100) summary(s))[["elapsed"]] / 100, 0))
tim <- data.frame(cell = "timing", read = c("update, 100 patients", "summary"), rate = c(tu, ts), smallest_group = NA)
print(tim)
S <- rbind(do.call(rbind, res), tim)
utils::write.csv(S, file.path(OUT, "_summary.csv"), row.names = FALSE)
cat("written", file.path(OUT, "_summary.csv"), "\n")
