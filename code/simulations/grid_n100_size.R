## grid_n100_size.R -- P1-3 addendum: EDGE size at n=100 with G=10 vs G=5 (review request:
## "a size row at n=100 with the drop-to-G=5 recommendation"). EDGE bases only (ms per call),
## single core, base design (quad-family null: x~U(-3,3), eta=0.6x, fit y~x). B=10,000.
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
suppressMessages(library(ebrahim.gof))
set.seed(20260721)
B <- 10000L; n <- 100L
res <- expand.grid(G = c(5L, 10L), basis = c("poly2","poly3","stukel"), stringsAsFactors = FALSE)
res$size <- NA_real_; res$na_rate <- NA_real_
for (i in seq_len(nrow(res))) {
  G <- res$G[i]; basis <- res$basis[i]
  p <- numeric(B)
  set.seed(20260721 + i * 1000L)
  for (b in seq_len(B)) {
    x <- runif(n, -3, 3)
    y <- rbinom(n, 1, plogis(0.6 * x))
    fit <- suppressWarnings(glm(y ~ x, family = binomial()))
    p[b] <- tryCatch(suppressWarnings(edge.gof(fit, G = G, basis = basis))$p_value,
                     error = function(e) NA_real_)
  }
  res$size[i]   <- mean(p < 0.05, na.rm = TRUE)
  res$na_rate[i] <- mean(is.na(p))
  cat(sprintf("G=%2d basis=%-7s size=%.4f na=%.4f\n", G, basis, res$size[i], res$na_rate[i]))
}
write.csv(res, edge_path("code/simulations/sim_n100_size.csv"), row.names = FALSE)
cat("DONE\n")
