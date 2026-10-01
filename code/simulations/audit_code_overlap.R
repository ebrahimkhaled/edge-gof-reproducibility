## audit_code_overlap.R -- code-lens audit (read-only, writes nothing): how many replicates two run-L cells SHARE
## (identical p-values of a partition-free test), and the paired McNemar rows behind claims C3, C5, C7.
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
setwd(edge_path("code/simulations"))
options(width = 220)
rd <- function(t) read.csv(sprintf("runL_%s_pvalues.csv", t), check.names = FALSE)
cat("=== replicate overlap between cells run with the same seeds (identical GiViTI and Stk.LR p-values) ===\n")
for (pr in list(c("base_asym", "base_asym_m10"), c("base_asym", "base_asym_m5"), c("base_asym_m10", "base_asym_m5"),
                c("base_cauchit", "base_cauchit_m10"), c("low_probit", "low_probit_m10"), c("s2_t4", "s2_t4_m10"),
                c("s2_probit", "s2_probit_m25sc"), c("s2_probit", "s2_probit_m10"), c("s2_probit_m25sc", "s2_probit_m10"))) {
  a <- rd(pr[1]); b <- rd(pr[2])
  for (n in intersect(unique(a$n), unique(b$n))) for (lk in unique(a$link)) {
    A <- a[a$n == n & a$link == lk, ]; B <- b[b$n == n & b$link == lk, ]; if (!nrow(A) || !nrow(B)) next
    cat(sprintf("  %-16s vs %-16s n=%-5d %-8s identical GiViTI %.3f  identical Stk.LR %.3f  (reps %d)\n", pr[1], pr[2], n, lk,
      mean(abs(A$GiViTI - B$GiViTI) < 1e-12), mean(abs(A$Stk.LR - B$Stk.LR) < 1e-12), nrow(A)))
  }
}
P <- read.csv("runL_paired.csv")
cat("\n=== paired rows: symmetric-tail links, g.sym / g.sym.sc vs Stk.LR and GiViTI ===\n")
z <- P[P$link %in% c("probit", "cauchit", "t4") & P$probe %in% c("g.sym", "g.sym.sc", "u.sym") & P$rival %in% c("Stk.LR", "GiViTI"), ]
print(z[order(z$link, z$tag, z$n, z$probe, z$rival), c("tag", "link", "c0", "s", "n", "probe", "rival", "pow_probe", "pow_rival", "diff", "discordant", "p", "verdict")], digits = 3, row.names = FALSE)
cat("\n=== paired rows: g.max3.sc vs rivals ===\n")
z <- P[P$probe == "g.max3.sc", ]
print(z[order(z$link, z$tag, z$n, z$rival), c("tag", "link", "n", "rival", "pow_probe", "pow_rival", "diff", "p", "verdict")], digits = 3, row.names = FALSE)
