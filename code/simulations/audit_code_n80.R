## audit_code_n80.R -- code-lens audit (read-only, writes nothing): the 80%-power sample size implied by the two simulated
## sample sizes of a run-L cell, for the one-df probes and Stukel's two-df tests, assuming local power
## 1 - F_chi2_k(q_k; c*n) and fitting c to both points. Compared with the scout's local-limit n80.
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
crit <- function(p) { p <- p[is.finite(p)]; as.numeric(quantile(p, .05, type = 1)) }
ncp_of <- function(pw, k) uniroot(function(l) pchisq(qchisq(.95, k), k, ncp = l, lower.tail = FALSE) - pw, c(1e-6, 300))$root
lam80 <- sapply(1:2, function(k) ncp_of(.8, k))
SC <- read.csv("scout_envelope.csv", stringsAsFactors = FALSE)
run <- function(tag, lk, tests, scout) {
  Q <- read.csv(sprintf("runL_%s_pvalues.csv", tag), check.names = FALSE); ns <- sort(unique(Q$n))
  s <- Q$s[1]; c0 <- Q$c0[1]
  for (t in tests) { k <- if (grepl("^Stk", t)) 2 else 1
    pw <- sapply(ns, function(n) { nul <- Q[Q$link == "logit" & Q$n == n, t]; a <- Q[Q$link == lk & Q$n == n, t]; mean(ifelse(is.finite(a), a <= crit(nul), FALSE)) })
    cc <- sapply(seq_along(ns), function(i) ncp_of(pw[i], k) / ns[i])
    sn <- scout[[t]]; sc <- if (is.null(sn)) NA else SC$n80[SC$s == s & SC$c0 == c0 & SC$link == lk & SC$test == sn]
    cat(sprintf("  %-15s %-8s %-14s power %s | implied n80 from each n: %s | pooled %.0f | scout local-limit n80 %s\n", tag, lk, t,
      paste(sprintf("%.3f@%d", pw, ns), collapse = " "), paste(sprintf("%.0f", lam80[k] / cc), collapse = " "), lam80[k] / mean(cc),
      ifelse(is.na(sc), "-", sprintf("%.0f", sc)))) }
}
SCN <- list(g.sym = "sym", u.sym = "sym", g.sym.sc = "sym", u.orc.probit = "matched", u.orc.cauchit = "matched", u.orc.t4 = "matched", Stk.LR = "stukel", Stk.joint = "stukel")
run("base_probit", "probit", c("g.sym", "u.sym", "u.orc.probit", "Stk.LR", "Stk.joint"), SCN)
run("base_cauchit", "cauchit", c("g.sym", "u.sym", "u.orc.cauchit", "Stk.LR", "Stk.joint"), SCN)
run("low_probit", "probit", c("g.sym", "u.sym", "u.orc.probit", "Stk.LR", "Stk.joint"), SCN)
run("s2_probit_m25sc", "probit", c("g.sym.sc", "u.sym", "u.orc.probit", "Stk.LR", "Stk.joint"), SCN)
run("s2_probit_m10", "probit", c("g.sym.sc", "u.sym", "u.orc.probit", "Stk.LR", "Stk.joint"), SCN)
run("s2_t4", "t4", c("g.sym", "u.sym", "u.orc.t4", "Stk.LR", "Stk.joint"), SCN)
