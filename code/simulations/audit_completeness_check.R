## audit_completeness_check.R  (completeness critic, 2026-09-13 night; read-only on run L files)
## (1) the symmetric probe against the paper's OWN bases (EDGE.stk, EDGE.poly3) and against the
##     better of Stukel's two proper forms, size-adjusted AND at the nominal 0.05 cut-off;
## (2) n80 for the probit cells under several interpolation rules (resolves the code/numbers
##     referees' disagreement on C6).
## Writes audit_completeness_check.csv; prints the n80 table. No cluster needed (CSV reads only).
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
tags <- c("base_probit", "base_cauchit", "base_cauchit_m10", "s2_cauchit", "s2_t4", "s2_t4_m10",
          "s2_probit", "s2_probit_m25sc", "s2_probit_m10", "low_probit", "low_probit_m10")

mcn <- function(a, b) {
  d1 <- sum(a & !b); d2 <- sum(!a & b)
  p <- if (d1 + d2 == 0) 1 else binom.test(d1, d1 + d2)$p.value
  c(d1 = d1, d2 = d2, p = p)
}

out <- list(); pw_store <- list()
for (tg in tags) {
  f <- sprintf("runL_%s_pvalues.csv", tg)
  if (!file.exists(f)) { cat("missing", f, "\n"); next }
  P <- read.csv(f, check.names = FALSE)
  alt <- unique(P[P$link != "logit", c("link", "s", "c0", "n")])
  for (i in seq_len(nrow(alt))) {
    a <- alt[i, ]
    A <- P[P$link == a$link & P$s == a$s & P$c0 == a$c0 & P$n == a$n, ]
    N <- P[P$link == "logit" & P$s == a$s & P$c0 == a$c0 & P$n == a$n, ]
    if (nrow(N) == 0) next
    orc <- paste0("u.orc.", a$link)
    tests <- intersect(c("g.sym", "g.sym.sc", "u.sym", orc, "g.max3", "g.max3.sc", "EDGE.poly3",
                         "EDGE.stk", "EDGE.stk.sc", "Stk.joint", "Stk.LR", "GiViTI"), names(P))
    adj <- list(); nom <- list(); pw <- c(); pn <- c(); sz <- c()
    for (t in tests) {
      crit <- quantile(N[[t]], 0.05, type = 1, na.rm = TRUE)
      adj[[t]] <- !is.na(A[[t]]) & A[[t]] <= crit
      nom[[t]] <- !is.na(A[[t]]) & A[[t]] <= 0.05
      pw[t] <- mean(adj[[t]]); pn[t] <- mean(nom[[t]]); sz[t] <- mean(N[[t]] <= 0.05, na.rm = TRUE)
    }
    probe <- if ("g.sym.sc" %in% tests) "g.sym.sc" else "g.sym"
    bestS <- if (pw["Stk.LR"] >= pw["Stk.joint"]) "Stk.LR" else "Stk.joint"
    bestSn <- if (pn["Stk.LR"] >= pn["Stk.joint"]) "Stk.LR" else "Stk.joint"
    rivals <- intersect(c("EDGE.stk", "EDGE.poly3", "EDGE.stk.sc", "u.sym", "GiViTI"), tests)
    for (rv in c(rivals, "bestStukel")) {
      ra <- if (rv == "bestStukel") bestS else rv
      rn <- if (rv == "bestStukel") bestSn else rv
      ma <- mcn(adj[[probe]], adj[[ra]]); mn <- mcn(nom[[probe]], nom[[rn]])
      out[[length(out) + 1]] <- data.frame(tag = tg, link = a$link, s = a$s, c0 = a$c0, n = a$n,
        B = nrow(A), probe = probe, rival = rv, rival_used_adj = ra,
        probe_adj = pw[probe], rival_adj = pw[ra], diff_adj = pw[probe] - pw[ra], p_adj = ma["p"],
        probe_nom = pn[probe], rival_nom = pn[rn], diff_nom = pn[probe] - pn[rn], p_nom = mn["p"],
        size_probe = sz[probe], size_rival = sz[rn])
    }
    pw_store[[paste(tg, a$link, a$s, a$c0, a$n)]] <- list(tag = tg, link = a$link, s = a$s, c0 = a$c0,
                                                          n = a$n, pw = pw)
  }
}
res <- do.call(rbind, out); rownames(res) <- NULL
write.csv(res, "audit_completeness_check.csv", row.names = FALSE)
options(width = 200)
cat("\n=== symmetric probe minus rival (adj = size-adjusted; nom = p<=0.05) ===\n")
print(transform(res[, c("tag", "link", "s", "c0", "n", "probe", "rival", "rival_used_adj", "probe_adj",
                        "rival_adj", "diff_adj", "p_adj", "diff_nom", "p_nom")],
                p_adj = signif(p_adj, 2), p_nom = signif(p_nom, 2),
                probe_adj = round(probe_adj, 3), rival_adj = round(rival_adj, 3),
                diff_adj = round(diff_adj, 3), diff_nom = round(diff_nom, 3)), row.names = FALSE)

## ---- n80 by interpolation rule, where a design has two sample sizes -----------------------
lam80 <- function(k) uniroot(function(l) pchisq(qchisq(0.95, k), k, ncp = l, lower.tail = FALSE) - 0.8,
                             c(0.1, 60))$root
ncp_of <- function(pow, k) uniroot(function(l) pchisq(qchisq(0.95, k), k, ncp = l, lower.tail = FALSE) - pow,
                                   c(1e-6, 80))$root
n80_rules <- function(n1, n2, p1, p2, k) {
  z1 <- qnorm(p1); z2 <- qnorm(p2)
  sq <- (sqrt(n1) + (qnorm(0.8) - z1) / (z2 - z1) * (sqrt(n2) - sqrt(n1)))^2
  lin <- n1 + (0.8 - p1) / (p2 - p1) * (n2 - n1)
  if (is.na(k)) return(c(probit_sqrt_n = sq, linear_n = lin, local_pt1 = NA, local_pt2 = NA, local_geo = NA))
  l1 <- lam80(k) / (ncp_of(p1, k) / n1); l2 <- lam80(k) / (ncp_of(p2, k) / n2)
  c(probit_sqrt_n = sq, linear_n = lin, local_pt1 = l1, local_pt2 = l2, local_geo = sqrt(l1 * l2))
}
dfk <- c(g.sym = 1, g.sym.sc = 1, u.sym = 1, u.orc.probit = 1, Stk.LR = 2, Stk.joint = 2,
         EDGE.stk = NA, EDGE.poly3 = NA, GiViTI = NA, g.max3.sc = NA, EDGE.stk.sc = 2)
cat("\n=== n80 by interpolation rule (probit designs) ===\n")
for (tg in c("base_probit", "s2_probit_m25sc", "s2_probit_m10", "low_probit")) {
  keys <- names(pw_store)[sapply(pw_store, function(z) z$tag == tg && z$link == "probit")]
  if (length(keys) < 2) next
  ns <- sapply(pw_store[keys], `[[`, "n"); o <- order(ns); keys <- keys[o]
  k1 <- pw_store[[keys[1]]]; k2 <- pw_store[[keys[length(keys)]]]
  for (t in intersect(names(dfk), names(k1$pw))) {
    r <- n80_rules(k1$n, k2$n, k1$pw[t], k2$pw[t], dfk[t])
    cat(sprintf("%-16s %-12s n=%5d/%5d pow %.3f/%.3f | sqrt-probit %6.0f | linear %6.0f | local pt1 %6.0f pt2 %6.0f geo %6.0f\n",
                tg, t, k1$n, k2$n, k1$pw[t], k2$pw[t], r[1], r[2], r[3], r[4], r[5]))
  }
}
cat("\nlambda80(1) =", round(lam80(1), 4), " lambda80(2) =", round(lam80(2), 4), "\n")
## local-limit power of the scout's matched direction and sym probe at n = 16,000 (base design)
for (nn in c(15485, 15558, 14955, 15870)) cat(sprintf("n80=%5d -> local power at n=16000: %.3f\n", nn,
  pchisq(qchisq(0.95, 1), 1, ncp = lam80(1) * 16000 / nn, lower.tail = FALSE)))
