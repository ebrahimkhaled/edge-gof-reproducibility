## analyse_block9c.R -- block 9c judged exactly as its pre-declaration fixed in advance.
## Contract: paper_EDGE/theory/PREDECLARATION_block9c_groupsize.md (sha256 6320ff35...), sections 3 and 4.
## Every rejection is at p <= 0.05; no p-value counts as no rejection; 3 nominal SE at B = 1000 is 0.0207.
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
suppressMessages(library(data.table))
SIMDIR <- edge_path("code/simulations")
DIR <- edge_battery("9c"); OUT <- file.path(DIR, "analysis")
source(file.path(SIMDIR, "_battery_tests.R")); source(file.path(SIMDIR, "_block9_contam.R"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ALPHA <- 0.05; B <- 1000L; SE3 <- 3 * sqrt(ALPHA * (1 - ALPHA) / B)

cells <- data.table(n = rep(c(1000L, 5000L), each = 5), k = rep(c(0L, 5L, 10L, 25L, 50L), 2))
cells[, cell := sprintf("g_n%d_k%02d", n, k)]
cells[, G_extra := ifelse(n == 1000L, "20,25", "25,50,100")]
miss <- cells$cell[!file.exists(file.path(DIR, paste0(cells$cell, "_pvalues.csv.gz")))]
if (length(miss)) stop("block 9c is not finished; missing: ", paste(miss, collapse = ", "))

## long table: one row per cell x G x test
L <- rbindlist(lapply(seq_len(nrow(cells)), function(i) {
  ce <- as.list(cells[i, ]); ce$type <- "full"; ce$ao <- FALSE
  D <- as.data.frame(fread(file.path(DIR, paste0(ce$cell, "_pvalues.csv.gz"))))
  arms <- bt_arms(ce)
  rbindlist(lapply(names(arms), function(a) rbindlist(lapply(
    c("EDGE.poly3.u", "EDGE.sym.u", "EDGE.poly3.sc", "EDGE.sym.sc"), function(base) {
      col <- paste0(base, ".", a)
      if (!col %in% names(D)) return(NULL)
      p <- D[[col]]
      data.table(cell = ce$cell, n = ce$n, k = ce$k, arm = a, G = unname(arms[[a]]),
                 m = ce$n / unname(arms[[a]]), f = ce$k / (ce$n / unname(arms[[a]])),
                 test = base, rejection = mean(is.finite(p) & p <= ALPHA),
                 in_top = mean(D[[paste0("corrupt_in_top.", a)]], na.rm = TRUE))
    }))))
}))
## Stukel uses no partition: one row per cell
St <- rbindlist(lapply(seq_len(nrow(cells)), function(i) {
  ce <- cells[i, ]; D <- as.data.frame(fread(file.path(DIR, paste0(ce$cell, "_pvalues.csv.gz"))))
  data.table(cell = ce$cell, n = ce$n, k = ce$k, test = "Stk.joint",
             rejection = mean(is.finite(D$Stk.joint) & D$Stk.joint <= ALPHA))
}))
fwrite(L, file.path(OUT, "cell_G_test.csv")); fwrite(St, file.path(OUT, "stukel.csv"))
E <- L[test == "EDGE.poly3.u"]
r <- function(nn, kk, GG) E[n == nn & k == kk & G == GG, rejection][1]

cat("3 nominal SE at B = 1000:", sprintf("%.4f", SE3), "\n")

cat("\n=====================  the surface: EDGE-poly3 unit, false-alarm rate  =====================\n")
print(dcast(E, n + k ~ G, value.var = "rejection"), row.names = FALSE)
cat("\n(columns are G; group size m = n/G; every rejection here is a FALSE ALARM)\n")
cat("\nStukel's joint score, which uses no partition, on the same data:\n")
print(dcast(St, n ~ k, value.var = "rejection"), row.names = FALSE)

cat("\n=====================  C9c.5 (control): k = 0 holds size at every G  =====================\n")
ctl <- E[k == 0][, .(n, G, m, rejection, inside = abs(rejection - 0.05) <= SE3)]
print(ctl, row.names = FALSE)
bad_G <- ctl[inside == FALSE]
cat("C9c.5 VERDICT:", if (nrow(bad_G) == 0) "HOLDS at every G" else
      paste("FAILS at", nrow(bad_G), "G values; those G are dropped from the other claims"), "\n")
drop <- if (nrow(bad_G)) bad_G[, paste(n, G)] else character(0)

cat("\n=====================  C9c.1: at n = 5000, k = 25, G = 200 beats G = 10 by >= 0.05  ==========\n")
a <- r(5000, 25, 200); b <- r(5000, 25, 10)
cat(sprintf("G = 200: %.3f   G = 10: %.3f   gap: %+.3f\n", a, b, a - b))
cat("C9c.1 VERDICT:", if ((a - b) >= 0.05) "HOLDS" else "FAILS", "\n")

cat("\n=====================  C9c.2: the rate k/n does not govern  =====================\n")
a <- r(1000, 5, 40); b <- r(5000, 25, 200)
cat(sprintf("both are k/n = 0.005 at the rule G:  n=1000,k=5: %.3f   n=5000,k=25: %.3f   |diff| = %.3f\n",
            a, b, abs(a - b)))
cat("C9c.2 VERDICT:", if (abs(a - b) > SE3) "HOLDS (the same rate behaves differently)" else "FAILS", "\n")

cat("\n=====================  C9c.3: non-decreasing in G  =====================\n")
mono <- rbindlist(lapply(unique(E[k > 0, .(n, k)])[, paste(n, k)], function(s) {
  z <- as.integer(strsplit(s, " ")[[1]]); x <- E[n == z[1] & k == z[2]][order(G)]
  data.table(n = z[1], k = z[2], worst_fall = min(diff(x$rejection)),
             path = paste(sprintf("%.3f", x$rejection), collapse = " -> "))
}))
print(mono, row.names = FALSE)
cat("C9c.3 VERDICT:", if (all(mono$worst_fall >= -SE3)) "HOLDS" else "FAILS", "\n")

cat("\n=====================  C9c.4: H-f (fraction only) against H-fm (fraction and sqrt(m))  =======\n")
prs <- list(list(f = 0.25, a = c(1000, 25, 10), b = c(5000, 25, 50), note = "same f, same m = 100, different n"),
            list(f = 0.50, a = c(1000, 25, 20), b = c(5000, 50, 50), note = "same f, m 50 -> 100"),
            list(f = 1.00, a = c(1000, 25, 40), b = c(5000, 50, 100), note = "same f, m 25 -> 50"))
P <- rbindlist(lapply(prs, function(p) {
  x <- r(p$a[1], p$a[2], p$a[3]); y <- r(p$b[1], p$b[2], p$b[3])
  data.table(f = p$f, note = p$note,
             left = sprintf("n=%d k=%d G=%d (m=%d)", p$a[1], p$a[2], p$a[3], p$a[1] / p$a[3]),
             right = sprintf("n=%d k=%d G=%d (m=%d)", p$b[1], p$b[2], p$b[3], p$b[1] / p$b[3]),
             rate_left = x, rate_right = y, diff = y - x, agree = abs(y - x) <= SE3)
}))
print(P[, .(f, left, right, rate_left, rate_right, diff = round(diff, 3), agree)], row.names = FALSE)
cat("\n", P$note[1], "\n", sep = "")
bigger_m_higher <- P[f > 0.25, all(diff > SE3)]
cat("\nreading fixed in advance:\n")
if (all(P$agree)) cat("  every pair agrees within 3 SE -> H-f SUPPORTED: the corrupted FRACTION of a group governs\n")
if (bigger_m_higher) cat("  the larger-m member is higher in both m-doubling pairs -> H-fm SUPPORTED: tolerance falls as f sqrt(m)\n")
if (!all(P$agree) && !bigger_m_higher) cat("  mixed -> NEITHER is claimed; the measured surface is reported as it stands\n")

cat("\n=====================  C9c.6 (reported): the other forms, and the top group  ==================\n")
print(dcast(L[k %in% c(0, 25) & test %in% c("EDGE.poly3.u", "EDGE.sym.u", "EDGE.poly3.sc", "EDGE.sym.sc")],
            n + k + G ~ test, value.var = "rejection"), row.names = FALSE)
cat("\nmean corrupted records in the highest-risk group:\n")
print(dcast(E[k > 0], n + k ~ G, value.var = "in_top"), row.names = FALSE)
cat("\nwritten:", OUT, "\n")
