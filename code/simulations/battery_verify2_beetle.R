## battery_verify2_beetle.R -- check (3): the E8.5 row orders on the beetle logit fit, against the as-built driver's output
## (Rscript run_M_battery.R --block 6 --cells real_beetle_logit,real_beetle_cloglog --B 1000 --root battery/_review/v2/drv).
##   (a) the stored-order row (rep 0) = the harness on the data as read, and = the dry run's rows 0-10;
##   (b) every order k = 0..1000 recomputed here with the installed package (def.gof, gof_hl, gof_ef, gof_stukel) on the rows
##       in the order sample.int(481) drawn after set.seed(20260914 + k), with tie boundaries counted independently;
##   (c) _row_orders.csv: median and 5%/95% quantiles recomputed per test over the random orders;
##   (d) the range over orders of the ungrouped tests;
##   (e) HL, HL_F and EDGE-poly3 at G = 8 and 10.
## Output: battery/_review/v2/beetle_check.log, beetle_orders_pkg.csv.gz
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
OUT <- edge_battery("_review", "v2")
DRV <- file.path(OUT, "drv", "6")
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))
C <- battery_cells()
cell <- as.list(C[C$block == "6" & C$cell == "real_beetle_logit", ])
P  <- as.data.frame(fread(file.path(DRV, "real_beetle_logit_pvalues.csv.gz")))
PC <- as.data.frame(fread(file.path(DRV, "real_beetle_cloglog_pvalues.csv.gz")))
RO <- as.data.frame(fread(file.path(DRV, "_row_orders.csv")))
META <- c("rep", "seed", "n", "events", "fit_ok", "glm_conv")
tests <- setdiff(names(P), META); tests <- tests[!grepl("^(flag\\.|lam[0-9]|chk\\.|stat\\.|info\\.)", tests)]
mx <- function(a, b) { f <- is.finite(a) & is.finite(b); c(max_abs_diff = if (any(f)) max(abs(a[f] - b[f])) else NA_real_,
                                                           na_mismatch = sum(is.finite(a) != is.finite(b))) }

sink(file.path(OUT, "beetle_check.log"), split = TRUE)
cat("battery_verify2_beetle.R |", format(Sys.time()), "| ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "\n")
cat(sprintf("driver file: %d rows, reps %d-%d consecutive %s; rep 0 seed NA %s; rep k seed = 20260914 + k: %s; rep errors %d\n",
            nrow(P), min(P$rep), max(P$rep), identical(as.integer(P$rep), 0:1000), is.na(P$seed[1]),
            all(P$seed[-1] == 20260914 + 1:1000), sum(P$flag.rep_error)))
cat(sprintf("cloglog refit file: %d rows, same seeds as the logit fit: %s; tie boundaries at G = 10 equal to the logit fit's in every order: %s\n",
            nrow(PC), identical(PC$seed, P$seed), identical(PC$info.tied_boundaries.G10, P$info.tied_boundaries.G10)))

## ---- (a) stored order ----------------------------------------------------------------------------------------------
be <- data.frame(dose = c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839),   # Bliss (1935), 8 doses
                 n = c(59, 60, 62, 56, 63, 59, 62, 60), k = c(6, 13, 18, 28, 52, 53, 61, 60))
d0 <- do.call(rbind, lapply(seq_len(nrow(be)), function(i)
  data.frame(dose = be$dose[i], y = c(rep(1, be$k[i]), rep(0, be$n[i] - be$k[i])))))
cat("\n(a) data typed here identical to the harness's beetle rows:", identical(d0, bt_real_data("beetle_logit")$d), "; n =", nrow(d0), "\n")
cx <- cell; cx$n <- nrow(d0)
h0 <- battery_rep(list(d = d0, f = y ~ dose), cx)
r0 <- unlist(P[P$rep == 0, names(h0)])
cat("    rep 0 of the driver against battery_rep() on the stored order:", paste(names(mx(r0, h0)), mx(r0, h0), collapse = ", "), "\n")
DR <- as.data.frame(fread(edge_battery("dryrun", "6", "real_beetle_logit_pvalues.csv.gz")))
cm <- intersect(names(DR), names(P))
dd <- mx(as.matrix(DR[, cm]), as.matrix(P[match(DR$rep, P$rep), cm]))
cat(sprintf("    dry run (8 workers, reps %d-%d) against the driver (12 workers): %s\n", min(DR$rep), max(DR$rep),
            paste(names(dd), dd, collapse = ", ")))
cat("    stored order, for information beside tab:concord (other implementations, G of the table not re-derived here):\n")
show <- c("EDGE.poly2.u.G8", "EDGE.poly3.u.G8", "EDGE.stk.u.G8", "HLF.G8", "HL.G8", "EDGE.poly2.u.G10", "EDGE.poly3.u.G10",
          "EDGE.stk.u.G10", "HLF.G10", "HL.G10", "HL_w", "Tsiatis", "Stk.LR", "Stk.joint")
cat(paste(sprintf("      %-18s %.4f", show, unlist(P[1, show])), collapse = "\n"), "\n")

## ---- (b) every order with the package ---------------------------------------------------------------------------------
RNGkind("L'Ecuyer-CMRG")
ns <- asNamespace("ebrahim.gof")
Gs <- c(G8 = 8, G10 = 10, Grule = bt_rule_G(nrow(d0)))
tied <- function(ph, grp) { s <- split(ph, grp); k <- length(s); if (k < 2) return(0)
  sum(vapply(s, max, 0)[-k] == vapply(s, min, 0)[-1]) }
rows <- vector("list", 1001)
t0 <- Sys.time()
for (k in 0:1000) {
  o <- if (k == 0) seq_len(nrow(d0)) else { set.seed(20260914 + k); sample.int(nrow(d0)) }
  dk <- d0[o, , drop = FALSE]; rownames(dk) <- NULL
  fit <- suppressWarnings(stats::glm(y ~ dose, family = stats::binomial(), data = dk))
  v <- c(rep = k)
  for (a in names(Gs)) {
    for (b in c("poly3", "poly2", "stukel", "sym")) for (w in c("unit", "score"))
      v[paste0("EDGE.", if (b == "stukel") "stk" else b, ".", if (w == "unit") "u" else "sc", ".", a)] <-
        tryCatch(suppressWarnings(ebrahim.gof::def.gof(fit, G = Gs[[a]], basis = b, weights = w)$p_value), error = function(e) NA_real_)
    ctx <- get(".gof_context", envir = ns)(fit, G = Gs[[a]])
    v[paste0("HL.", a)]  <- as.numeric(get("gof_hl", envir = ns)(ctx)$p_value)
    v[paste0("HLF.", a)] <- as.numeric(get("gof_ef", envir = ns)(ctx)$p_value)
    v[paste0("ties.", a)] <- tied(ctx$ph, get(".gof_groups_ef", envir = ns)(ctx$ph, Gs[[a]]))
  }
  ctx <- get(".gof_context", envir = ns)(fit, G = 10)
  v["Stk.joint"] <- as.numeric(get("gof_stukel", envir = ns)(ctx, list(form = "joint"))$p_value)
  v["Stk.LR"]    <- as.numeric(get("gof_stukel", envir = ns)(ctx, list(form = "lr"))$p_value)
  rows[[k + 1]] <- v
}
K <- as.data.frame(do.call(rbind, rows))
fwrite(K, file.path(OUT, "beetle_orders_pkg.csv.gz"))
cat(sprintf("\n(b) 1001 orders recomputed with the package in %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
pt <- setdiff(names(K), c("rep", grep("^ties\\.", names(K), value = TRUE)))
B <- t(sapply(pt, function(t) mx(K[[t]], P[[t]])))
print(B)
cat("    tests not ok at 1e-8:", sum(!(B[, "na_mismatch"] == 0 & (is.na(B[, "max_abs_diff"]) | B[, "max_abs_diff"] <= 1e-8))), "\n")
for (a in names(Gs))
  cat(sprintf("    tie boundaries %s: independent count = harness info.tied_boundaries in every order: %s (values %s)\n", a,
              identical(as.numeric(K[[paste0("ties.", a)]]), as.numeric(P[[paste0("info.tied_boundaries.", a)]])),
              paste(sort(unique(K[[paste0("ties.", a)]])), collapse = ",")))
cat(sprintf("    distinct fitted risks per order: %s\n", paste(sort(unique(P$info.distinct_risks)), collapse = ",")))
cat(sprintf("    distinct row orders among the 1000 random ones: %d\n",
            length(unique(lapply(1:1000, function(k) { set.seed(20260914 + k); sample.int(nrow(d0)) })))))

## ---- (c) _row_orders.csv ----------------------------------------------------------------------------------------------
R1 <- RO[RO$cell == "real_beetle_logit", ]
cat(sprintf("\n(c) _row_orders.csv: %d rows for real_beetle_logit, one per test: %s; orders column = 1000: %s\n", nrow(R1),
            setequal(R1$test, tests) && !anyDuplicated(R1$test), all(R1$orders == 1000)))
Q <- P[P$rep > 0, ]
rc <- t(sapply(tests, function(t) {
  q <- Q[[t]]; qf <- q[is.finite(q)]; r <- R1[R1$test == t, ]
  own <- c(if (length(qf)) median(qf) else NA, if (length(qf)) quantile(qf, 0.05) else NA, if (length(qf)) quantile(qf, 0.95) else NA,
           P[[t]][P$rep == 0], length(qf), mean(is.finite(q) & q <= 0.05))
  theirs <- c(r$p_median, r$p_q05, r$p_q95, r$p_stored, r$orders_with_p, r$share_p_le_05)
  f <- is.finite(own) & is.finite(theirs)
  c(max_abs_diff = if (any(f)) max(abs(own[f] - theirs[f])) else 0, na_mismatch = sum(is.finite(own) != is.finite(theirs)))
}))
cat("    recomputed median, q05, q95 (type 7, random orders only), stored value, orders with p, share p <= 0.05: max |diff| over all tests",
    max(rc[, 1]), "; NA mismatches", sum(rc[, 2]), "\n")
wq <- t(sapply(tests, function(t) { q <- P[[t]][is.finite(P[[t]])]; r <- R1[R1$test == t, ]
  c(if (length(q)) abs(median(q) - r$p_median) else NA) }))
cat("    (the same median including the stored order would differ by up to", format(max(wq, na.rm = TRUE), digits = 3), ")\n")

## ---- (d) ungrouped tests ----------------------------------------------------------------------------------------------
cat("\n(d) range over all 1001 rows (max - min) and the 5%-95% width in _row_orders.csv\n")
ung <- c("Stk.joint", "Stk.LR", "Stk.sym1", "Stk.marg", "GiViTI", "GiViTI.t50", "Cubic.LR", "HL_w", "PH", "Tsiatis", "Xie", "PR")
for (t in ung) {
  x <- P[[t]]; r <- R1[R1$test == t, ]
  cat(sprintf("    %-11s finite %4d  stored %-12s  max-min %-10s  q95-q05 %s\n", t, sum(is.finite(x)),
              format(x[1], digits = 8), if (any(is.finite(x))) format(diff(range(x, na.rm = TRUE)), digits = 3) else "NA",
              format(r$p_q95 - r$p_q05, digits = 3)))
}
cat("  cloglog refit (link-agnostic tests):\n")
for (t in c("GiViTI", "GiViTI.t50", "HL_w", "PH", "Tsiatis", "HL.G8", "HL.G10", "HLF.G8", "HLF.G10")) {
  x <- PC[[t]]
  cat(sprintf("    %-11s finite %4d  stored %-12s  max-min %s\n", t, sum(is.finite(x)), format(x[1], digits = 8),
              if (any(is.finite(x))) format(diff(range(x, na.rm = TRUE)), digits = 3) else "NA"))
}

## ---- (e) the requested numbers -----------------------------------------------------------------------------------------
cat("\n(e) beetle logit fit, 1000 random orders\n")
rep_t <- c("HL.G8", "HL.G10", "HLF.G8", "HLF.G10", "EDGE.poly3.u.G8", "EDGE.poly3.u.G10", "EDGE.poly3.sc.G8", "EDGE.poly3.sc.G10")
E <- do.call(rbind, lapply(rep_t, function(t) { r <- R1[R1$test == t, ]; q <- Q[[t]]
  data.frame(test = t, stored = r$p_stored, median = r$p_median, q05 = r$p_q05, q95 = r$p_q95, min = min(q, na.rm = TRUE),
             max = max(q, na.rm = TRUE), share_p_le_05 = r$share_p_le_05, orders_with_p = r$orders_with_p,
             tied_boundaries_stored = r$tied_boundaries_stored, share_orders_tied = r$share_orders_tied_boundary) }))
print(E, digits = 4, row.names = FALSE)
fwrite(E, file.path(OUT, "beetle_row_order_table.csv"))
sink()
