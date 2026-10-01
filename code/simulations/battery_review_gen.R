## battery_review_gen.R -- generator checks (review item d).
##  1. Deterministic: BT_LINKINV links against own closed forms of Stukel's (1988) h-family and the E1 links on a grid;
##     plateau generators against inv_stukel() parsed from _master.R (the July definition) and from _dgp_library.R.
##  2. 200,000 draws through bt_data() for every block 2 (design, link), the exact Stukel links, the plateau links,
##     the logistic nulls, the block 4 skew and sparse cloglog cells: event rate and AUC against nplan.csv (E7) and
##     against exact population values; a 20-bin calibration chi-square of y against own p (does y follow the named link?).
##  3. Plateau cells at n = 1000 re-run for 5 replicates against the stored July p-values (sim_power_broad_pvalues.csv).
## Writes battery/_review/review_gen.log and review_gen.csv.
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
OUT <- edge_battery("_review")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(data.table))
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_battery_cells.R"))
sink(file.path(OUT, "review_gen.log"), split = TRUE)
cat("battery_review_gen.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
RNGkind("L'Ecuyer-CMRG")
Cells <- battery_cells()
NP <- read.csv(edge_battery("nplan.csv"), stringsAsFactors = FALSE)

## ---- own link definitions -------------------------------------------------------------------------------------
t4cdf <- function(t) { u <- t^2 / 4; 0.5 + 0.375 * (t / sqrt(1 + u)) * (1 - (1 / 12) * t^2 / (1 + u)) }
h_stukel <- function(eta, a1, a2) {                  # Stukel (1988) Sec. 2, written from the paper's cases
  out <- numeric(length(eta))
  for (i in seq_along(eta)) {
    e <- eta[i]
    out[i] <- if (e >= 0) {
      if (a1 > 0) (exp(a1 * e) - 1) / a1 else if (a1 < 0) -log(1 - a1 * e) / a1 else e
    } else {
      if (a2 > 0) -(exp(a2 * abs(e)) - 1) / a2 else if (a2 < 0) log(1 - a2 * abs(e)) / a2 else e
    }
  }
  out
}
MYF <- list(logit = function(e) 1 / (1 + exp(-e)), probit = function(e) pnorm(e), cauchit = function(e) 0.5 + atan(e) / pi,
            t4 = t4cdf, loglog = function(e) exp(-exp(-e)), cloglog = function(e) -expm1(-exp(e)),
            stk_long = function(e) 1 / (1 + exp(-h_stukel(e, -1, -1))), stk_short = function(e) 1 / (1 + exp(-h_stukel(e, 1, 1))),
            stk_asym = function(e) 1 / (1 + exp(-h_stukel(e, -1, 1))))

## ---- 1. deterministic -----------------------------------------------------------------------------------------
cat("1. links on a grid eta in [-8, 8] (step 0.001): max |BT_LINKINV - own|\n")
eg <- seq(-8, 8, by = 0.001)
for (lk in names(MYF)) cat(sprintf("  %-9s %.2e\n", lk, max(abs(BT_LINKINV[[lk]](eg) - MYF[[lk]](eg)))))
envM <- new.env()
for (e in parse(file.path(SIMDIR, "_master.R")))
  if (is.call(e) && (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("="))) && is.name(e[[2]]) && as.character(e[[2]]) == "inv_stukel")
    eval(e, envM)
stopifnot(exists("inv_stukel", envir = envM, inherits = FALSE))
eg2 <- seq(-6, 6, by = 1e-4)
own_plateau <- function(e, a1, a2) {                 # the July capped quadratic inverse, written out
  z <- e; pos <- e >= 0
  z[pos] <- (-1 + sqrt(pmax(0, 1 + 2 * a1 * e[pos]))) / a1
  z[!pos] <- (-1 + sqrt(pmax(0, 1 + 2 * a2 * e[!pos]))) / a2
  1 / (1 + exp(-z))
}
for (pl in list(c("plateau_upper", -1, -1), c("plateau_lower", 1, 1), c("plateau_both", -1, 1))) {
  a1 <- as.numeric(pl[2]); a2 <- as.numeric(pl[3])
  cat(sprintf("  %-13s identical to _master.R inv_stukel(%g, %g): %s; max |diff| vs own capped form %.1e; range of p %.4f-%.4f\n",
              pl[1], a1, a2, identical(BT_LINKINV[[pl[1]]](eg2), envM$inv_stukel(eg2, a1, a2)),
              max(abs(BT_LINKINV[[pl[1]]](eg2) - own_plateau(eg2, a1, a2))),
              min(BT_LINKINV[[pl[1]]](eg2)), max(BT_LINKINV[[pl[1]]](eg2))))
}

## ---- exact population values on a midpoint grid (uniform-x designs) --------------------------------------------------
M <- 200000L; xg <- -3 + 6 * (seq_len(M) - 0.5) / M; X2 <- c(xg, xg); D2 <- rep(0:1, each = M)
auc_w <- function(p) {
  o <- order(p); p <- p[o]; key <- match(p, unique(p))
  cs <- tapply(p, key, sum); cn <- tapply(1 - p, key, sum); below <- cumsum(cn) - cn
  sum(cs * (below + cn / 2)) / (sum(p) * sum(1 - p))
}
exact_design <- function(link, s, c0) {
  F <- if (link %in% names(MYF)) MYF[[link]] else BT_LINKINV[[link]]
  p <- F(c0 + s * (0.6 * X2 + 0.5 * D2)); c(event = mean(p), auc = auc_w(p)) }
auc_emp <- function(score, y) {                      # Mann-Whitney with mid-ranks
  r <- rank(score); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0))
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---- 2. draws ---------------------------------------------------------------------------------------------------------
cat("\n2. 200,000 draws through bt_data()\n")
NDRAW <- 200000L
sel <- rbind(
  Cells[Cells$block == "2", ][!duplicated(paste(Cells$design[Cells$block == "2"], Cells$link[Cells$block == "2"])), ],
  Cells[Cells$block == "1b" & Cells$cell %in% c("null_base_n1000", "null_auc_n3500", "null_e12_n1000", "null_skew_n2400", "sparse49_n200"), ],
  Cells[Cells$block == "3" & Cells$cell %in% c("stk_long_n1000", "stk_short_n1000", "stk_asym_n1000", "cauchit_n1000", "t4_n1000", "loglog_n1000",
                                               "link_plateau_upper_n1000", "link_plateau_lower_n1000", "link_plateau_both_n1000"), ],
  Cells[Cells$block == "4" & Cells$cell %in% c("probit_skew_n9400", "cauchit_skew_n2400", "sparse49_cloglog_n200"), ])
rows <- list()
for (i in seq_len(nrow(sel))) {
  ce <- as.list(sel[i, ]); ce$n <- NDRAW
  set.seed(424242L + i); dat <- bt_data(ce); D <- dat$d; y <- D$y
  ## own p from the drawn covariates
  if (ce$generator == "design") {
    lk <- ce$link; eta <- ce$c0 + ce$s * (0.6 * D$x + 0.5 * D$d)
  } else if (ce$generator == "sparse") {
    lk <- ce$link; eta <- ce$intercept + ce$slope * D$x
  } else {                                           # dgp_alt link family: plateau links on the base design
    lk <- PLATEAU_NAME[[ce$param]]; eta <- 0.6 * D$x + 0.5 * D$d
  }
  p <- if (lk %in% names(MYF)) MYF[[lk]](eta) else BT_LINKINV[[lk]](eta)
  own_link <- lk %in% names(MYF)
  ev <- mean(y); se_ev <- sqrt(ev * (1 - ev) / NDRAW)
  au <- auc_emp(p, y)
  ## 20 equal-frequency bins of own p
  br <- unique(quantile(p, seq(0, 1, length.out = 21))); gb <- cut(p, br, include.lowest = TRUE, labels = FALSE)
  O <- tapply(y, gb, sum); E <- tapply(p, gb, sum); V <- tapply(p * (1 - p), gb, sum)
  X2cal <- sum((O - E)^2 / V); pcal <- stats::pchisq(X2cal, length(O), lower.tail = FALSE)
  np <- NP[NP$design == ifelse(is.na(ce$design), "", ce$design) & NP$link == lk, ]
  ex <- if (ce$generator %in% c("design", "dgp_alt") && (is.na(ce$xdist) || ce$xdist == "uniform"))
          exact_design(lk, if (is.na(ce$s)) 1 else ce$s, if (is.na(ce$c0)) 0 else ce$c0) else c(event = NA, auc = NA)
  rows[[i]] <- data.frame(block = ce$block, cell = ce$cell, design = ce$design, link = lk, own_link = own_link,
    event_draw = ev, event_mean_p = mean(p), event_exact = ex[["event"]], event_nplan = if (nrow(np)) np$event_rate else NA,
    z_event_nplan = if (nrow(np)) (ev - np$event_rate) / se_ev else NA,
    z_event_exact = (ev - ex[["event"]]) / se_ev,
    auc_draw = au, auc_exact = ex[["auc"]], auc_nplan = if (nrow(np)) np$auc else NA,
    calib_chisq20 = X2cal, calib_p = pcal, stringsAsFactors = FALSE)
}
R <- do.call(rbind, rows)
fwrite(R, file.path(OUT, "review_gen.csv"))
cat("block cell                        link           event: draw   exact   nplan   z(nplan) | AUC: draw   exact   nplan | calib chi2(20) p\n")
for (i in seq_len(nrow(R))) with(R[i, ], cat(sprintf("%-5s %-27s %-13s %.4f  %s  %s  %6s  | %.4f  %s  %s | %.3f\n",
  block, cell, link, event_draw, ifelse(is.na(event_exact), "  -   ", sprintf("%.4f", event_exact)),
  ifelse(is.na(event_nplan), "  -   ", sprintf("%.4f", event_nplan)), ifelse(is.na(z_event_nplan), "-", sprintf("%.2f", z_event_nplan)),
  auc_draw, ifelse(is.na(auc_exact), "  -   ", sprintf("%.4f", auc_exact)), ifelse(is.na(auc_nplan), "  -   ", sprintf("%.4f", auc_nplan)), calib_p)))
cat(sprintf("\nmax |z| event rate vs nplan %.2f; max |AUC draw - nplan| %.4f; min calibration p %.3f (%d draws)\n",
            max(abs(R$z_event_nplan), na.rm = TRUE), max(abs(R$auc_draw - R$auc_nplan), na.rm = TRUE), min(R$calib_p), nrow(R)))

## ---- 3. plateau cells against the stored July p-values ------------------------------------------------------------------
cat("\n3. plateau cells at n = 1000, replicates 1-5, against sim_power_broad_pvalues.csv\n")
S <- fread(file.path(SIMDIR, "sim_power_broad_pvalues.csv"), select = c("family", "param", "n", "G", "rep", "seed", "test", "p_value"),
           colClasses = list(character = "param"))[family == "link" & grepl("^stukel_", param) & n == 1000L & G == 10L & rep <= 5L]
pairs <- rbind(c("EDGE.poly2.u.G10", "DEF.poly2"), c("EDGE.poly3.u.G10", "DEF.poly3"), c("EDGE.stk.u.G10", "DEF.stukel"),
               c("HL.G10", "HL"), c("HLF.G10", "EF"), c("Stk.marg", "Stukel"))
for (pn in names(PLATEAU_NAME)) {
  ce <- as.list(Cells[Cells$block == "3" & Cells$cell == sprintf("link_%s_n1000", PLATEAU_NAME[[pn]]), ])
  H <- as.data.frame(do.call(rbind, lapply(1:5, battery_one, cell = ce)))
  Sp <- dcast(S[param == pn], rep + seed ~ test, value.var = "p_value")[order(rep)]
  md <- vapply(seq_len(nrow(pairs)), function(j) {
    a <- H[[pairs[j, 1]]]; b <- Sp[[pairs[j, 2]]]; ok <- is.finite(a) & is.finite(b)
    if (any(ok)) max(abs(a[ok] - b[ok])) else NA_real_ }, numeric(1))
  cat(sprintf("  %-13s (%s) seeds equal %s; max|diff| per test: %s\n", PLATEAU_NAME[[pn]], pn, identical(as.numeric(H$seed), as.numeric(Sp$seed)),
              paste(sprintf("%s %.1e", pairs[, 2], md), collapse = ", ")))
}
sink()
