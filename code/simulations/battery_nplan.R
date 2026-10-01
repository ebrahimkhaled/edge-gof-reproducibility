## battery_nplan.R -- sample sizes of the new battery cells, from the test-agnostic rule of the pre-declaration
## (PREDECLARATION_restructure_battery.md, E0.7). No test is run here and no test's power enters the rule.
##
## For each design x link: draw a population of N = 400,000 from the design, take pi* = the limit of the fitted
## working logistic model y ~ x + d (fitted with the true probabilities as weights, i.e. the fractional response
## p with a quasi-binomial family; the stacked form with weights p and 1 - p is checked to give the same
## coefficients), and D = E[(p - pi*)^2 / (pi* (1 - pi*))]. n_env = 7.85 / D is the 80% power size of an oracle
## one-degree-of-freedom test at alpha = 0.05. n_lo = 0.6 n_env and n_hi = n_env, each rounded to two significant
## figures and kept within [300, 20,000]. A cell with n_env > 40,000 is run once at n = 20,000 as a control.
##
## Also reported: realised event rate and AUC of the true probabilities per design and link, the logistic nulls,
## and whether the E1 intercepts give 12% +/- 0.5% events.
## Writes battery/nplan.csv and battery/nplan.log.
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
OUT    <- edge_battery()
dir.create(OUT, showWarnings = FALSE)
logf <- edge_battery("nplan.log")
sink(logf, split = TRUE)

NPLAN_SEED <- 20260913L
N <- 400000L
N_ORACLE <- 7.85                                   # (qnorm(.975) + qnorm(.8))^2 = 7.849
N_MIN <- 300; N_MAX <- 20000; N_CONTROL_ABOVE <- 40000

LINKS <- list(logit = plogis, probit = pnorm, cauchit = pcauchy, t4 = function(e) pt(e, 4),
              loglog = function(e) exp(-exp(-e)), cloglog = function(e) 1 - exp(-exp(e)))
## E1: intercept giving exactly 12% events at s = 1, per link
C0_E12 <- c(logit = -2.64, probit = -2.01, t4 = -2.18, cauchit = -3.15, loglog = -1.60, cloglog = -2.76)
DESIGNS <- list(base = c(s = 1, c0 = NA), auc = c(s = 2, c0 = NA), e12 = c(s = 1, c0 = NA))
BLOCK2_LINKS <- c("probit", "cauchit", "t4", "loglog")

auc_pop <- function(p) {                           # P(score_case > score_control), score = p, outcome ~ Bern(p)
  o <- order(p); p <- p[o]; w0 <- 1 - p
  cw0 <- cumsum(w0) - w0 / 2                       # ties at identical p count one half
  sum(p * cw0) / (sum(p) * sum(w0))
}
round2 <- function(v) pmin(pmax(signif(v, 2), N_MIN), N_MAX)

cat(sprintf("battery_nplan.R  %s  seed %d  N = %d\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), NPLAN_SEED, N))
cat(sprintf("R %s\n\n", R.version$version.string))

set.seed(NPLAN_SEED)
x <- runif(N, -3, 3); d <- rbinom(N, 1, 0.5)      # one population draw, shared by every design and link
X <- cbind(1, x, d); lp <- 0.6 * x + 0.5 * d

## check once that the fractional-response fit equals the stacked fit with weights p and 1 - p
sub <- 1:20000
p_chk <- pnorm(lp[sub])
b_frac <- suppressWarnings(glm.fit(X[sub, ], p_chk, family = quasibinomial()))$coefficients
b_stack <- suppressWarnings(glm.fit(rbind(X[sub, ], X[sub, ]), c(rep(1, 20000), rep(0, 20000)),
                                    weights = c(p_chk, 1 - p_chk), family = binomial()))$coefficients
cat(sprintf("working-model check (probit, 20,000 points): max |b_frac - b_stack| = %.2e\n\n", max(abs(b_frac - b_stack))))
stopifnot(max(abs(b_frac - b_stack)) < 1e-8)

rows <- list()
for (dz in names(DESIGNS)) for (lk in names(LINKS)) {
  s  <- DESIGNS[[dz]][["s"]]
  c0 <- if (dz == "e12") C0_E12[[lk]] else 0
  if (dz != "e12" && lk == "cloglog") next          # cloglog enters only the e12 intercept check
  p  <- LINKS[[lk]](c0 + s * lp)
  ev <- mean(p); au <- auc_pop(p)
  if (lk == "logit") {
    D <- NA_real_; n_env <- NA_real_
  } else {
    fit <- suppressWarnings(glm.fit(X, p, family = quasibinomial()))
    pis <- fit$fitted.values
    D <- mean((p - pis)^2 / (pis * (1 - pis)))
    n_env <- N_ORACLE / D
  }
  in_b2 <- lk %in% BLOCK2_LINKS
  control <- in_b2 && n_env > N_CONTROL_ABOVE
  n_lo <- if (in_b2 && !control) round2(0.6 * n_env) else NA_real_
  n_hi <- if (in_b2 && !control) round2(n_env) else NA_real_
  n_run <- if (!in_b2) "" else if (control) "20000" else paste(unique(c(n_lo, n_hi)), collapse = ",")
  note <- ""
  if (in_b2 && !control && n_lo == n_hi) note <- "n_lo and n_hi coincide after rounding and clamping: one cell"
  if (in_b2 && !control && (signif(n_env, 2) > N_MAX || signif(0.6 * n_env, 2) < N_MIN))
    note <- paste(note, "clamped to [300, 20000]")
  rows[[length(rows) + 1]] <- data.frame(design = dz, link = lk, s = s, c0 = c0, xdist = "uniform", N = N,
    seed = NPLAN_SEED, event_rate = ev, auc = au, D = D, n_env = n_env, n_lo = n_lo, n_hi = n_hi,
    control = control, in_block2 = in_b2, n_run = n_run,
    e12_ok = if (dz == "e12") abs(ev - 0.12) <= 0.005 else NA, note = trimws(note), stringsAsFactors = FALSE)
}
NP <- do.call(rbind, rows)
write.csv(NP, edge_battery("nplan.csv"), row.names = FALSE)

cat("design  link     s    c0   event   AUC        D      n_env   n_lo   n_hi  run at\n")
for (i in seq_len(nrow(NP))) with(NP[i, ], cat(sprintf("%-6s  %-7s %2g %5.2f  %.4f  %.3f  %.3e  %8.0f  %5s  %5s  %s%s\n",
  design, link, s, c0, event_rate, auc, ifelse(is.na(D), NA, D), ifelse(is.na(n_env), NA, n_env),
  ifelse(is.na(n_lo), "-", format(n_lo)), ifelse(is.na(n_hi), "-", format(n_hi)),
  ifelse(control, "20000 (control)", n_run), ifelse(nzchar(note), paste0("  [", note, "]"), ""))))

cat("\nE1 check, 12% events at s = 1 (target 0.120 +/- 0.005):\n")
e12 <- NP[NP$design == "e12", ]
for (i in seq_len(nrow(e12))) cat(sprintf("  %-7s c0 = %5.2f  event rate %.4f  %s\n", e12$link[i], e12$c0[i],
                                         e12$event_rate[i], ifelse(e12$e12_ok[i], "ok", "OUTSIDE")))
if (!all(e12$e12_ok)) cat("  WARNING: at least one E1 intercept misses 12% +/- 0.5%\n")

cat("\nMatched logistic nulls (block 1b) need these n per design:\n")
for (dz in names(DESIGNS)) {
  z <- NP[NP$design == dz & NP$in_block2, ]
  ns <- sort(unique(as.numeric(unlist(strsplit(z$n_run, ",")))))
  cat(sprintf("  %-5s %s\n", dz, paste(ns, collapse = ", ")))
}

cat("\nText for E7 (to be appended to the pre-declaration by the author; this script does not edit it):\n\n")
cat("| Design | Link | Event rate | AUC | D | n_env | n_lo | n_hi | Cells run at |\n|---|---|---|---|---|---|---|---|---|\n")
for (i in which(NP$in_block2)) with(NP[i, ], cat(sprintf("| %s | %s | %.3f | %.3f | %.3g | %.0f | %s | %s | %s |\n",
  design, link, event_rate, auc, D, n_env, ifelse(is.na(n_lo), "-", format(n_lo)), ifelse(is.na(n_hi), "-", format(n_hi)),
  ifelse(control, "20,000 (control)", n_run))))
cat(sprintf("\nPopulation draw N = %d, seed %d, x ~ U(-3, 3), d ~ Bernoulli(0.5), one draw shared by all designs and links.\n", N, NPLAN_SEED))
sink()
