## null_calibration_checks.R  — EXPERIMENT 2.6 (SPECKIT §2.6)
## Goodness-of-null (VALID: fixed POPULATION lambda + Imhof PIT) + parametric-bootstrap
## vs Omega-calibration agreement + Imhof-vs-Satterthwaite tail comparison.
##
## Contract: source _harness.R then source(DGP). Battery/DEF go through the installed
## ebrahim.gof package (workers library() it). We re-derive the DEF Omega/Z/S/lambda
## EXACTLY as in ebrahim.gof::def.gof (R/def_gof.R: equal-frequency groups by fitted p,
## r=(o-e)/sqrt(Vg), Omega=I-U(X'WX)^-1U', S=(Z'r)'(Z'Z)^-1(Z'r),
## lam=eigen((Z'Z)^-1 Z'Omega Z)) so the population lambda spectrum and the per-sample S
## are computed from the same construction the package ships.
##
## NOTE ON lambda (SPECKIT §2.6a / REVIEWER #8): the reference eigenvalues lambda_pop are
## POPULATION values from a SINGLE n_cal (=2e5 full / 1e4 smoke) correctly-specified draw,
## NOT per-sample-estimated. The PIT u_b = 1 - imhof(S_b, lambda_pop)$Qq is therefore valid
## Uniform(0,1) under H0 (no parameters re-estimated per sample), so KS/AD are honest.
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

suppressMessages({
  library(parallel)
})
SIMDIR <- edge_path("code/simulations")
source(file.path(SIMDIR, "_harness.R"))   # SIMDIR, DGP, EK_NCORES, ek_reps, ek_cluster, ek_run_cell, ek_append, ek_mcse
source(DGP)                               # dgp_null/dgp_alt, sqb/sib/scb, inv_stukel, EK_CURATED, ...

suppressMessages({
  library(ebrahim.gof)
  library(CompQuadForm)                   # imhof
})

STOP_LOG <- file.path(SIMDIR, "_STOP_TRIGGERS.log")
ek_stop_trigger <- function(trigger, detail) {
  line <- sprintf("%s\t%s\t%s\t%s", format(Sys.time()), "null_calibration_checks.R", trigger, detail)
  cat(line, "\n", file = STOP_LOG, append = TRUE)
  cat("STOP-TRIGGER:", line, "\n")
}

## smoke knobs (self-test uses these): REPS_SCALE tiny + N_CAL_SMOKE flag
SMOKE   <- as.numeric(Sys.getenv("REPS_SCALE", "1")) < 0.1
N_CAL   <- if (SMOKE) 1e4 else 2e5      # capped 1e6->2e5 for memory; eigenvalue spectrum is population-scale + identical
B_KS    <- ek_reps(2000)                # SPECKIT: B=2000 null samples per (n,G,basis); smoke ~ B=100 via REPS_SCALE
BASES   <- c("poly2", "poly3", "stukel")
N_GRID  <- if (SMOKE) c(500L) else c(500L, 1000L, 5000L)
G_FIX   <- 10L

## ------------------------------------------------------------------------------------
## Core: replicate ebrahim.gof::def.gof internals to get (S, lambda) for a fitted glm.
## Returns list(S=, lam=) with lam = eigenvalues of (Z'Z)^-1 Z'Omega Z (lam>1e-9), and
## the group-mean-prob vector pbar so a FIXED reference basis/lambda can be reused.
## This mirrors def_gof.R verbatim in construction (Omega, Z, S, lam).
## ------------------------------------------------------------------------------------
.def_basis_local <- function(pbar, basis) {
  if (basis %in% c("poly2", "poly3")) {
    deg <- if (basis == "poly2") 2L else 3L
    if (length(unique(round(pbar, 8))) < deg + 1) return(NULL)
    as.matrix(stats::poly(pbar, deg))
  } else {
    e <- stats::qlogis(pbar)
    cbind(e, e^2 * (e >= 0), -e^2 * (e < 0))
  }
}

## Compute (S, lam) for one fitted glm and one basis, exactly as def.gof does.
def_S_lam <- function(fit, basis, G = 10L) {
  y   <- as.numeric(fit$y)
  ph  <- pmin(pmax(as.numeric(stats::fitted(fit)), 1e-6), 1 - 1e-6)
  eta <- as.numeric(stats::predict(fit, type = "link"))
  dmu <- fit$family$mu.eta(eta)
  X   <- stats::model.matrix(fit)
  n   <- length(y)
  V   <- ph * (1 - ph); w <- dmu^2 / V
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  idx <- split(seq_len(n), grp); Gn <- length(idx)
  og   <- vapply(idx, function(I) sum(y[I]),   numeric(1))
  eg   <- vapply(idx, function(I) sum(ph[I]),  numeric(1))
  Vg   <- vapply(idx, function(I) sum(V[I]),   numeric(1))
  pbar <- vapply(idx, function(I) mean(ph[I]), numeric(1))
  r    <- (og - eg) / sqrt(Vg)
  U    <- t(vapply(idx, function(I) colSums(dmu[I] * X[I, , drop = FALSE]), numeric(ncol(X)))) / sqrt(Vg)
  Omega <- diag(Gn) - U %*% solve(crossprod(X, w * X)) %*% t(U)
  Z <- .def_basis_local(pbar, basis)
  if (is.null(Z)) return(NULL)
  Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]
  if (ncol(Z) < 1) return(NULL)
  ZtZ <- crossprod(Z); Zr <- crossprod(Z, r)
  S   <- as.numeric(t(Zr) %*% solve(ZtZ) %*% Zr)
  lam <- Re(eigen(solve(ZtZ) %*% (t(Z) %*% Omega %*% Z), only.values = TRUE)$values)
  lam <- lam[lam > 1e-9]
  list(S = S, lam = lam)
}

## Imhof upper-tail P(Q > q) with fixed lambda; PIT u = 1 - Qq (lower-tail CDF at S).
imhof_upper <- function(q, lam) {
  p <- tryCatch(CompQuadForm::imhof(q, lam)$Qq, error = function(e) NA_real_)
  min(max(p, 0), 1)
}
imhof_pit <- function(q, lam) 1 - imhof_upper(q, lam)   # F_ref(S) = lower-tail CDF

## Satterthwaite scaled-chi2 upper tail (def_gof.R .def_pvalue satterthwaite branch).
satter_upper <- function(q, lam) {
  cc <- sum(lam^2) / sum(lam); nu <- sum(lam)^2 / sum(lam^2)
  stats::pchisq(q / cc, df = nu, lower.tail = FALSE)
}

## ------------------------------------------------------------------------------------
## Anderson-Darling statistic + p-value of {u} against Uniform(0,1).
## goftest is not installed here, so implement the standard uniformity AD directly:
##   A2 = -m - (1/m) sum_{i=1}^m (2i-1)[ln u_(i) + ln(1 - u_(m+1-i)) ]
## p-value via Marsaglia & Marsaglia (2004) AD CDF (parameter-free / fully specified case).
## ------------------------------------------------------------------------------------
ad_unif_stat <- function(u) {
  u <- sort(u); m <- length(u)
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  i <- seq_len(m)
  A2 <- -m - sum((2 * i - 1) * (log(u) + log(1 - rev(u)))) / m
  A2
}
## Marsaglia-Marsaglia adinf() + errfix() : CDF of A2 for the fully-specified case.
.ad_adinf <- function(z) {
  if (z < 2) return(exp(-1.2337141 / z) / sqrt(z) *
      (2.00012 + (0.247105 - (0.0649821 - (0.0347962 - (0.011672 - 0.00168691 * z) * z) * z) * z) * z))
  exp(-exp(1.0776 - (2.30695 - (0.43424 - (0.082433 - (0.008056 - 0.0003146 * z) * z) * z) * z) * z))
}
.ad_errfix <- function(n, x) {
  if (x > 0.8) return((-130.2137 + (745.2337 - (1705.091 - (1950.646 - (1116.360 - 255.7844 * x) * x) * x) * x) * x) / n)
  cc <- 0.01265 + 0.1757 / n
  if (x < cc) { t <- x / cc; t <- sqrt(t) * (1 - t) * (49 * t - 102)
    return(t * (0.0037 / n^2 + 0.00078 / n + 0.00006) / n) }
  t <- (x - cc) / (0.8 - cc); t <- -0.00022633 + (6.54034 - (14.6538 - (14.458 - (8.259 - 1.91864 * t) * t) * t) * t) * t
  t * (0.04213 + 0.01365 / n) / n
}
ad_unif_pvalue <- function(A2, m) {
  cdf <- .ad_adinf(A2) + .ad_errfix(m, .ad_adinf(A2))
  p <- 1 - cdf
  min(max(p, 0), 1)
}

## ====================================================================================
## (a) VALID goodness-of-null: fixed POPULATION lambda + Imhof PIT -> KS + AD vs Uniform.
## ====================================================================================
run_goodness_of_null <- function(cl) {
  ref_path <- file.path(SIMDIR, "null_ref_lambda.csv")
  ks_path  <- file.path(SIMDIR, "null_ks_table.csv")
  if (file.exists(ref_path)) file.remove(ref_path)
  if (file.exists(ks_path))  file.remove(ks_path)

  ## family used for the null geometry: quad-null (y~x, logit) is the canonical DEF null.
  ## lambda_pop is computed ONCE per (n,G,basis) from ONE n_cal correctly-specified draw.
  cell_id <- 0L
  for (n in N_GRID) {
    for (basis in BASES) {
      cell_id <- cell_id + 1L
      ## ---- (a.1) POPULATION lambda from ONE big correctly-specified draw ----
      set.seed(20260713L + cell_id)      # deterministic reference draw (population spectrum)
      G0  <- dgp_null("quad", N_CAL)
      fit0 <- suppressWarnings(glm(G0$f, data = G0$d, family = binomial()))
      sl0 <- def_S_lam(fit0, basis, G = G_FIX)
      if (is.null(sl0) || length(sl0$lam) == 0) {
        ek_stop_trigger("null_ref_lambda_degenerate",
                        sprintf("n=%d basis=%s: degenerate basis/lambda on n_cal draw", n, basis))
        next
      }
      lam_pop <- sl0$lam
      ek_append(data.frame(n = n, G = G_FIX, basis = basis,
                           j = seq_along(lam_pop), lambda_pop = lam_pop),
                ref_path)

      ## ---- (a.2) B null samples at target n -> S_b -> u_b = imhof PIT at FIXED lam_pop ----
      seed_base <- 20260713L + cell_id * 100000L
      one_rep <- function(rep, cell) {
        Gs  <- dgp_null("quad", cell$n)
        fit <- suppressWarnings(glm(Gs$f, data = Gs$d, family = binomial()))
        sl  <- def_S_lam(fit, cell$basis, G = cell$G)
        S   <- if (is.null(sl)) NA_real_ else sl$S
        ## return >=2 named elements so ek_run_cell preserves column names (length-1
        ## returns get simplified by parSapply to an unnamed vector -> column "V1").
        c(S = S, ok = as.numeric(is.finite(S)))
      }
      clusterExport(cl, c("def_S_lam", ".def_basis_local"), envir = environment())
      cell <- list(n = n, basis = basis, G = G_FIX)
      df   <- ek_run_cell(cl, B_KS, seed_base, one_rep, cell)
      Sb   <- df$S[is.finite(df$S)]
      ## PIT with the FIXED population lambda (computed on the master, exported for speed)
      ub   <- vapply(Sb, function(s) imhof_pit(s, lam_pop), numeric(1))
      ub   <- ub[is.finite(ub)]
      Beff <- length(ub)
      if (Beff < 10) {
        ek_stop_trigger("null_pit_empty",
                        sprintf("n=%d basis=%s: only %d finite PIT values", n, basis, Beff))
        next
      }
      ks   <- suppressWarnings(stats::ks.test(ub, "punif"))
      A2   <- ad_unif_stat(ub); adp <- ad_unif_pvalue(A2, Beff)
      ek_append(data.frame(n = n, G = G_FIX, basis = basis,
                           ks_stat = as.numeric(ks$statistic), ks_p = as.numeric(ks$p.value),
                           ad_stat = A2, ad_p = adp, B = Beff,
                           lambda_source = "pop_2e5"),
                ks_path)
      cat(sprintf("[gon] n=%d basis=%s B=%d  KS=%.4f (p=%.3f)  AD=%.4f (p=%.3f)  min/maxlam=%.3g\n",
                  n, basis, Beff, ks$statistic, ks$p.value, A2, adp, min(lam_pop) / max(lam_pop)))
    }
  }
  invisible(list(ref = ref_path, ks = ks_path))
}

## ====================================================================================
## (b) parametric-bootstrap vs Omega-calibration |Delta p| on >=4 sim cells + 3 real sets.
##  Omega p  = def.gof(fit)$p_value (analytic weighted-chi2 at per-sample estimated lambda).
##  Boot  p  = parametric bootstrap: refit under H0, simulate y* ~ Bernoulli(phat), recompute S*,
##             p_boot = mean(S* >= S_obs).  |Delta p| = |p_omega - p_boot|.
## ====================================================================================
run_bootstrap_agreement <- function(cl) {
  out_path <- file.path(SIMDIR, "null_bootstrap_agreement.csv")
  if (file.exists(out_path)) file.remove(out_path)
  Bboot <- if (SMOKE) 100L else 999L

  ## def_S_lam variant that takes an explicit (y, model.matrix) so the bootstrap can
  ## refit via glm.fit on the already-evaluated design (handles transformed terms like
  ## log(Volume) in real data without re-applying the transform).
  def_S_lam_mm <- function(y, ph, dmu, X, basis, G = 10L) {
    n <- length(y); V <- ph * (1 - ph); w <- dmu^2 / V
    grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
    idx <- split(seq_len(n), grp); Gn <- length(idx)
    og   <- vapply(idx, function(I) sum(y[I]),   numeric(1))
    eg   <- vapply(idx, function(I) sum(ph[I]),  numeric(1))
    Vg   <- vapply(idx, function(I) sum(V[I]),   numeric(1))
    pbar <- vapply(idx, function(I) mean(ph[I]), numeric(1))
    r    <- (og - eg) / sqrt(Vg)
    U    <- t(vapply(idx, function(I) colSums(dmu[I] * X[I, , drop = FALSE]), numeric(ncol(X)))) / sqrt(Vg)
    Omega <- diag(Gn) - U %*% solve(crossprod(X, w * X)) %*% t(U)
    Z <- .def_basis_local(pbar, basis); if (is.null(Z)) return(NULL)
    Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]; if (ncol(Z) < 1) return(NULL)
    ZtZ <- crossprod(Z); Zr <- crossprod(Z, r); S <- as.numeric(t(Zr) %*% solve(ZtZ) %*% Zr)
    lam <- Re(eigen(solve(ZtZ) %*% (t(Z) %*% Omega %*% Z), only.values = TRUE)$values)
    lam <- lam[lam > 1e-9]
    list(S = S, lam = lam)
  }

  ## one boot-vs-Omega |Delta p| for a given fitted glm + basis
  boot_vs_omega <- function(fit, basis, G, Bboot) {
    sl <- def_S_lam(fit, basis, G = G)
    if (is.null(sl)) return(c(p_omega = NA_real_, p_boot = NA_real_, S_obs = NA_real_))
    ## analytic Omega p (Imhof at per-sample lambda -- the object def.gof reports)
    p_om <- imhof_upper(sl$S, sl$lam)
    ## parametric bootstrap under the fitted (null) model: resimulate y* ~ Bern(phat),
    ## refit via glm.fit on the frozen design X (mu.eta from the binomial-logit family).
    X   <- stats::model.matrix(fit)
    fam <- binomial()
    Sstar <- numeric(Bboot)
    for (b in seq_len(Bboot)) {
      ystar <- rbinom(nrow(X), 1, pmin(pmax(as.numeric(stats::fitted(fit)), 1e-6), 1 - 1e-6))
      fb <- suppressWarnings(tryCatch(stats::glm.fit(X, ystar, family = fam),
                                      error = function(e) NULL))
      if (is.null(fb) || !isTRUE(fb$converged)) { Sstar[b] <- NA_real_; next }
      phb <- pmin(pmax(fb$fitted.values, 1e-6), 1 - 1e-6)
      etab <- fam$linkfun(phb); dmub <- fam$mu.eta(etab)
      slb <- def_S_lam_mm(ystar, phb, dmub, X, basis, G = G)
      Sstar[b] <- if (is.null(slb)) NA_real_ else slb$S
    }
    Sstar <- Sstar[is.finite(Sstar)]
    p_bt <- mean(Sstar >= sl$S)
    c(p_omega = p_om, p_boot = p_bt, S_obs = sl$S)
  }

  add_row <- function(source, label, n, basis, res) {
    ek_append(data.frame(source = source, dataset = label, n = n, basis = basis,
                         p_omega = res["p_omega"], p_boot = res["p_boot"],
                         abs_delta_p = abs(res["p_omega"] - res["p_boot"]),
                         Bboot = Bboot, row.names = NULL), out_path)
    cat(sprintf("[boot] %-10s %-14s n=%d %-6s  p_omega=%.4f p_boot=%.4f |dp|=%.4f\n",
                source, label, n, basis, res["p_omega"], res["p_boot"],
                abs(res["p_omega"] - res["p_boot"])))
  }

  ## --- >=4 simulation cells (correctly-specified nulls across families/n) ---
  sim_cells <- list(
    list(fam = "quad",    n = 500L),
    list(fam = "binint",  n = 500L),
    list(fam = "contint", n = 1000L),
    list(fam = "link",    n = 1000L)
  )
  for (i in seq_along(sim_cells)) {
    cc <- sim_cells[[i]]
    set.seed(20260714L + i)
    G0 <- dgp_null(cc$fam, cc$n)
    fit <- suppressWarnings(glm(G0$f, data = G0$d, family = binomial()))
    res <- boot_vs_omega(fit, "poly3", G_FIX, Bboot)
    add_row("sim", cc$fam, cc$n, "poly3", res)
  }

  ## --- 3 real datasets (birthwt / vaso / GLOW; datasets::infert as a fallback) ---
  ## Lazy-loaded package datasets are pulled via data(..., envir=e), NOT asNamespace get.
  load_ds <- function(name, pkg) {
    e <- new.env()
    ok <- tryCatch({ utils::data(list = name, package = pkg, envir = e); TRUE },
                   error = function(er) FALSE)
    if (!ok || !exists(name, envir = e)) return(NULL)
    get(name, envir = e)
  }
  real_fits <- list()
  if (requireNamespace("MASS", quietly = TRUE)) {
    bw <- load_ds("birthwt", "MASS")
    if (!is.null(bw)) { bw$low <- as.integer(bw$low)
      real_fits[["birthwt"]] <- suppressWarnings(glm(low ~ age + lwt + smoke + ht + ui, data = bw, family = binomial())) }
  }
  if (requireNamespace("robustbase", quietly = TRUE)) {
    vs <- load_ds("vaso", "robustbase")
    if (!is.null(vs))
      real_fits[["vaso"]] <- suppressWarnings(glm(Y ~ log(Volume) + log(Rate), data = vs, family = binomial()))
  }
  if (requireNamespace("aplore3", quietly = TRUE)) {
    gl <- load_ds("glow500", "aplore3")
    if (!is.null(gl)) { gl$fracture01 <- as.integer(gl$fracture == "Yes")
      real_fits[["glow"]] <- suppressWarnings(glm(fracture01 ~ age + weight + priorfrac + premeno + raterisk,
                                                  data = gl, family = binomial())) }
  }
  ## fallback 3rd real set if a source is absent: datasets::infert (case ~ ...)
  if (length(real_fits) < 3) {
    inf <- load_ds("infert", "datasets")
    if (!is.null(inf))
      real_fits[["infert"]] <- suppressWarnings(glm(case ~ age + parity + induced + spontaneous, data = inf, family = binomial()))
  }
  for (nm in names(real_fits)) {
    fit <- real_fits[[nm]]
    res <- boot_vs_omega(fit, "poly3", G_FIX, Bboot)
    add_row("real", nm, length(fit$y), "poly3", res)
  }
  invisible(out_path)
}

## ====================================================================================
## (c) Imhof vs Satterthwaite p-value scatter over realistic lambda spectra (incl one
##     near-zero eigenvalue). max|Delta p| in the small-p tail -> fallback rule.
## ====================================================================================
run_imhof_satterthwaite <- function() {
  out_path <- file.path(SIMDIR, "imhof_satterthwaite.csv")
  if (file.exists(out_path)) file.remove(out_path)

  ## realistic DEF spectra: 1-3 df. Include a near-zero eigenvalue case.
  spectra <- list(
    "poly2_typical"   = c(1.00, 0.62),
    "poly3_typical"   = c(1.00, 0.55, 0.28),
    "stukel_typical"  = c(1.00, 0.47, 0.12),
    "near_zero_1"     = c(1.00, 0.40, 0.004),   # min/max = 0.004
    "near_zero_2"     = c(1.00, 0.008),         # min/max = 0.008 (2-df, tiny 2nd)
    "well_conditioned"= c(1.00, 0.90, 0.80),
    "one_df"          = c(1.00)
  )
  ## grid of statistic values spanning the small-p tail
  ## choose S so that p ranges roughly (1e-4 .. 0.5); scale by sum(lam)
  for (nm in names(spectra)) {
    lam <- spectra[[nm]]
    ratio <- min(lam) / max(lam)
    ## quantile-driven S grid: use satterthwaite inverse to hit target p's, plus dense linear grid
    target_p <- c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001, 5e-4, 1e-4)
    cc <- sum(lam^2) / sum(lam); nu <- sum(lam)^2 / sum(lam^2)
    Sgrid <- qchisq(1 - target_p, df = nu) * cc
    for (k in seq_along(Sgrid)) {
      S  <- Sgrid[k]
      pi <- imhof_upper(S, lam)
      ps <- satter_upper(S, lam)
      ek_append(data.frame(spectrum = nm, min_over_max = ratio, k_df = length(lam),
                           S = S, p_imhof = pi, p_satter = ps,
                           abs_delta_p = abs(pi - ps),
                           tail = target_p[k] <= 0.05,
                           row.names = NULL), out_path)
    }
  }
  ## fallback-rule calibration: smallest ratio at which imhof & satter still agree
  ## (max|dp| in the tail <= 0.005) => "use Imhof when min(lam)/max(lam) < c".
  tab <- utils::read.csv(out_path, stringsAsFactors = FALSE)
  tail_tab <- tab[tab$tail, ]
  agg <- aggregate(abs_delta_p ~ spectrum + min_over_max, data = tail_tab, FUN = max)
  agg <- agg[order(agg$min_over_max), ]
  ## c = the largest ratio among spectra whose tail max|dp| still exceeds 0.005
  bad <- agg[agg$abs_delta_p > 0.005, ]
  c_rule <- if (nrow(bad) == 0) NA_real_ else max(bad$min_over_max)
  cat(sprintf("[imhof-vs-satter] tail max|dp| by spectrum (sorted by min/max ratio):\n"))
  print(agg)
  cat(sprintf("[imhof-vs-satter] FALLBACK RULE: use Imhof when min(lam)/max(lam) < c, c = %s\n",
              if (is.na(c_rule)) "NA (satter agrees everywhere tested)" else formatC(c_rule, format = "g", digits = 3)))
  ## persist the rule alongside the scatter (one extra tiny file for the paper)
  ek_append(data.frame(rule = "use_imhof_when_min_over_max_below_c",
                       c_threshold = c_rule,
                       tail_agree_tol = 0.005, row.names = NULL),
            file.path(SIMDIR, "imhof_satterthwaite_rule.csv"))
  invisible(out_path)
}

## ====================================================================================
## MAIN
## ====================================================================================
cat(sprintf("null_calibration_checks.R  SMOKE=%s  N_CAL=%g  B_KS=%d  N_GRID={%s}\n",
            SMOKE, N_CAL, B_KS, paste(N_GRID, collapse = ",")))

## EK_NCORES can be overridden by the caller via EK_NCORES_OVERRIDE (self-test uses 2L to
## avoid oversubscribing while sibling agents run). Otherwise the harness default is used.
.nc_ovr <- suppressWarnings(as.integer(Sys.getenv("EK_NCORES_OVERRIDE", "")))
if (!is.na(.nc_ovr) && .nc_ovr >= 1L) EK_NCORES <- .nc_ovr
cl <- ek_cluster(20260713L)
on.exit(stopCluster(cl), add = TRUE)
## make the DEF helpers + imhof PIT available on workers
clusterExport(cl, c("def_S_lam", ".def_basis_local", "imhof_pit", "imhof_upper"))
invisible(clusterEvalQ(cl, suppressMessages(library(CompQuadForm))))

t0 <- Sys.time()
run_goodness_of_null(cl)          # (a) -> null_ref_lambda.csv + null_ks_table.csv
run_bootstrap_agreement(cl)       # (b) -> null_bootstrap_agreement.csv
run_imhof_satterthwaite()         # (c) -> imhof_satterthwaite.csv (+ _rule.csv)
cat(sprintf("DONE in %.1f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
