## _block9d_robust.R -- one replicate of block 9d: does robust estimation solve the problem instead?
## Contract: paper_EDGE/theory/PREDECLARATION_block9d_robust.md (sha256 799a85aa...).
##
## One replicate draws ONE data set and fits it TWICE -- maximum likelihood and robust -- so the two
## estimators always see the same data, as section 1 requires. The seed stream is therefore indexed by
## the data configuration (truth x k, 8 of them), and each replicate yields two result rows, one per
## estimator: 8 seed streams, 16 declared cells.
##
## The battery's tests all read one fitted object, fq (bt_fit in _battery_tests.R). For the robust arm
## this file builds the same object from robustbase::glmrob and hands it to the battery's own code, so
## every test statistic is computed by exactly the code blocks 0-9c used. Nothing in
## _battery_tests.R is edited: its bt_fit is replaced for the duration of one call and restored.
##
## What that means, stated here so the analysis cannot forget it: the estimation adjustment in every
## test is the ML one, (X'WX)^{-1}, now evaluated at the robust estimate. That mismatch is exactly what
## claim D9d.1 checks before anything else is read. Stukel's likelihood-ratio refit has no robust
## analogue (a robust fit has no deviance): under the robust arm it is computed as an ML refit against a
## robust base and is NOT a likelihood-ratio test; it carries no claim and is reported separately.

B9D_SEED0 <- 600000000

b9d_configs <- function() {
  C <- expand.grid(k = c(0L, 5L, 10L, 25L), truth = c("logit", "probit"),
                   stringsAsFactors = FALSE)[, c("truth", "k")]
  C$config_id <- seq_len(nrow(C))
  C$seed_base <- B9D_SEED0 + C$config_id * 10000
  C$n <- 1000L
  C
}

## the 16 declared cells: every configuration under each estimator
b9d_cells <- function() {
  C <- b9d_configs()
  X <- rbind(transform(C, estimator = "ML"), transform(C, estimator = "robust"))
  X$cell <- sprintf("%s_k%02d_%s", X$truth, X$k, X$estimator)
  X$type <- "full"; X$ao <- FALSE; X$G_extra <- ""; X$generator <- "b9d"
  X$corruption <- ifelse(X$k == 0L, "clean", "C1"); X$rate <- X$k / X$n
  X[order(X$config_id, X$estimator), ]
}

## The fitted object the battery's tests read, built from a robust fit. Same fields as bt_fit().
b9d_fit_robust <- function(dat) {
  fit <- tryCatch(suppressWarnings(robustbase::glmrob(dat$f, data = dat$d, family = binomial(),
                                                      method = "Mqle")),
                  error = function(e) NULL)
  if (is.null(fit) || !isTRUE(fit$converged)) return(NULL)
  X   <- stats::model.matrix(dat$f, data = dat$d)
  eta <- drop(X %*% stats::coef(fit))
  p   <- stats::plogis(eta)
  ph  <- pmin(pmax(p, 1e-6), 1 - 1e-6)                       # the package's clamp, as bt_fit applies it
  dmu <- stats::binomial()$mu.eta(eta)
  w   <- dmu^2 / (ph * (1 - ph))
  fit$family <- stats::binomial()                            # the tests ask fq$fit$family$link
  ## bt_stukel's likelihood-ratio branch reads fit$rank and fit$deviance, which a robust fit does not
  ## carry. Without them battery_rep() stops -- it does not guard that call -- and the WHOLE robust row
  ## is lost, every test with it (found by block9d_selftest.R). Both have well-defined values: the rank
  ## is the number of coefficients, and the binomial deviance is evaluated at the robust estimate.
  ## Stk.LR then compares an ML augmented fit with a robust base; it is not a likelihood-ratio test and
  ## carries no claim (pre-declaration D9d.2 names Stukel's two SCORE tests, not the refit).
  y <- as.numeric(dat$d$y)
  fit$rank     <- ncol(X)
  fit$deviance <- -2 * sum(y * log(ph) + (1 - y) * log(1 - ph))
  list(fit = fit, y = as.numeric(dat$d$y), eta = eta, ph = ph, p_raw = p, dmu = dmu, X = X,
       A = crossprod(X, w * X), n = nrow(dat$d))
}

## Run the battery's own battery_rep() on a pre-built fit, by standing its bt_fit aside for one call.
b9d_battery_on <- function(dat, cell, fq) {
  env <- environment(battery_rep)
  orig <- get("bt_fit", envir = env)
  assign("bt_fit", function(d) fq, envir = env)
  on.exit(assign("bt_fit", orig, envir = env), add = TRUE)
  battery_rep(dat, cell)
}

b9d_extra <- function(cell) c("k_corrupt", "b.intercept", "b.x", "b.d", "robust_converged",
                              paste0("corrupt_in_top.", names(bt_arms(cell))), "flag.b9d_error")
b9d_names <- function(cell) c(battery_names(cell), b9d_extra(cell))

## One replicate: one data set, two fits, two rows.
b9d_one <- function(rep, config) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  seed <- config$seed_base + rep
  set.seed(seed)
  g   <- b9_gen(config$truth, config$n, if (config$k == 0L) "clean" else "C1", config$k / config$n)
  dat <- list(d = g$d, f = g$f)
  out <- list()
  for (est in c("ML", "robust")) {
    cell <- list(type = "full", ao = FALSE, G_extra = "", n = config$n, generator = "b9d")
    nm   <- b9d_names(cell)
    v <- tryCatch({
      z <- setNames(rep(NA_real_, length(nm)), nm)
      fq <- if (est == "ML") bt_fit(dat) else b9d_fit_robust(dat)
      z["robust_converged"] <- if (est == "robust") as.numeric(!is.null(fq)) else NA_real_
      if (!is.null(fq)) {
        r <- if (est == "ML") battery_rep(dat, cell) else b9d_battery_on(dat, cell, fq)
        z[battery_names(cell)] <- r[battery_names(cell)]
        cf <- stats::coef(fq$fit)
        z["b.intercept"] <- unname(cf[1]); z["b.x"] <- unname(cf[2]); z["b.d"] <- unname(cf[3])
        for (a in names(bt_arms(cell)))
          z[paste0("corrupt_in_top.", a)] <- b9_top_group(fq$ph, g$corrupt, bt_arms(cell)[[a]])
      }
      z["k_corrupt"] <- length(g$corrupt)
      z["flag.b9d_error"] <- 0
      z
    }, error = function(e) {
      z <- setNames(rep(NA_real_, length(nm)), nm); z["flag.b9d_error"] <- 1; z
    })
    out[[est]] <- c(rep = rep, seed = seed, v[nm])
  }
  out
}
