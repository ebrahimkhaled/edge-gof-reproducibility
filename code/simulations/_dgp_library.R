## _dgp_library.R  — shared DGPs, the curated battery caller, and diagnostics for the
## EDGE (Phase A) simulations. Generators lifted verbatim-in-spirit from the verified
## sources: Paper_Simu_Inshallah.R (solve_* betas, inv_stukel), _master.R (sqb/sib/scb,
## linkp, gen families), _separation.R (stuk_diag). Curated battery goes through the
## package (dogfoods run.all.gof; verified 31.5 ms/call at n=1000, tests= subsets, no NAs).
## Source AFTER devtools::load_all() of the ebrahim.gof dev tree in each (worker) process.

logit <- function(p) log(p / (1 - p)); lg <- logit
inv_logit <- function(eta) 1 / (1 + exp(-eta))

## --- beta solvers (verified: Paper_Simu_Inshallah.R L204/232/258; _master.R L8-10) ---
sqb <- function(J) solve(rbind(c(1, -1.5, 2.25), c(1, 3, 9), c(1, -3, 9)),
                         c(lg(.05), lg(.95), lg(J)))
sib <- function(I) { e0 <- lg(.1); e1 <- lg(.2); e2 <- lg(.2 + I)
  b0 <- (e0 + e1) / 2; b1 <- (e1 - e0) / 6; b3 <- (e2 - (b0 + 3 * b1)) / 6; c(b0, b1, 3 * b3, b3) }
scb <- function(K) solve(rbind(c(1, -3, -2, 6), c(1, -3, 0, 0), c(1, 3, 0, 0), c(1, 3, 2, 6)),
                         c(lg(.1), lg(.1), lg(.2), lg(.2 + K)))

## --- Stukel inverse link (verified: _master.R L5-7) ---
inv_stukel <- function(ev, a1, a2) { z <- numeric(length(ev)); pos <- ev >= 0
  if (any(pos))  z[pos]  <- if (abs(a1) < 1e-12) ev[pos]  else (-1 + sqrt(pmax(0, 1 + 2 * a1 * ev[pos])))  / a1
  if (any(!pos)) z[!pos] <- if (abs(a2) < 1e-12) ev[!pos] else (-1 + sqrt(pmax(0, 1 + 2 * a2 * ev[!pos]))) / a2
  plogis(z) }

## --- link inverse-cdf dispatcher (verified: _master.R L11-13) ---
linkp <- function(scn, eta) switch(scn,
  logit = plogis(eta), cloglog = 1 - exp(-exp(eta)), loglog = exp(-exp(-eta)), probit = pnorm(eta),
  cauchit = pcauchy(eta), robit_t4 = pt(eta, df = 4), scobit2 = plogis(eta)^2, scobit_half = plogis(eta)^0.5,
  stukel_heavy = inv_stukel(eta, -1, -1), stukel_light = inv_stukel(eta, 1, 1), stukel_asym = inv_stukel(eta, -1, 1),
  stukel_long = inv_stukel(eta, -1, -1), stukel_short = inv_stukel(eta, 1, 1),
  stop("unknown link: ", scn))

EK_LINKS <- c("probit", "cloglog", "stukel_heavy", "stukel_light", "stukel_asym")
EK_ROUGH <- c("osc2", "osc4", "sawtooth", "bump")

## --- ALTERNATIVE (misspecified) DGPs. Returns list(d=data, f=formula, cat=<catvar or NA>) ---
## family in {quad, binint, contint, link, rough}; param = numeric severity or scenario name.
dgp_alt <- function(family, param, n) {
  if (family == "quad") { b <- sqb(param); x <- runif(n, -3, 3)
    list(d = data.frame(x = x, y = rbinom(n, 1, plogis(b[1] + b[2] * x + b[3] * x^2))), f = y ~ x, cat = NA)
  } else if (family == "binint") { b <- sib(param); x <- runif(n, -3, 3); d <- rbinom(n, 1, .5)
    list(d = data.frame(x = x, d = d, y = rbinom(n, 1, plogis(b[1] + b[2] * x + b[3] * d + b[4] * x * d))), f = y ~ x + d, cat = "d")
  } else if (family == "contint") { b <- scb(param); x <- runif(n, -3, 3); z <- rnorm(n)
    list(d = data.frame(x = x, z = z, y = rbinom(n, 1, plogis(b[1] + b[2] * x + b[3] * z + b[4] * x * z))), f = y ~ x + z, cat = NA)
  } else if (family == "link") { x <- runif(n, -3, 3); d <- rbinom(n, 1, .5)
    list(d = data.frame(x = x, d = d, y = rbinom(n, 1, linkp(param, 0.6 * x + 0.5 * d))), f = y ~ x + d, cat = "d")
  } else if (family == "rough") { x <- runif(n, -3, 3)
    eta <- switch(param, osc2 = 0.8 * x + 1.5 * sin(2 * x), osc4 = 0.8 * x + 1.5 * sin(4 * x),
      sawtooth = 0.8 * x + 1.2 * (2 * (x / 1.5 - floor(x / 1.5 + 0.5))),
      bump = 0.3 + 0.7 * x + 1.0 * exp(-x^2 / 2), stop("rough: ", param))
    list(d = data.frame(x = x, y = rbinom(n, 1, plogis(eta))), f = y ~ x, cat = NA)
  } else stop("unknown family: ", family)
}

## --- correctly-specified NULL of each family (generate correct model, fit correct model) ---
dgp_null <- function(family, n) {
  if (family == "quad") { x <- runif(n, -3, 3); list(d = data.frame(x = x, y = rbinom(n, 1, plogis(0.6 * x))), f = y ~ x, cat = NA)
  } else if (family == "binint") { x <- runif(n, -3, 3); d <- rbinom(n, 1, .5)
    list(d = data.frame(x = x, d = d, y = rbinom(n, 1, plogis(0.6 * x + 0.5 * d))), f = y ~ x + d, cat = "d")
  } else if (family == "contint") { x <- runif(n, -3, 3); z <- rnorm(n)
    list(d = data.frame(x = x, z = z, y = rbinom(n, 1, plogis(0.6 * x + 0.5 * z))), f = y ~ x + z, cat = NA)
  } else if (family == "link") { x <- runif(n, -3, 3); d <- rbinom(n, 1, .5)
    list(d = data.frame(x = x, d = d, y = rbinom(n, 1, plogis(0.6 * x + 0.5 * d))), f = y ~ x + d, cat = "d")
  } else if (family == "rough") { x <- runif(n, -3, 3); list(d = data.frame(x = x, y = rbinom(n, 1, plogis(0.8 * x))), f = y ~ x, cat = NA)
  } else stop("unknown family: ", family)
}

## --- benchmark well-spec DGP with p estimated coefficients (Exp 1); k=p-4 nuisance N(0,1) coef 0 ---
gen_bench <- function(n, p = 4L) {
  x1 <- runif(n, -3, 3); x2 <- rnorm(n); d <- rbinom(n, 1, .5)
  df <- data.frame(x1 = x1, x2 = x2, d = d, y = rbinom(n, 1, plogis(-0.2 + 0.6 * x1 + 0.4 * x2 + 0.5 * d)))
  k <- as.integer(p) - 4L
  if (k > 0) { Z <- matrix(rnorm(n * k), n, k); colnames(Z) <- paste0("z", seq_len(k)); df <- cbind(df, Z) }
  f <- as.formula(paste("y ~ x1 + x2 + d", if (k > 0) paste("+", paste0("z", seq_len(k), collapse = " + ")) else ""))
  list(d = df, f = f, cat = "d")
}

## --- sparse DGP for the Stukel-failure experiment (Exp 1C); low intercept => rare events ---
gen_sparse_link <- function(n, intercept = -4.9, beta = 1.0, cov = "chisq4") {
  x <- if (cov == "chisq1") rchisq(n, 1) else rchisq(n, 4)
  x <- as.numeric(scale(x))                    # standardize so intercept controls the event rate
  list(d = data.frame(x = x, y = rbinom(n, 1, plogis(intercept + beta * x))), f = y ~ x, cat = NA)
}

## --- Stukel-link alternative with chosen (a1,a2) for the A(delta)>0 "EDGE loses" case (Exp 2.3.2) ---
gen_unfavorable_asym <- function(n, a1, a2) {
  x <- runif(n, -3, 3); d <- rbinom(n, 1, .5)
  list(d = data.frame(x = x, d = d, y = rbinom(n, 1, inv_stukel(0.6 * x + 0.5 * d, a1, a2))), f = y ~ x + d, cat = "d")
}

## --- Stukel separation diagnostic (verified verbatim: _separation.R L10-16). fit MUST have data= ---
stuk_diag <- function(fit) { e <- predict(fit); d <- fit$data; d$za <- 0.5 * e^2 * (e >= 0); d$zb <- -0.5 * e^2 * (e < 0); sep <- FALSE
  fa <- tryCatch(withCallingHandlers(glm(update(formula(fit), . ~ . + za + zb), data = d, family = binomial()),
    warning = function(w) { if (grepl("did not converge|numerically 0 or 1", conditionMessage(w))) sep <<- TRUE; invokeRestart("muffleWarning") }),
    error = function(er) NULL)
  if (is.null(fa)) return(c(p = NA, fail = 1))
  if (!isTRUE(fa$converged)) sep <- TRUE; pf <- fitted(fa); if (any(pf < 1e-7 | pf > 1 - 1e-7)) sep <- TRUE
  c(p = pchisq(deviance(fit) - deviance(fa), 2, lower.tail = FALSE), fail = as.numeric(sep))
}

## --- alignment functional A(delta)=sum (1-2 pibar_g)(o_g-e_g)/V_g on one large calibration draw ---
A_delta <- function(gen_fn, n_cal = 2e5, G = 10) {
  G0 <- gen_fn(n_cal); fit <- suppressWarnings(glm(G0$f, data = G0$d, family = binomial()))
  ph <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6); y <- G0$d$y
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n_cal / G)), G)
  o <- as.numeric(tapply(y, grp, sum)); e <- as.numeric(tapply(ph, grp, sum))
  pb <- as.numeric(tapply(ph, grp, mean)); V <- as.numeric(tapply(ph * (1 - ph), grp, sum))
  sum((1 - 2 * pb) * (o - e) / V)
}

## --- the curated battery through the package (returns p AND statistic per test, curated order) ---
EK_CURATED <- c("DEF.poly2", "DEF.poly3", "DEF.stukel", "EF", "HL", "HL-equalwidth",
                "Pigeon-Heyse", "Stukel", "Tsiatis", "Xie")
run_curated_full <- function(fit, G = 10) {
  r <- tryCatch(run.all.gof(fit, tests = EK_CURATED, G = G, include_slow = FALSE, install = "no"), error = function(e) NULL)
  if (is.null(r)) { out <- rep(NA_real_, 2 * length(EK_CURATED))
    names(out) <- c(paste0("p.", EK_CURATED), paste0("s.", EK_CURATED)); return(out) }
  idx <- match(EK_CURATED, as.character(r$Test))
  out <- c(as.numeric(r$p_value)[idx], as.numeric(r$Statistic)[idx])
  names(out) <- c(paste0("p.", EK_CURATED), paste0("s.", EK_CURATED)); out
}
