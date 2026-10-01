## _battery_tests.R -- one replicate of the EDGE restructure battery (paper_EDGE/theory/PREDECLARATION_restructure_battery.md,
## sections B and E). battery_one(rep, cell) seeds set.seed(cell$seed_base + rep), draws the cell's data, fits the working
## logistic model once and returns one named numeric vector with every test and flag of E2.
##
## Conventions (E0.3, E0.5, E0.6):
##   * every test, both G arms (G = 10 and G = max(10, round(n/25))), both weighting forms and all four Stukel forms are
##     computed on the same data set;
##   * a test that errors, declines or returns NaN is NA here and is counted as no rejection by the driver;
##     GiViTI warnings are suppressed, only an error gives NA;
##   * HL_w, Pigeon-Heyse, Tsiatis, Xie and Pulkstenis-Robinson run at G = 10 only and not at n >= 10,000.
## No package is loaded here. EDGE and the joint Stukel score follow ebrahim.gof 2.8.0 (def.gof, gof_stukel form "joint");
## HL, HL_F, HL_w, Pigeon-Heyse, Tsiatis, Xie and Pulkstenis-Robinson are line-by-line ports of the package's internal
## wrappers (R/run_all_gof.R, R/ebrahim_farrington_test.R). Block 0 of run_M_battery.R checks both against the installed
## package to 1e-8 and against the stored July p-values to 1e-6.
##
## Real one-off fits (block 6, E8.5): rep 0 is the stored row order and rep k >= 1 the rows in the order sample.int(n)
## drawn after set.seed(20260914 + k), so fits of the same data share their orders. Each row also records the number of
## distinct fitted risks and, per G, the group boundaries that split tied risks.
##
## RNG: workers use L'Ecuyer-CMRG, as ek_cluster() did through clusterSetRNGStream(); set.seed(seed_base + rep) under that
## kind is what makes the July grid cells byte-identical.
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

BT_SIMDIR <- edge_path("code/simulations")
source(file.path(BT_SIMDIR, "_dgp_library.R"))    # dgp_null, dgp_alt, inv_stukel, sqb, sib, scb, gen_sparse_link

bt_rule_G <- function(n) max(10, round(n / 25))
BT_BASES  <- c("poly3", "poly2", "stk", "sym")
BT_EXT_TESTS <- c("EDGE.stk4", "EDGE.poly4", "SlopeLRT", "EDGE.slope2", "HL.ext", "GiViTI.ext", "EDGE.sym3")
BT_RIVAL_NMAX <- 10000
BT_ORDER_SEED <- 20260914                      # E8.5: row order k of a real one-off fit is drawn after set.seed(BT_ORDER_SEED + k)

## ---- links -------------------------------------------------------------------------------------------------
## Stukel's (1988, Sec. 2) exact h-family: logit(mu) = h(eta). a1 acts for eta >= 0, a2 for eta < 0; a = 0 is the logistic.
stukel_h <- function(eta, a1, a2) {
  h <- eta
  pos <- which(eta >= 0); neg <- which(eta < 0)
  if (length(pos) && a1 != 0) {
    e <- eta[pos]
    h[pos] <- if (a1 > 0) (exp(a1 * e) - 1) / a1 else -log(1 - a1 * e) / a1
  }
  if (length(neg) && a2 != 0) {
    e <- -eta[neg]                                 # |eta|
    h[neg] <- if (a2 > 0) -(exp(a2 * e) - 1) / a2 else log(1 - a2 * e) / a2
  }
  h
}
## Hosmer, Hosmer, le Cessie and Lemeshow (1997), Table V shapes
BT_STUKEL_SHAPES <- list(stk_long = c(-1, -1), stk_short = c(1, 1), stk_asym = c(-1, 1))

BT_LINKINV <- list(
  logit = plogis, probit = pnorm, cauchit = pcauchy, t4 = function(e) pt(e, 4),
  loglog = function(e) exp(-exp(-e)), cloglog = function(e) 1 - exp(-exp(e)),
  stk_long  = function(e) plogis(stukel_h(e, -1, -1)),
  stk_short = function(e) plogis(stukel_h(e, 1, 1)),
  stk_asym  = function(e) plogis(stukel_h(e, -1, 1)),
  ## the capped generators of the earlier versions (E0.2), = old stukel_heavy / stukel_light / stukel_asym
  plateau_upper = function(e) inv_stukel(e, -1, -1),
  plateau_lower = function(e) inv_stukel(e, 1, 1),
  plateau_both  = function(e) inv_stukel(e, -1, 1))

## ---- data generators ---------------------------------------------------------------------------------------
## Each returns list(d = data.frame with y, f = fitted formula), or list(y, p) for the external-validation cells.
## Draw order is kept exactly as in the script each generator comes from.
bt_x <- function(n, xdist) switch(xdist, uniform = runif(n, -3, 3),
                                  skewed = 1.5 * as.numeric(scale(rchisq(n, 4))), stop("xdist: ", xdist))

bt_gen_pstar <- function(scenario, n, pstar) {                    # run_A_pstar_giviti.R
  x1 <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); eta_base <- 0.6 * x1 + 0.5 * d
  pr <- switch(scenario, null = plogis(eta_base), cloglog = linkp("cloglog", eta_base),
    quad = { b <- sqb(0.5); plogis(b[1] + b[2] * x1 + b[3] * x1^2) },
    binint = { b <- sib(0.4); plogis(b[1] + b[2] * x1 + b[3] * d + b[4] * x1 * d) }, stop("scenario: ", scenario))
  df <- data.frame(y = rbinom(n, 1, pr), x1 = x1, d = d)
  k <- pstar - 3L
  if (k > 0) { Z <- matrix(rnorm(n * k), n, k); colnames(Z) <- paste0("z", seq_len(k)); df <- cbind(df, Z) }
  f <- stats::as.formula(paste("y ~ x1 + d", if (k > 0) paste("+", paste0("z", seq_len(k), collapse = " + ")) else ""))
  list(d = df, f = f)
}

bt_gen_runG <- function(scn, n) {                                  # run_G_validate_rule.R (non-sparse scenarios)
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5)
  if (scn == "osc4") return(list(d = data.frame(y = rbinom(n, 1, plogis(0.8 * x + 1.5 * sin(4 * x))), x = x), f = y ~ x))
  pr <- switch(scn, null = plogis(0.6 * x + 0.5 * d),
    quad = { b <- sqb(0.5); plogis(b[1] + b[2] * x + b[3] * x^2) },
    binint = { b <- sib(0.4); plogis(b[1] + b[2] * x + b[3] * d + b[4] * x * d) }, stop(scn))
  list(d = data.frame(y = rbinom(n, 1, pr), x = x, d = d), f = y ~ x + d)
}

bt_gen_runK <- function(scn, n) {                                  # run_K_basis_score.R gen + run_K2_ao_basis.R gen2
  x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); eta <- 0.6 * x + 0.5 * d
  if (scn == "osc4") return(list(d = data.frame(y = rbinom(n, 1, plogis(0.8 * x + 1.5 * sin(4 * x))), x = x), f = y ~ x))
  pr <- switch(scn, null = plogis(eta), cloglog = 1 - exp(-exp(eta)), probit = pnorm(eta), loglog = exp(-exp(-eta)),
    quad = { b <- sqb(0.5); plogis(b[1] + b[2] * x + b[3] * x^2) },
    binint = { b <- sib(0.4); plogis(b[1] + b[2] * x + b[3] * d + b[4] * x * d) }, stop(scn))
  list(d = data.frame(y = rbinom(n, 1, pr), x = x, d = d), f = y ~ x + d)
}

bt_gen_census <- function(scn, n) {                                # _master_add.R, the 8 census scenarios
  if (scn == "logx") { x <- runif(n, 0.3, 6); list(d = data.frame(x = x, y = rbinom(n, 1, plogis(-1 + 1.5 * log(x)))), f = y ~ x) }
  else if (scn == "int_binbin") { d1 <- rbinom(n, 1, .5); d2 <- rbinom(n, 1, .5)
    list(d = data.frame(d1 = d1, d2 = d2, y = rbinom(n, 1, plogis(-0.6 + 0.6 * d1 + 1.0 * d2 + 1.4 * d1 * d2))), f = y ~ d1 + d2) }
  else if (scn == "skew") { x <- rchisq(n, 4); xs <- as.numeric(scale(x))
    list(d = data.frame(x = x, y = rbinom(n, 1, plogis(-0.3 + 0.7 * xs + 0.3 * (xs^2 - 1)))), f = y ~ x) }
  else if (scn == "joint") { x <- runif(n, -2.5, 2.5); d <- rbinom(n, 1, .5)
    list(d = data.frame(x = x, y = rbinom(n, 1, plogis(-0.2 + 0.8 * x + 1.0 * d + 0.8 * x * d))), f = y ~ x) }
  else if (scn == "corr") { x1 <- rnorm(n); x2 <- 0.5 * x1 + sqrt(0.75) * rnorm(n)
    list(d = data.frame(x1 = x1, x2 = x2, y = rbinom(n, 1, plogis(-0.2 + 0.7 * x1 + 0.7 * x2 + 0.6 * x1 * x2))), f = y ~ x1 + x2) }
  else if (scn == "omit_x2x3") { x <- runif(n, -2.5, 2.5)
    list(d = data.frame(x = x, y = rbinom(n, 1, plogis(0.3 + 0.7 * x + 0.2 * (x^2 - 2.083) + 0.09 * x^3))), f = y ~ x) }
  else if (scn == "omit_2cov") { x <- runif(n, -2.5, 2.5); z1 <- rnorm(n); z2 <- rnorm(n)
    list(d = data.frame(x = x, y = rbinom(n, 1, plogis(0.2 + 0.6 * x + 0.8 * z1 + 0.8 * z2))), f = y ~ x) }
  else if (scn == "omit_2int") { x <- runif(n, -2.5, 2.5); d <- rbinom(n, 1, .5); z <- rnorm(n)
    list(d = data.frame(x = x, d = d, z = z, y = rbinom(n, 1, plogis(0.2 + 0.6 * x + 0.5 * d + 0.4 * z + 0.7 * x * d + 0.7 * x * z))),
         f = y ~ x + d + z) }
  else stop("census scenario: ", scn)
}

## matched logistic nulls of the census and crossover designs (pre-launch addition, review F1/F2): the covariates are drawn
## exactly as in the alternative, and y from the population limit pi* of the fitted working model under that alternative
## (fractional logistic fit to the true p on a draw of 400,000, seed 20260913, as E0.7), coefficients rounded to 0.01.
BT_NULL_COEF <- list(crossover = c(0, 0.62, 0), logx = c(-1.42, 0.61), int_binbin = c(-0.88, 1.12, 1.54),
                     skew = c(-1.60, 0.30), joint = c(0.18, 1.06), corr = c(-0.05, 0.60, 0.60), omit_x2x3 = c(0.15, 0.93),
                     omit_2cov = c(0.16, 0.48), omit_2int = c(0.12, 0.70, 0.31, 0.22))
bt_gen_design_null <- function(scn, n) {
  b <- BT_NULL_COEF[[scn]]
  if (scn == "crossover") { x <- runif(n, -2.5, 2.5); d <- rbinom(n, 1, 0.5)
    return(list(d = data.frame(x = x, d = d, y = rbinom(n, 1, plogis(b[1] + b[2] * x + b[3] * d))), f = y ~ x + d)) }
  if (scn == "logx") { x <- runif(n, 0.3, 6); return(list(d = data.frame(x = x, y = rbinom(n, 1, plogis(b[1] + b[2] * x))), f = y ~ x)) }
  if (scn == "int_binbin") { d1 <- rbinom(n, 1, .5); d2 <- rbinom(n, 1, .5)
    return(list(d = data.frame(d1 = d1, d2 = d2, y = rbinom(n, 1, plogis(b[1] + b[2] * d1 + b[3] * d2))), f = y ~ d1 + d2)) }
  if (scn == "skew") { x <- rchisq(n, 4); return(list(d = data.frame(x = x, y = rbinom(n, 1, plogis(b[1] + b[2] * x))), f = y ~ x)) }
  if (scn == "joint") { x <- runif(n, -2.5, 2.5); d <- rbinom(n, 1, .5)
    return(list(d = data.frame(x = x, y = rbinom(n, 1, plogis(b[1] + b[2] * x))), f = y ~ x)) }
  if (scn == "corr") { x1 <- rnorm(n); x2 <- 0.5 * x1 + sqrt(0.75) * rnorm(n)
    return(list(d = data.frame(x1 = x1, x2 = x2, y = rbinom(n, 1, plogis(b[1] + b[2] * x1 + b[3] * x2))), f = y ~ x1 + x2)) }
  if (scn == "omit_x2x3") { x <- runif(n, -2.5, 2.5); return(list(d = data.frame(x = x, y = rbinom(n, 1, plogis(b[1] + b[2] * x))), f = y ~ x)) }
  if (scn == "omit_2cov") { x <- runif(n, -2.5, 2.5); z1 <- rnorm(n); z2 <- rnorm(n)
    return(list(d = data.frame(x = x, y = rbinom(n, 1, plogis(b[1] + b[2] * x))), f = y ~ x)) }
  if (scn == "omit_2int") { x <- runif(n, -2.5, 2.5); d <- rbinom(n, 1, .5); z <- rnorm(n)
    return(list(d = data.frame(x = x, d = d, z = z, y = rbinom(n, 1, plogis(b[1] + b[2] * x + b[3] * d + b[4] * z))), f = y ~ x + d + z)) }
  stop("design null: ", scn)
}

## GLOW design of _realdata_power.R, built once per process
bt_glow_env <- new.env()
bt_glow <- function() {
  if (is.null(bt_glow_env$GD)) {
    gl <- aplore3::glow500
    GD <- data.frame(y = as.integer(gl$fracture == "Yes"), age = gl$age, weight = gl$weight,
                     priorfrac = as.integer(gl$priorfrac == "Yes"), momfrac = as.integer(gl$momfrac == "Yes"))
    fit0 <- stats::glm(y ~ age + weight + priorfrac + momfrac, data = GD, family = stats::binomial())
    eta0 <- as.numeric(stats::predict(fit0)); etac <- eta0 - mean(eta0)
    acll <- stats::uniroot(function(a) mean(1 - exp(-exp(a + 1.9 * etac))) - mean(GD$y), c(-15, 15))$root
    bt_glow_env$GD <- GD; bt_glow_env$eta0 <- eta0; bt_glow_env$etac <- etac; bt_glow_env$acll <- acll
    bt_glow_env$z <- as.numeric(scale(GD$age))
  }
  bt_glow_env
}
bt_gen_glow <- function(kind) {
  e <- bt_glow()
  p <- switch(kind, null = plogis(e$eta0), link = 1 - exp(-exp(e$acll + 1.9 * e$etac)),
              nonlin = plogis(e$eta0 + 0.85 * (e$z^2 - 1)), int = plogis(e$eta0 + 1.6 * e$GD$priorfrac * e$z), stop(kind))
  d <- e$GD; d$y <- rbinom(length(p), 1, p)
  list(d = d, f = y ~ age + weight + priorfrac + momfrac)
}

## run_H_overconfidence.R: the frozen model reports a distorted p, y is drawn from the true p
bt_distort <- function(p, kind, s) {
  e <- stats::qlogis(p)
  switch(kind, null = p, temp = stats::plogis(e / s), asym_hi = stats::plogis(ifelse(e >= 0, e / s, e)),
         asym_lo = stats::plogis(ifelse(e < 0, e / s, e)), probit_comp = stats::pnorm(e * 0.588 / s), stop(kind))
}

## one-off real data sets (block 6); no random numbers except the fixed Diabetes-130 split. Read once per process.
bt_real_env <- new.env()
bt_real_data <- function(name) {
  if (is.null(bt_real_env[[name]])) bt_real_env[[name]] <- bt_real_data0(name)
  bt_real_env[[name]]
}
bt_real_data0 <- function(name) {
  if (name == "beetle_cloglog") { z <- bt_real_data0("beetle_logit"); z$link <- "cloglog"; return(z) }   # tab:concord refit row
  if (name == "beetle_logit") {
    be <- data.frame(dose = c(1.6907, 1.7242, 1.7552, 1.7842, 1.8113, 1.8369, 1.8610, 1.8839),
                     n = c(59, 60, 62, 56, 63, 59, 62, 60), k = c(6, 13, 18, 28, 52, 53, 61, 60))
    d <- do.call(rbind, Map(function(ds, n, k) data.frame(dose = ds, y = c(rep(1, k), rep(0, n - k))), be$dose, be$n, be$k))
    return(list(d = d, f = y ~ dose))
  }
  if (name %in% c("lbw_additive", "lbw_fix")) {
    bw <- MASS::birthwt
    d <- data.frame(y = bw$low, age = bw$age, lwt = bw$lwt, race = factor(bw$race), smoke = bw$smoke)
    if (name == "lbw_additive") return(list(d = d, f = y ~ age + lwt + race + smoke))
    d$al <- d$age * d$lwt; d$sl <- d$smoke * d$lwt
    return(list(d = d, f = y ~ age + lwt + race + smoke + al + sl))
  }
  if (name %in% c("vaso_linear", "vaso_log")) {
    e <- new.env(); utils::data("vaso", package = "catdata", envir = e)
    d <- data.frame(y = e$vaso$vaso - min(e$vaso$vaso), vol = e$vaso$vol, rate = e$vaso$rate); d$y <- ifelse(d$y > 0, 1, 0)
    return(list(d = d, f = if (name == "vaso_linear") y ~ vol + rate else y ~ log(vol) + log(rate)))
  }
  if (name == "uis_model19") {
    d <- readRDS(file.path(BT_SIMDIR, "uis_data.rds")); d$y <- as.integer(d$DFREE == "no"); d$NDRGFP1 <- 10 / (d$NDRGTX + 1)
    return(list(d = d, f = y ~ AGE + NDRGTX + IVHX + RACE + TREAT + SITE + AGE:NDRGFP1 + RACE:SITE))
  }
  if (name %in% c("diabetes_dev", "diabetes_val")) {                # run_I_bigdata.R
    D <- utils::read.csv(file.path(BT_SIMDIR, "..", "data_large", "diabetic_data.csv"), stringsAsFactors = FALSE,
                         na.strings = c("?", ""))
    D <- D[D$gender %in% c("Female", "Male"), ]
    D$y <- as.integer(D$readmitted == "<30")
    age_mid <- c("[0-10)" = 5, "[10-20)" = 15, "[20-30)" = 25, "[30-40)" = 35, "[40-50)" = 45, "[50-60)" = 55,
                 "[60-70)" = 65, "[70-80)" = 75, "[80-90)" = 85, "[90-100)" = 95)
    D$age_n <- unname(age_mid[D$age]); D$female <- as.integer(D$gender == "Female")
    D$insulin <- factor(D$insulin, levels = c("No", "Steady", "Up", "Down"))
    D$change <- as.integer(D$change == "Ch"); D$diabetesMed <- as.integer(D$diabetesMed == "Yes")
    f <- y ~ age_n + female + time_in_hospital + num_lab_procedures + num_procedures + num_medications +
      number_outpatient + number_emergency + number_inpatient + number_diagnoses + insulin + change + diabetesMed
    D <- D[stats::complete.cases(D[, all.vars(f)]), ]
    set.seed(20260917); dev <- sample(nrow(D), floor(nrow(D) / 2))
    if (name == "diabetes_dev") return(list(d = D[dev, ], f = f))
    fit <- stats::glm(f, data = D[dev, ], family = stats::binomial())
    Dv <- D[-dev, ]
    return(list(y = Dv$y, p = as.numeric(stats::predict(fit, newdata = Dv, type = "response"))))
  }
  stop("real data set: ", name)
}
## a real data set in another row order (E8.5): the rows of the data frame, or y and p together in external mode
bt_nrow <- function(dat) if (is.null(dat$d)) length(dat$y) else nrow(dat$d)
bt_reorder <- function(dat, o) {
  if (is.null(dat$d)) { dat$y <- dat$y[o]; dat$p <- dat$p[o] }
  else { dat$d <- dat$d[o, , drop = FALSE]; rownames(dat$d) <- NULL }
  dat
}

bt_data <- function(cell) {
  n <- cell$n
  switch(cell$generator,
    dgp_null = dgp_null(cell$family, n),
    dgp_alt = dgp_alt(cell$family, if (cell$family %in% c("quad", "binint", "contint")) as.numeric(cell$param) else cell$param, n),
    design_null = bt_gen_design_null(cell$scen, n),
    crossover = {                                                  # grid_edge_loses.R gen_crossover
      x <- runif(n, -2.5, 2.5); d <- rbinom(n, 1, 0.5)
      list(d = data.frame(x = x, d = d, y = rbinom(n, 1, plogis(0.9 * x - 0.5 * x * d))), f = y ~ x + d) },
    omit_2cov = {                                                  # grid_proj_power.R gen_cell("omit_2cov")
      x <- runif(n, -2.5, 2.5); z1 <- rnorm(n); z2 <- rnorm(n)
      list(d = data.frame(x = x, y = rbinom(n, 1, plogis(0.2 + 0.6 * x + 0.8 * z1 + 0.8 * z2))), f = y ~ x) },
    design = {                                                     # E1: eta = c0 + s(0.6x + 0.5d), fitted y ~ x + d
      x <- bt_x(n, cell$xdist); d <- rbinom(n, 1, 0.5); eta <- cell$c0 + cell$s * (0.6 * x + 0.5 * d)
      list(d = data.frame(x = x, d = d, y = rbinom(n, 1, BT_LINKINV[[cell$link]](eta))), f = y ~ x + d) },
    sparse = {                                                     # gen_sparse_link covariates: x = scale(chi-square(4))
      x <- rchisq(n, 4); x <- as.numeric(scale(x))
      list(d = data.frame(x = x, y = rbinom(n, 1, BT_LINKINV[[cell$link]](cell$intercept + cell$slope * x))), f = y ~ x) },
    pstar = bt_gen_pstar(cell$scen, n, as.integer(cell$pstar)),
    runG = bt_gen_runG(cell$scen, n),
    runK = bt_gen_runK(cell$scen, n),
    census = bt_gen_census(cell$scen, n),
    glow = bt_gen_glow(cell$scen),
    external = { x <- runif(n, -3, 3); d <- rbinom(n, 1, .5); p <- stats::plogis(0.6 * x + 0.5 * d); y <- rbinom(n, 1, p)
      list(y = y, p = bt_distort(p, cell$kind, cell$dist_s)) },
    real = bt_real_data(cell$dataset),
    stop("unknown generator: ", cell$generator))
}

## ---- the fitted model, once per replicate ---------------------------------------------------------------------
bt_fit <- function(dat) {
  fam <- if (is.null(dat[["link"]])) stats::binomial() else stats::binomial(link = dat[["link"]])   # a link other than logit: real_beetle_cloglog only
  fit <- tryCatch(suppressWarnings(stats::glm(dat$f, data = dat$d, family = fam)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  eta <- as.numeric(fit$linear.predictors)
  ph  <- pmin(pmax(as.numeric(stats::fitted(fit)), 1e-6), 1 - 1e-6)     # the package's clamp
  dmu <- fit$family$mu.eta(eta)
  X   <- stats::model.matrix(fit)
  w   <- dmu^2 / (ph * (1 - ph))
  list(fit = fit, y = as.numeric(fit$y), eta = eta, ph = ph, p_raw = as.numeric(stats::fitted(fit)), dmu = dmu, X = X,
       A = crossprod(X, w * X), n = length(fit$y))
}

## ---- EDGE: one grouping per G, shared by every basis and both forms -----------------------------------------------
## The group sums, Omega and both statistics use def.gof's own arithmetic, operation for operation (vapply sums, explicit
## Omega = I - U (X'WX)^-1 U'). An algebraically equal shortcut (rowsum, Z'Omega Z without Omega) differed from the package by
## up to 1.5e-7 on ill-conditioned sparse samples, beyond the 1e-8 of E2.
bt_groups <- function(fq, G) {
  n <- fq$n; ph <- fq$ph; y <- fq$y; dmu <- fq$dmu; X <- fq$X
  if (G > n) return(NULL)
  V <- ph * (1 - ph)
  grp  <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)     # def.gof's equal-frequency grouping
  idx  <- split(seq_len(n), grp)
  og   <- vapply(idx, function(I) sum(y[I]),   numeric(1))
  eg   <- vapply(idx, function(I) sum(ph[I]),  numeric(1))
  Vg   <- vapply(idx, function(I) sum(V[I]),   numeric(1))
  pbar <- vapply(idx, function(I) mean(ph[I]), numeric(1))
  r    <- (og - eg) / sqrt(Vg)
  U    <- t(vapply(idx, function(I) colSums(dmu[I] * X[I, , drop = FALSE]), numeric(ncol(X)))) / sqrt(Vg)
  Omega <- diag(length(idx)) - U %*% solve(fq$A) %*% t(U)
  list(G = G, grp = grp, ng = lengths(idx), og = og, eg = eg, Vg = Vg, pbar = pbar, r = r, U = U, Omega = Omega)
}
bt_zoz <- function(gs, A, Z) t(Z) %*% gs$Omega %*% Z

bt_basis <- function(pbar, basis) {
  if (basis %in% c("poly2", "poly3")) {
    deg <- if (basis == "poly2") 2L else 3L
    if (length(unique(round(pbar, 8))) < deg + 1) return(NULL)          # def.gof stops here
    Z <- as.matrix(stats::poly(pbar, deg))
  } else {
    e <- stats::qlogis(pbar)
    Z <- switch(basis, stk = cbind(e, e^2 * (e >= 0), -e^2 * (e < 0)), sym = cbind(e * abs(e)),
                ao = cbind(1 - log1p(exp(e)) / plogis(e)), stop("basis: ", basis))   # ao = run K's ao_col
  }
  ## def.gof's column rule (E8.3): a column with sum |z| <= 1e-8 is dropped, as in 2.7.0, and the kept columns are scaled
  ## to unit length before any solve. A Stukel half-column that reaches one group whose mean risk is a hair above or below
  ## 0.5 is tiny but kept. Both statistics and the eigenvalues are scale-free; solve() is not.
  Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]
  if (!ncol(Z)) return(NULL)
  Z / rep(sqrt(colSums(Z^2)), each = nrow(Z))
}

## unit form: S = r'P_Z r against sum_j lambda_j chi2_1, Satterthwaite (def.gof's default method)
bt_edge_unit <- function(gs, A, Z) {
  ZtZ <- crossprod(Z); Zr <- crossprod(Z, gs$r)
  S <- as.numeric(t(Zr) %*% solve(ZtZ) %*% Zr)
  lam <- Re(eigen(solve(ZtZ) %*% (t(Z) %*% gs$Omega %*% Z), only.values = TRUE)$values)
  lam <- lam[lam > 1e-9]
  if (!length(lam)) return(list(p = NA_real_, S = S, lam = numeric(0)))
  cc <- sum(lam^2) / sum(lam); nu <- sum(lam)^2 / sum(lam^2)
  list(p = stats::pchisq(S / cc, nu, lower.tail = FALSE), S = S, lam = sort(lam, decreasing = TRUE))
}

## score form: columns times sqrt(V_g), u'I^-1 u on the rank of I read after scaling I to a correlation matrix.
## A column whose post-fit information is below 1e-10 of its unadjusted information (Zs'Zs) is one the model already spans
## (for example every column of a sample whose fitted logit is constant); it is left out, and no column left means no p-value.
bt_edge_score <- function(gs, A, Z) {
  Zs <- Z * sqrt(gs$Vg)
  u <- drop(crossprod(Zs, gs$r))
  I <- crossprod(Zs, gs$Omega %*% Zs)
  d <- sqrt(pmax(diag(I), 0)); ok <- d > 0 & diag(I) > 1e-10 * colSums(Zs^2)
  guard <- any(!ok & is.finite(diag(I)))                            # a column left out by this rule, zero or negative information included
  if (!any(ok)) return(list(p = NA_real_, S = NA_real_, k = 0L, guard = guard))
  R <- I[ok, ok, drop = FALSE] / outer(d[ok], d[ok])
  ev <- eigen((R + t(R)) / 2, symmetric = TRUE)
  pos <- ev$values > 1e-8; k <- sum(pos)
  if (!k) return(list(p = NA_real_, S = NA_real_, k = 0L, guard = guard))
  S <- sum(drop(crossprod(ev$vectors[, pos, drop = FALSE], u[ok] / d[ok]))^2 / ev$values[pos])
  list(p = stats::pchisq(S, k, lower.tail = FALSE), S = S, k = k, guard = guard)
}

## ---- Hosmer-Lemeshow family (ports of the package wrappers) ---------------------------------------------------------
bt_hl_stat <- function(y, ph, grp) {                               # .gof_hl_stat
  idx <- split(seq_along(y), grp)
  O  <- vapply(idx, function(I) sum(y[I]), numeric(1))
  E  <- vapply(idx, function(I) sum(ph[I]), numeric(1))
  ng <- vapply(idx, length, numeric(1))
  keep <- ng > 0 & E > 1e-8 & E < ng - 1e-8
  hl <- sum((O[keep] - E[keep])^2 / (E[keep] * (1 - E[keep] / ng[keep])))
  df <- sum(keep) - 2
  if (df < 1) NA_real_ else stats::pchisq(hl, df, lower.tail = FALSE)
}
bt_hlw <- function(y, ph, G) {                                     # gof_hlw
  br <- seq(0, 1, length.out = G + 1); br[1] <- -Inf; br[length(br)] <- Inf
  bt_hl_stat(y, ph, cut(ph, breaks = br, labels = FALSE))
}
bt_hlf <- function(y, ph, G) {                                     # ef.gof -> .ebrahim_farrington_grouped, chisq reference
  n <- length(y)
  if (G < 2 || G > n) return(NA_real_)
  o <- order(ph); ys <- y[o]; ps <- ph[o]
  ga <- rep(1:G, each = floor(n / G)); if (n %% G > 0) ga <- c(ga, rep(G, n %% G))
  yg <- as.numeric(tapply(ys, ga, sum)); mg <- as.numeric(tapply(ys, ga, length)); pg <- as.numeric(tapply(ps, ga, mean))
  pb <- pmax(pmin(pg, 0.999999), 0.000001)
  v <- mg * pb * (1 - pb)
  Tst <- sum((yg - mg * pb)^2 / v) - sum(((1 - 2 * pb) / v) * (yg - mg * pb))
  stats::pchisq(Tst, df = G - 2, lower.tail = FALSE)
}
bt_ph <- function(y, ph, g) {                                      # gof_ph_test (Pigeon-Heyse J2)
  br  <- stats::quantile(ph, probs = seq(0, 1, length.out = g + 1))
  grp <- tryCatch(cut(ph, breaks = br, labels = FALSE, include.lowest = TRUE), error = function(e) NULL)
  if (is.null(grp) || length(unique(grp[!is.na(grp)])) < 2) return(NA_real_)
  idx  <- split(seq_along(y), grp)
  Ok   <- vapply(idx, function(I) sum(y[I]), numeric(1))
  nk   <- vapply(idx, length, numeric(1))
  pbar <- vapply(idx, function(I) mean(ph[I]), numeric(1))
  Vk   <- nk * pbar * (1 - pbar)
  phik <- vapply(idx, function(I) sum(ph[I] * (1 - ph[I])), numeric(1)) / Vk
  ok   <- is.finite(Vk) & Vk > 0 & is.finite(phik) & phik > 0
  J2   <- sum(((Ok[ok] - nk[ok] * pbar[ok])^2 / Vk[ok]) / phik[ok])
  stats::pchisq(J2, length(idx) - 1, lower.tail = FALSE)
}

## ---- covariate-space rivals (ports) --------------------------------------------------------------------------------
bt_ginv <- function(M, tol = sqrt(.Machine$double.eps)) {           # .gof_ginv
  s <- svd(M); pos <- s$d > max(tol * s$d[1], 0)
  if (!any(pos)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}
bt_kmeans <- function(mat, centers, nstart = 1L) {                 # .gof_kmeans: seed 123, caller's RNG state restored
  has <- exists(".Random.seed", envir = .GlobalEnv)
  old <- if (has) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(123)
  cl <- stats::kmeans(mat, centers = centers, nstart = nstart)$cluster
  if (has) assign(".Random.seed", old, envir = .GlobalEnv)
  cl
}
bt_tsiatis <- function(fq, G = 10) {                               # gof_tsiatis
  X <- fq$X; ph <- fq$ph; y <- fq$y
  cov_mat <- X[, -1, drop = FALSE]
  if (ncol(cov_mat) < 1) return(NA_real_)
  cl <- tryCatch(bt_kmeans(cov_mat, G), error = function(e) NULL)
  if (is.null(cl)) return(NA_real_)
  Xc <- stats::model.matrix(~ factor(cl))[, -1, drop = FALSE]
  if (ncol(Xc) < 1) return(NA_real_)
  W <- ph * (1 - ph); U <- colSums(Xc * (y - ph))
  V11 <- crossprod(X, W * X); V12 <- crossprod(X, W * Xc); V22 <- crossprod(Xc, W * Xc)
  V <- V22 - t(V12) %*% bt_ginv(V11) %*% V12
  Tstat <- as.numeric(t(U) %*% bt_ginv(V) %*% U)
  rankV <- sum(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 1e-8)
  stats::pchisq(Tstat, rankV, lower.tail = FALSE)
}
bt_xie <- function(fq) {                                           # gof_xie
  X <- fq$X; ph <- fq$ph; y <- fq$y
  k <- ncol(X) - 1; cov_mat <- X[, -1, drop = FALSE]
  if (ncol(cov_mat) < 1) return(NA_real_)
  G <- if (k < 5) 10 else k + 5
  cl <- tryCatch(bt_kmeans(cov_mat, G, nstart = 25L), error = function(e) NULL)
  if (is.null(cl)) return(NA_real_)
  stat <- 0
  for (I in split(seq_along(y), cl)) {
    ng <- length(I); pbar <- mean(ph[I])
    if (pbar > 0 && pbar < 1) stat <- stat + (sum(y[I]) - ng * pbar)^2 / (ng * pbar * (1 - pbar))
  }
  df <- G - k / 2 - 1
  if (df <= 0) NA_real_ else stats::pchisq(stat, df, lower.tail = FALSE)
}
bt_pr <- function(fq) {                                            # gof_pr (categorical covariates found automatically)
  data <- tryCatch(stats::model.frame(fq$fit), error = function(e) NULL)
  if (is.null(data)) return(NA_real_)
  cand <- names(data)[-1]
  maxlev <- getOption("ebrahim.gof.pr.maxlev", 6)
  cat_vars <- cand[vapply(cand, function(v) {
    col <- data[[v]]
    is.factor(col) || is.character(col) || is.logical(col) || (is.numeric(col) && length(unique(col)) <= maxlev)
  }, logical(1))]
  if (length(cat_vars) == 0) return(NA_real_)
  y <- fq$y; ph <- fq$ph
  patt <- do.call(paste, c(lapply(cat_vars, function(v) as.character(data[[v]])), sep = "_"))
  M <- length(unique(patt))
  lev <- character(length(y))
  for (pp in unique(patt)) { ix <- which(patt == pp); lev[ix] <- ifelse(ph[ix] <= stats::median(ph[ix]), "low", "high") }
  idx <- split(seq_along(y), paste(patt, lev, sep = "::"))
  os <- vapply(idx, function(I) sum(y[I] == 1), numeric(1)); of <- vapply(idx, function(I) sum(y[I] == 0), numeric(1))
  es <- vapply(idx, function(I) sum(ph[I]), numeric(1));     ef <- vapply(idx, function(I) sum(1 - ph[I]), numeric(1))
  keep <- es > 0 & ef > 0
  chisq <- sum((os[keep] - es[keep])^2 / es[keep] + (of[keep] - ef[keep])^2 / ef[keep])
  df <- 2 * M - length(cat_vars) - 2
  if (df <= 0) NA_real_ else stats::pchisq(chisq, df, lower.tail = FALSE)
}

## ---- Stukel: joint score, likelihood-ratio refit, one-parameter symmetric score, marginal sum ------------------------
## joint = gof_stukel(form = "joint"): u'I^-1 u on the non-zero half-columns (1 df when one half is identically zero).
## As in the package, W, u and I use the raw fitted risks (not clamped to [1e-6, 1 - 1e-6]) and aliased model columns
## are left out; the clamped risks only decide the side of 0.5, which the clamp cannot change.
bt_stukel_joint_stat <- function(fq) {
  y <- fq$y; ph <- fq$ph; eta <- fq$eta; pr <- fq$p_raw; W <- pr * (1 - pr)
  X <- fq$X[, !is.na(stats::coef(fq$fit)), drop = FALSE]
  Z <- cbind(za = 0.5 * eta^2 * (ph >= 0.5), zb = -0.5 * eta^2 * (ph < 0.5))
  keep <- colSums(Z != 0) > 0
  if (!any(keep)) return(list(chi = NA_real_, df = 0L, half0 = NA_real_, Z = Z, keep = keep, guard = FALSE))
  Zk <- Z[, keep, drop = FALSE]
  u <- colSums(Zk * (y - pr))
  df <- ncol(Zk); guard <- FALSE
  chi <- tryCatch({
    ZWX <- crossprod(Zk, W * X)
    I <- crossprod(Zk, W * Zk) - ZWX %*% solve(crossprod(X, W * X), t(ZWX))
    ## a half-column the model already spans (post-fit information below 1e-10 of its unadjusted information) is left out
    ok <- diag(I) > 1e-10 * colSums(W * Zk^2)
    guard <- any(!ok); df <- sum(ok)
    if (!df) NA_real_ else as.numeric(crossprod(u[ok], solve(I[ok, ok, drop = FALSE], u[ok])))
  }, error = function(e) NA_real_)
  if (!is.finite(chi) || chi < 0) chi <- NA_real_
  list(chi = chi, df = df, half0 = as.numeric(!all(keep)), Z = Z, keep = keep, guard = guard)
}

bt_stukel <- function(fq) {
  out <- c(Stk.joint = NA_real_, Stk.LR = NA_real_, Stk.sym1 = NA_real_, Stk.marg = NA_real_,
           flag.stk_half0 = NA_real_, flag.stk_diag = NA_real_, flag.info_guard = 0)
  js <- bt_stukel_joint_stat(fq)
  out["flag.stk_half0"] <- js$half0
  if (isTRUE(js$guard)) out["flag.info_guard"] <- 1
  if (is.finite(js$chi)) out["Stk.joint"] <- stats::pchisq(js$chi, js$df, lower.tail = FALSE)

  ## marginal sum: the 2.6.0 / LogisticDx statistic, statmod::glm.scoretest on both half-columns (NaN when one is zero)
  zm <- tryCatch(abs(statmod::glm.scoretest(fq$fit, js$Z)), error = function(e) NULL)
  if (!is.null(zm)) { pm <- stats::pchisq(sum(zm^2), 2, lower.tail = FALSE); if (is.finite(pm)) out["Stk.marg"] <- pm }

  ## one-parameter symmetric score (alpha1 = alpha2), ungrouped column eta|eta|
  ## raw fitted risks and the non-aliased columns, as in the joint score (E10.1)
  y <- fq$y; pr <- fq$p_raw; W <- pr * (1 - pr)
  X <- fq$X[, !is.na(stats::coef(fq$fit)), drop = FALSE]
  zs <- fq$eta * abs(fq$eta)
  if (any(zs != 0)) {
    I1 <- tryCatch({ zWX <- crossprod(X, W * zs); sum(W * zs^2) - as.numeric(crossprod(zWX, solve(crossprod(X, W * X), zWX))) },
                   error = function(e) NA_real_)
    thr <- 1e-10 * sum(W * zs^2)
    if (is.finite(I1) && !(I1 > thr)) out["flag.info_guard"] <- 1
    if (is.finite(I1) && I1 > thr) out["Stk.sym1"] <- stats::pchisq(sum(zs * (y - pr))^2 / I1, 1, lower.tail = FALSE)
  }

  ## likelihood-ratio refit (gof_stukel form "lr": glm.fit, default control); the same refit gives stuk_diag's flag
  ## (non-convergence, a "did not converge" / "numerically 0 or 1" warning, or a fitted value within 1e-7 of 0 or 1)
  sep <- FALSE
  f1 <- tryCatch(withCallingHandlers(
      stats::glm.fit(cbind(X, js$Z[, js$keep, drop = FALSE]), y, family = stats::binomial()),
      warning = function(w) {
        if (grepl("did not converge|numerically 0 or 1", conditionMessage(w))) sep <<- TRUE
        invokeRestart("muffleWarning")
      }), error = function(e) NULL)
  if (is.null(f1)) {
    out["flag.stk_diag"] <- 1
  } else {
    pf <- f1$fitted.values
    out["flag.stk_diag"] <- as.numeric(sep || !isTRUE(f1$converged) || any(pf < 1e-7 | pf > 1 - 1e-7))
    if (isTRUE(f1$converged) && any(js$keep)) {
      k <- f1$rank - fq$fit$rank
      if (k >= 1) out["Stk.LR"] <- stats::pchisq(max(fq$fit$deviance - f1$deviance, 0), k, lower.tail = FALSE)
    }
  }
  out
}

## ---- GiViTI and the fixed-degree cubic calibration LR --------------------------------------------------------------------
bt_giviti <- function(y, p, thres) tryCatch(suppressWarnings(
  as.numeric(givitiR::givitiCalibrationTest(y, p, devel = "internal", thres = thres)$p.value)), error = function(e) NA_real_)

bt_cubic <- function(y, eta) tryCatch({
  d1 <- suppressWarnings(stats::glm.fit(cbind(1, eta), y, family = stats::binomial()))$deviance
  d3 <- suppressWarnings(stats::glm.fit(cbind(1, eta, eta^2, eta^3), y, family = stats::binomial()))$deviance
  stats::pchisq(max(d1 - d3, 0), 2, lower.tail = FALSE) }, error = function(e) NA_real_)

## ---- external-validation mode (run H): frozen predictions, Omega = I ----------------------------------------------------
bt_ext_grp <- function(p, G) {
  br <- stats::quantile(p, probs = seq(0, 1, length.out = G + 1), type = 1); br[1] <- -Inf; br[length(br)] <- Inf
  as.integer(cut(p, breaks = unique(br), include.lowest = TRUE))
}
bt_ext_resid <- function(y, p, G) {
  g <- bt_ext_grp(p, G); o <- tapply(y, g, sum); e <- tapply(p, g, sum); v <- tapply(p * (1 - p), g, sum)
  list(r = as.numeric((o - e) / sqrt(v)), eta = as.numeric(stats::qlogis(tapply(p, g, mean))), K = max(g))
}
bt_ext_edge_s <- function(z, basis) {                             # statistic, df and p-value
  eta <- z$eta
  Z <- switch(basis, slope = cbind(1, eta), poly3 = cbind(1, eta, eta^2, eta^3),
              stukel = cbind(1, eta, eta^2 * (eta >= 0), -eta^2 * (eta < 0)),
              sym = cbind(1, eta, eta * abs(eta)))                   # EDGE-sym external column (new, E3 block 5)
  Z <- Z[, apply(Z, 2, function(c) stats::sd(c) > 1e-10 | all(c == 1)), drop = FALSE]
  Q <- qr(Z); if (Q$rank < ncol(Z)) Z <- Z[, Q$pivot[seq_len(Q$rank)], drop = FALSE]
  S <- sum(qr.fitted(qr(Z), z$r)^2)
  c(S = S, df = ncol(Z), p = stats::pchisq(S, df = ncol(Z), lower.tail = FALSE))
}
bt_ext_edge <- function(z, basis) bt_ext_edge_s(z, basis)[["p"]]
BT_EXT_GROUPED <- c(EDGE.stk4 = "stukel", EDGE.poly4 = "poly3", EDGE.slope2 = "slope", EDGE.sym3 = "sym")
bt_ext_G <- function(cell) {                                       # extra G of an external cell (tab:diabetes validation half)
  if (is.null(cell$G_extra) || is.na(cell$G_extra) || !nzchar(cell$G_extra)) integer(0)
  else as.integer(strsplit(as.character(cell$G_extra), ",")[[1]])
}
bt_ext_slope <- function(y, p) {
  lp <- stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  f1 <- tryCatch(stats::glm(y ~ lp, family = stats::binomial()), error = function(e) NULL)
  if (is.null(f1)) return(NA_real_)
  ll0 <- sum(y * log(p) + (1 - y) * log(1 - p))
  stats::pchisq(2 * (as.numeric(stats::logLik(f1)) - ll0), df = 2, lower.tail = FALSE)
}

battery_rep_external <- function(dat, cell) {
  y <- dat$y; p <- dat$p; n <- length(y); G <- max(6L, as.integer(round(n / 25)))
  out <- setNames(rep(NA_real_, length(BT_EXT_TESTS)), BT_EXT_TESTS)
  z <- tryCatch(bt_ext_resid(y, p, G), error = function(e) NULL)
  if (!is.null(z)) {
    out["EDGE.stk4"]   <- tryCatch(bt_ext_edge(z, "stukel"), error = function(e) NA_real_)
    out["EDGE.poly4"]  <- tryCatch(bt_ext_edge(z, "poly3"), error = function(e) NA_real_)
    out["EDGE.slope2"] <- tryCatch(bt_ext_edge(z, "slope"), error = function(e) NA_real_)
    out["EDGE.sym3"]   <- tryCatch(bt_ext_edge(z, "sym"), error = function(e) NA_real_)
    out["HL.ext"]      <- stats::pchisq(sum(z$r^2), df = z$K, lower.tail = FALSE)
  }
  out["SlopeLRT"] <- tryCatch(bt_ext_slope(y, p), error = function(e) NA_real_)
  out["GiViTI.ext"] <- tryCatch(suppressWarnings(
    as.numeric(givitiR::givitiCalibrationTest(y, p, devel = "external")$p.value)), error = function(e) NA_real_)
  ext <- c()
  for (Gx in bt_ext_G(cell)) {                                     # the grouped tests again at each extra G, with statistics
    zx <- tryCatch(bt_ext_resid(y, p, Gx), error = function(e) NULL)
    for (t in names(BT_EXT_GROUPED)) {
      s <- if (is.null(zx)) c(S = NA_real_, df = NA_real_, p = NA_real_) else
        tryCatch(bt_ext_edge_s(zx, BT_EXT_GROUPED[[t]]), error = function(e) c(S = NA_real_, df = NA_real_, p = NA_real_))
      ext[paste0(t, ".G", Gx)] <- s[["p"]]; ext[paste0("stat.", t, ".G", Gx)] <- s[["S"]]
    }
    hs <- if (is.null(zx)) NA_real_ else sum(zx$r^2)
    ext[paste0("HL.ext.G", Gx)] <- if (is.null(zx)) NA_real_ else stats::pchisq(hs, df = zx$K, lower.tail = FALSE)
    ext[paste0("stat.HL.ext.G", Gx)] <- hs
  }
  if (identical(cell$generator, "real")) {                         # E8.5: ties in the frozen risks (cut() never splits them)
    ext["info.distinct_risks"] <- length(unique(p))
    gl <- c(Gdef = G, setNames(bt_ext_G(cell), paste0("G", bt_ext_G(cell))))
    for (a in unique(names(gl)))
      ext[paste0("info.tied_boundaries.", a)] <- tryCatch(bt_tied_boundaries(p, bt_ext_grp(p, gl[[a]])), error = function(e) NA_real_)
  }
  c(n = n, events = sum(y), out, ext, flag.rep_error = 0)
}

## ---- ties in the fitted risks (real one-off fits, E8.5) ---------------------------------------------------------------
## the number of boundaries between consecutive groups with the same fitted risk on both sides (groups in risk order)
bt_tied_boundaries <- function(ph, grp) {
  s <- split(ph, grp)
  k <- length(s)
  if (k < 2) return(0)
  mx <- vapply(s, max, numeric(1)); mn <- vapply(s, min, numeric(1))
  sum(mx[-k] == mn[-1])
}
## the same for bt_hlf's blocks of the sorted risks
bt_tied_boundaries_hlf <- function(ph, G) {
  n <- length(ph)
  if (G < 2 || G > n) return(NA_real_)
  ps <- sort(ph); ga <- rep(1:G, each = floor(n / G)); if (n %% G > 0) ga <- c(ga, rep(G, n %% G))
  b <- which(diff(ga) != 0)
  sum(ps[b] == ps[b + 1])
}
bt_order_names <- function(cell) {
  if (!identical(cell$generator, "real")) return(character(0))
  if (cell$type == "external")
    return(c("info.distinct_risks", paste0("info.tied_boundaries.", unique(c("Gdef", paste0("G", bt_ext_G(cell)))))))
  arms <- names(bt_arms(cell))
  c("info.distinct_risks", paste0("info.tied_boundaries.", arms), paste0("info.tied_boundaries_hlf.", arms))
}

## ---- names of the output vector -------------------------------------------------------------------------------------------
bt_arms <- function(cell) {
  ex <- if (is.null(cell$G_extra) || is.na(cell$G_extra) || !nzchar(cell$G_extra)) integer(0)
        else as.integer(strsplit(as.character(cell$G_extra), ",")[[1]])
  a <- c(G10 = 10L, Grule = as.integer(bt_rule_G(cell$n)))
  ex <- setdiff(ex, 10L)
  if (length(ex)) a <- c(a, setNames(ex, paste0("G", ex)))
  a
}

battery_names <- function(cell) {
  if (cell$type == "external") {
    ex <- unlist(lapply(bt_ext_G(cell), function(Gx) c(rbind(paste0(c(names(BT_EXT_GROUPED), "HL.ext"), ".G", Gx),
                                                            paste0("stat.", c(names(BT_EXT_GROUPED), "HL.ext"), ".G", Gx)))))
    return(c("n", "events", BT_EXT_TESTS, ex, bt_order_names(cell), "flag.rep_error"))
  }
  arms <- names(bt_arms(cell))
  ed <- character(0)
  for (a in arms) {
    for (b in BT_BASES) ed <- c(ed, paste0("EDGE.", b, ".u.", a), paste0("EDGE.", b, ".sc.", a))
    if (isTRUE(cell$ao)) ed <- c(ed, paste0("EDGE.ao.u.", a))
    if (cell$type == "edgeonly") ed <- c(ed, paste0("lam", 1:3, ".poly3.", a))
  }
  if (cell$type == "edgeonly")
    return(c("n", "events", "fit_ok", "glm_conv", ed, paste0("flag.evlt.", arms), "flag.degenerate", "flag.info_guard", "flag.rep_error"))
  c("n", "events", "fit_ok", "glm_conv", ed, paste0("HL.", arms), paste0("HLF.", arms),
    "Stk.joint", "Stk.LR", "Stk.sym1", "Stk.marg", "GiViTI", "GiViTI.t50", "Cubic.LR",
    "HL_w", "PH", "Tsiatis", "Xie", "PR", bt_order_names(cell),
    paste0("flag.evlt.", arms), "flag.stk_half0", "flag.stk_diag", "flag.rivals_run", "flag.degenerate", "flag.info_guard", "flag.rep_error")
}

## ---- one replicate ------------------------------------------------------------------------------------------------------------
bt_edge_arm <- function(fq, G, cell) {
  v <- c(); guard <- FALSE
  gs <- tryCatch(bt_groups(fq, G), error = function(e) NULL)
  ## EDGE is defined for a logit fit; for another link (real_beetle_cloglog) only the link-agnostic tests run, as in
  ## _realdata_direct.R
  bases <- if (identical(fq$fit$family$link, "logit")) c(BT_BASES, if (isTRUE(cell$ao)) "ao") else character(0)
  for (b in setdiff(c(BT_BASES, if (isTRUE(cell$ao)) "ao"), bases)) {
    v[paste0("EDGE.", b, ".u")] <- NA_real_
    if (b != "ao") v[paste0("EDGE.", b, ".sc")] <- NA_real_
  }
  for (b in bases) {
    Z <- if (is.null(gs)) NULL else tryCatch(bt_basis(gs$pbar, b), error = function(e) NULL)
    un <- if (is.null(Z)) NULL else tryCatch(bt_edge_unit(gs, fq$A, Z), error = function(e) NULL)
    v[paste0("EDGE.", b, ".u")] <- if (is.null(un)) NA_real_ else un$p
    if (b != "ao") {
      sc <- if (is.null(Z)) NULL else tryCatch(bt_edge_score(gs, fq$A, Z), error = function(e) NULL)
      v[paste0("EDGE.", b, ".sc")] <- if (is.null(sc)) NA_real_ else sc$p
      if (isTRUE(sc$guard)) guard <- TRUE
    }
    if (b == "poly3" && cell$type == "edgeonly") {
      lam <- if (is.null(un)) numeric(0) else un$lam
      v[paste0("lam", 1:3, ".poly3")] <- c(lam, NA, NA, NA)[1:3]
    }
  }
  if (cell$type == "full") {
    v["HL"]  <- if (is.null(gs)) NA_real_ else tryCatch(bt_hl_stat(fq$y, fq$ph, gs$grp), error = function(e) NA_real_)
    v["HLF"] <- tryCatch(bt_hlf(fq$y, fq$ph, G), error = function(e) NA_real_)
  }
  attr(v, "guard") <- guard
  v
}

battery_rep <- function(dat, cell) {
  nm <- battery_names(cell)
  out <- setNames(rep(NA_real_, length(nm)), nm)
  n <- nrow(dat$d); ev <- sum(dat$d$y)
  out["n"] <- n; out["events"] <- ev; out["fit_ok"] <- 0; out["flag.rep_error"] <- 0
  arms <- bt_arms(cell)
  for (a in names(arms)) out[paste0("flag.evlt.", a)] <- as.numeric(min(ev, n - ev) < arms[[a]])
  ## no event or no non-event: the working model has no maximum-likelihood fit, so no test has a p-value (E0.5
  ## clarification, review D1); the sample stays in every denominator and is counted in the declined column
  out["flag.degenerate"] <- as.numeric(min(ev, n - ev) == 0)
  if (out[["flag.degenerate"]] == 1) return(out)
  fq <- bt_fit(dat)
  if (is.null(fq)) return(out)
  out["fit_ok"] <- 1; out["glm_conv"] <- as.numeric(isTRUE(fq$fit$converged))

  cache <- list()
  for (a in names(arms)) {
    key <- as.character(arms[[a]])
    if (is.null(cache[[key]])) cache[[key]] <- bt_edge_arm(fq, arms[[a]], cell)
    v <- cache[[key]]
    out[paste0(names(v), ".", a)] <- v
  }
  ## flag.info_guard: a score statistic left out a column whose post-fit information is below 1e-10 of its unadjusted one
  eg <- as.numeric(any(vapply(cache, function(z) isTRUE(attr(z, "guard")), logical(1))))
  out["flag.info_guard"] <- eg
  if (identical(cell$generator, "real")) {                        # E8.5: ties in the fitted risks, in this row order
    out["info.distinct_risks"] <- length(unique(fq$ph))
    for (a in names(arms)) {
      G <- arms[[a]]
      grp <- pmin(ceiling(rank(fq$ph, ties.method = "first") / (fq$n / G)), G)     # def.gof's grouping, as bt_groups
      out[paste0("info.tied_boundaries.", a)] <- if (G > fq$n) NA_real_ else bt_tied_boundaries(fq$ph, grp)
      out[paste0("info.tied_boundaries_hlf.", a)] <- bt_tied_boundaries_hlf(fq$ph, G)
    }
  }
  if (cell$type == "edgeonly") return(out)

  if (identical(fq$fit$family$link, "logit")) {                   # Stukel and the cubic calibration LR need a logit fit
    st <- bt_stukel(fq)
    out[names(st)] <- st
    out["flag.info_guard"] <- max(eg, st[["flag.info_guard"]])
    out["Cubic.LR"] <- bt_cubic(fq$y, fq$eta)
  }
  out["GiViTI"]     <- bt_giviti(fq$y, fq$p_raw, 0.95)
  out["GiViTI.t50"] <- bt_giviti(fq$y, fq$p_raw, 0.50)
  if (n < BT_RIVAL_NMAX) {
    out["flag.rivals_run"] <- 1
    out["HL_w"]    <- tryCatch(bt_hlw(fq$y, fq$ph, 10), error = function(e) NA_real_)
    out["PH"]      <- tryCatch(bt_ph(fq$y, fq$ph, 10), error = function(e) NA_real_)
    out["Tsiatis"] <- tryCatch(suppressWarnings(bt_tsiatis(fq, 10)), error = function(e) NA_real_)   # k-means warnings only
    out["Xie"]     <- tryCatch(suppressWarnings(bt_xie(fq)), error = function(e) NA_real_)
    out["PR"]      <- tryCatch(bt_pr(fq), error = function(e) NA_real_)
  } else out["flag.rivals_run"] <- 0
  out
}

battery_one <- function(rep, cell) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  real <- identical(cell$generator, "real")
  if (real && is.na(cell$n)) cell$n <- nrow(bt_data(cell)$d)
  nm <- battery_names(cell)
  ## a real one-off fit (E8.5): rep 0 is the stored row order, rep k the order drawn after set.seed(BT_ORDER_SEED + k); the
  ## data are read first (the Diabetes-130 split sets its own seed when the file is read)
  seed <- if (real) { if (rep > 0) BT_ORDER_SEED + rep else NA_real_ } else if (is.na(cell$seed_base)) NA_real_ else cell$seed_base + rep
  if (!real && !is.na(seed)) set.seed(seed)
  v <- tryCatch({
    dat <- bt_data(cell)
    if (real && rep > 0) { set.seed(seed); dat <- bt_reorder(dat, sample.int(bt_nrow(dat))) }
    if (cell$type == "external") battery_rep_external(dat, cell) else battery_rep(dat, cell)
  }, error = function(e) { z <- setNames(rep(NA_real_, length(nm)), nm); z["flag.rep_error"] <- 1; z })
  c(rep = rep, seed = seed, v[nm])
}

## a chunk of replicates for one worker
bt_chunk <- function(reps, cell) do.call(rbind, lapply(reps, function(r) battery_one(r, cell)))
