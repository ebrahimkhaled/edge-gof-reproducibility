## pkg280_verify_fix_engine.R -- one copy of ebrahim.gof, loaded with pkgload::load_all(), run on a fixed set of fits:
##   def.gof unit form (poly2 / poly3 / stukel / sym x Satterthwaite / Imhof x G = 10 and the "auto" number),
##   def.gof score form (same bases and G), the (y, ph, X) path, def.ensemble.gof, and the Stukel row of
##   run.all.gof in every form the copy has.
## Fits: the 9 fixtures of pkg280_reference.R, 7 fixtures of pkg280_verify_num_exact.R, 200 random fits each of the
## base (n 1000), c0 = -2 (n 500) and s = 2 (n 1000) designs, and 2000 more c0 = -2 fits for the stukel basis only.
## Run once per copy, then pkg280_verify_fix_compare.R. Verification of the package, not a paper result.
## Run: Rscript pkg280_verify_fix_engine.R <package dir> <label> <out dir>

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
PKG <- args[1]; LAB <- args[2]; ODIR <- args[3]
options(warn = 1)
suppressPackageStartupMessages(pkgload::load_all(PKG, export_all = FALSE, quiet = TRUE))
sha   <- system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE)
ver   <- as.character(utils::packageVersion("ebrahim.gof"))
has_w <- "weights" %in% names(formals(def.gof))
cat(LAB, "| ebrahim.gof", ver, "| HEAD", sha, "| weights argument:", has_w, "|", R.version.string, "|",
    format(Sys.time()), "\n")
t0 <- Sys.time()
lap <- function(label) cat(sprintf("  %-28s done at %.1f min\n", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))

## ---- fixtures of pkg280_reference.R (same code, same seeds) ----
make_fit <- function(seed = 1, n = 500, link = "logit") {
  set.seed(seed)
  x <- runif(n, -3, 3)
  eta <- 0.6 * x
  p <- if (link == "cloglog") 1 - exp(-exp(eta)) else 1 / (1 + exp(-eta))
  glm(rbinom(n, 1, p) ~ x, family = binomial())
}
gof_demo <- ebrahim.gof::gof_demo
wrong <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial())
right <- glm(outcome ~ poly(age, 2) + bmi + sex + treatment, data = gof_demo, family = binomial())
set.seed(42)
x1 <- rnorm(400); x2 <- runif(400, -2, 2)
fit42 <- glm(rbinom(400, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2)) ~ x1 + x2, family = binomial())
set.seed(5)
x1 <- rnorm(400); d <- factor(sample(c("A", "B"), 400, replace = TRUE))
fitcat <- glm(rbinom(400, 1, plogis(0.3 + 0.6 * x1 + ifelse(d == "B", 0.5, 0))) ~ x1 + d, family = binomial())
REF <- list(mf1_500 = make_fit(1, 500), mf1_600 = make_fit(1, 600), mf2_600 = make_fit(2, 600),
            mf3_600 = make_fit(3, 600), cloglog7 = make_fit(7, 1500, "cloglog"),
            seed42 = fit42, cat5 = fitcat, demo_wrong = wrong, demo_right = right)

## ---- fixtures of pkg280_verify_num_exact.R (same code, same seeds) ----
mk <- function(seed, n, c0 = 0, s = 1, truth = "logit", eta_fun = NULL) {
  set.seed(seed)
  dat <- data.frame(xa = runif(n, -3, 3), db = rbinom(n, 1, 0.5))
  eta <- if (is.null(eta_fun)) c0 + s * (0.6 * dat$xa + 0.5 * dat$db) else eta_fun(dat)
  dat$out <- rbinom(n, 1, if (truth == "logit") plogis(eta) else 1 - exp(-exp(eta)))
  suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat))
}
find_fit <- function(seeds, make, cond) {
  for (sd in seeds) { f <- make(sd); if (cond(f)) return(f) }
  stop("find_fit: no seed satisfied the condition")
}
EXA <- list(
  base      = mk(101, 1000),
  c0_m2     = mk(102, 500, c0 = -2),
  s_2       = mk(103, 1000, s = 2),
  cloglog   = mk(104, 1500, truth = "cloglog"),
  all_below = find_fit(5:500, function(sd) mk(sd, 400, eta_fun = function(d) -4 + 0.5 * d$xa),
                       function(f) max(fitted(f)) < 0.5),
  all_above = find_fit(6:500, function(sd) mk(sd, 400, eta_fun = function(d) 4 + 0.5 * d$xa),
                       function(f) min(fitted(f)) >= 0.5),
  offset    = local({ set.seed(108); dat <- data.frame(xa = runif(800, -3, 3), db = rbinom(800, 1, 0.5))
                      dat$out <- rbinom(800, 1, plogis(-0.5 + 0.7 * dat$xa + 0.4 * dat$db))
                      glm(out ~ xa + offset(0.4 * db), family = binomial(), data = dat) })
)

## ---- random data sets, all drawn before any package function runs ----
gen_data <- function(n, c0, s) {
  dat <- data.frame(xa = runif(n, -3, 3), db = rbinom(n, 1, 0.5))
  dat$out <- rbinom(n, 1, plogis(c0 + s * (0.6 * dat$xa + 0.5 * dat$db)))
  dat
}
set.seed(880001); D_base  <- replicate(200,  gen_data(1000, 0, 1), simplify = FALSE)
set.seed(880002); D_c0m2  <- replicate(200,  gen_data(500, -2, 1), simplify = FALSE)
set.seed(880003); D_s2    <- replicate(200,  gen_data(1000, 0, 2), simplify = FALSE)
set.seed(880004); D_extra <- replicate(2000, gen_data(500, -2, 1), simplify = FALSE)
fit_of <- function(dat) suppressWarnings(glm(out ~ xa + db, family = binomial(), data = dat))
lap("fixtures and data")

## ---- one fit ----
cap <- function(expr) tryCatch(as.list(suppressWarnings(expr)), error = function(e) paste("ERROR:", conditionMessage(e)))
near_half <- function(fit, G) {                      # the group mean risks closest to one half
  ph  <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6); n <- length(ph)
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  pbar <- as.numeric(tapply(ph, grp, mean))
  up <- pbar[pbar >= 0.5]; lo <- pbar[pbar < 0.5]
  c(n_up = length(up), min_up = if (length(up)) min(up) else NA_real_,
    max_up = if (length(up)) max(up) else NA_real_, max_lo = if (length(lo)) max(lo) else NA_real_)
}
one <- function(fit, full = TRUE) {
  res <- list()
  n <- length(fit$y); Ga <- max(10, round(n / 25))
  bases <- if (full) c("poly2", "poly3", "stukel", "sym") else "stukel"
  for (G in unique(c(10, Ga))) {
    res[[sprintf("pbar|G%d", G)]] <- near_half(fit, G)
    for (b in bases) {
      for (m in if (full) c("satterthwaite", "imhof") else "satterthwaite")
        res[[sprintf("unit|%s|G%d|%s", b, G, m)]] <- cap(def.gof(fit, G = G, basis = b, method = m))
      if (has_w) res[[sprintf("score|%s|G%d", b, G)]] <- cap(def.gof(fit, G = G, basis = b, weights = "score"))
    }
  }
  if (!full) return(res)
  y <- as.numeric(fit$y); ph <- as.numeric(fitted(fit)); X <- model.matrix(fit)
  for (b in c("poly3", "stukel")) {
    res[[sprintf("yphX unit|%s|G10", b)]] <- cap(def.gof(y, ph, X = X, G = 10, basis = b))
    if (has_w) res[[sprintf("yphX score|%s|G10", b)]] <- cap(def.gof(y, ph, X = X, G = 10, basis = b, weights = "score"))
  }
  res[["ensemble|cct|G10"]]    <- cap(def.ensemble.gof(fit))
  res[["ensemble|cct+EF|G10"]] <- cap(def.ensemble.gof(fit, add_ef = TRUE))
  stk <- function(ctl) cap(as.data.frame(run.all.gof(fit, tests = "Stukel", install = "no", control = ctl)))
  res[["stukel|default"]] <- stk(list())
  if (has_w) for (f in c("marginal", "joint", "lr")) res[[paste0("stukel|", f)]] <- stk(list(Stukel = list(form = f)))
  res
}

R <- list(meta = list(label = LAB, version = ver, head = sha, has_w = has_w, time = Sys.time(), R = R.version.string))
R$ref9  <- lapply(REF, one);                                lap("ref9")
R$exa7  <- lapply(EXA, one);                                lap("exa7")
R$base  <- lapply(D_base, function(d) one(fit_of(d)));      lap("base 200")
R$c0m2  <- lapply(D_c0m2, function(d) one(fit_of(d)));      lap("c0 = -2, 200")
R$s2    <- lapply(D_s2,   function(d) one(fit_of(d)));      lap("s = 2, 200")
R$extra <- lapply(D_extra, function(d) one(fit_of(d), full = FALSE)); lap("c0 = -2 stukel only, 2000")
file <- file.path(ODIR, sprintf("engine_%s.rds", LAB))
saveRDS(R, file)
cat("saved", file, "\n")
