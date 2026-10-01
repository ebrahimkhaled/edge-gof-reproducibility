## _block9b_rivals.R -- one replicate of block 9b: are the slow rivals robust to a corrupted record?
## Contract: paper_EDGE/theory/PREDECLARATION_block9b_slow_rivals_robustness.md (sha256 72d703a9...).
##
## Section 1, to the letter:
##   data   x ~ U(-3, 3), d ~ Bernoulli(0.5), eta = 0.6x + 0.5d, logistic truth, fitted as y ~ x + d, n = 500
##   cells  clean | x4 | x8 : in each replicate k = 3 records chosen at random have x multiplied by 4 or 8,
##          their outcome still drawn from the truth at the ORIGINAL x
##   tests  the projection test proj_pvalue(y, X, B = 250), X the model matrix; BAGofT at the PACKAGE
##          defaults (BAGofT 1.0.0: nsplits = 100, nsim = 100, ne = floor(5 sqrt(n)), parFun = parRF());
##          EDGE-poly3 and EDGE-sym unit at the rule G, Stukel's joint score, one-parameter score and LR
##          refit, GiViTI, the cubic calibration LR and HL at G = 10 -- all through the battery's own code
##   seeds  set.seed(400000000 + cell_id * 10000 + rep); the rivals' random numbers continue that stream
##
## The order inside a replicate is fixed and is part of the contract: data, then the battery's tests,
## then the projection test, then BAGofT. Seconds per data set are stored for the two slow rivals.
## BAGofT is the PACKAGE implementation, not the fast one: the author chose the declared form (2026-09-19).

B9B_SEED0 <- 400000000
B9B_N     <- 500L
B9B_K     <- 3L

b9b_cells <- function() {
  C <- data.frame(cell = c("clean", "x4", "x8"), mult = c(1, 4, 8), stringsAsFactors = FALSE)
  C$cell_id   <- seq_len(nrow(C))
  C$seed_base <- B9B_SEED0 + C$cell_id * 10000
  C$n <- B9B_N; C$k <- ifelse(C$mult == 1, 0L, B9B_K); C$B <- 100L
  C$type <- "full"; C$ao <- FALSE; C$G_extra <- ""; C$generator <- "b9b"
  C
}

b9b_gen <- function(mult, n = B9B_N, k = B9B_K) {
  x <- runif(n, -3, 3)
  d <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(0.6 * x + 0.5 * d))     # the outcome at the ORIGINAL x
  idx <- integer(0)
  if (mult != 1) { idx <- sample.int(n, k); x[idx] <- mult * x[idx] }
  list(d = data.frame(y = y, x = x, d = d), f = as.formula("y ~ x + d"), corrupt = idx)
}

b9b_names <- function(cell) c(battery_names(cell), "proj", "BAGofT", "BAGofT.p3",
                              "sec_proj", "sec_bagoft", "k_corrupt", "flag.b9b_error")

b9b_one <- function(rep, cell) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  nm   <- b9b_names(cell)
  seed <- cell$seed_base + rep
  set.seed(seed)
  v <- tryCatch({
    g   <- b9b_gen(cell$mult)
    dat <- list(d = g$d, f = g$f)
    out <- setNames(rep(NA_real_, length(nm)), nm)
    out[battery_names(cell)] <- battery_rep(dat, cell)[battery_names(cell)]
    out["k_corrupt"] <- length(g$corrupt)

    ## the projection test, on the model matrix of the fitted logistic model
    fit <- suppressWarnings(glm(g$f, data = g$d, family = binomial()))
    t0 <- proc.time()[["elapsed"]]
    out["proj"] <- tryCatch(proj_pvalue(g$d$y, stats::model.matrix(fit), B = 250)$p_value,
                            error = function(e) NA_real_)
    out["sec_proj"] <- proc.time()[["elapsed"]] - t0

    ## BAGofT at the package defaults; the model's own columns only
    t0 <- proc.time()[["elapsed"]]
    bg <- tryCatch(suppressWarnings(suppressMessages(BAGofT::BAGofT(
            testModel = BAGofT::testGlmBi(formula = y ~ x + d, link = "logit"),
            data = g$d, parFun = BAGofT::parRF()))), error = function(e) NULL)
    out["sec_bagoft"] <- proc.time()[["elapsed"]] - t0
    if (!is.null(bg)) { out["BAGofT"] <- bg$p.value; out["BAGofT.p3"] <- bg$p.value3 }

    out["flag.b9b_error"] <- 0
    out
  }, error = function(e) { z <- setNames(rep(NA_real_, length(nm)), nm); z["flag.b9b_error"] <- 1; z })
  c(rep = rep, seed = seed, v[nm])
}
