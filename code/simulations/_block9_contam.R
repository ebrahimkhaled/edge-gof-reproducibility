## _block9_contam.R -- one replicate of block 9 (paper_EDGE/theory/PREDECLARATION_block9_contamination.md,
## sha256 99eb4c27...). b9_one(rep, cell) seeds set.seed(cell$seed_base + rep), draws the base design, corrupts a
## handful of records as the declaration says, fits y ~ x + d once and returns the battery's own test vector plus
## the two extra quantities section 1 asks for: the fitted coefficients, and how many corrupted records land in
## the highest-risk group.
##
## Nothing here changes _battery_tests.R. The tests themselves are battery_rep(), used read-only, so block 9 and
## blocks 0-8 compute the same statistics by the same code. battery_rep() also returns HL_w, PH, Tsiatis, Xie and
## PR, which block 9 does not claim anything about; they are stored because they cost one pass over data that is
## already in memory, and the analysis uses only the declared list.
##
## Draw order is fixed and must not be changed: x, then d, then y from the truth at the ORIGINAL x, then the
## corruption. Corrupting after y is drawn is what makes C1 and C3 "the prediction is wrong, the outcome is not".

B9_TRUTHS <- c("logit", "probit", "cloglog")          # (a), (b), (c) of section 1
B9_RATES  <- c(0.001, 0.002, 0.005, 0.01)
B9_N      <- c(1000L, 5000L)
B9_SEED0  <- 300000000                                 # section 1: seed_base = 300000000 + cell_id * 10000

## The inverse links of the three truths, taken from the battery's own table so that "probit" here and "probit"
## in block 1 mean the same function.
b9_inv <- function(truth) {
  f <- BT_LINKINV[[truth]]
  if (is.null(f)) stop("block 9: unknown truth '", truth, "'")
  f
}

b9_k <- function(n, rate) as.integer(round(rate * n))

## One contaminated data set. Returns the battery's dat list plus the rows that were corrupted.
b9_gen <- function(truth, n, corruption, rate) {
  x <- runif(n, -3, 3)
  d <- rbinom(n, 1, 0.5)
  eta <- 0.6 * x + 0.5 * d
  p   <- b9_inv(truth)(eta)                       # the true risk, at the ORIGINAL x
  y   <- rbinom(n, 1, p)

  k   <- if (identical(corruption, "clean")) 0L else b9_k(n, rate)
  idx <- integer(0)
  if (k > 0L) {
    if (corruption == "C1") {                     # corrupted covariate: prediction wrong, outcome right
      idx <- sample.int(n, k)
      x[idx] <- 4 * x[idx]
    } else if (corruption == "C3") {              # unit error: the record's prediction becomes milder
      idx <- sample.int(n, k)
      x[idx] <- x[idx] / 10
    } else if (corruption == "C2") {              # missed events: prediction right, recorded outcome wrong
      idx <- order(p, decreasing = TRUE)[seq_len(k)]
      y[idx] <- 0L                                # 0L, not 0: a double would retype the whole outcome vector
    } else if (corruption == "C4") {              # duplicated high-risk records, over randomly chosen rows
      src <- order(p, decreasing = TRUE)[seq_len(k)]
      idx <- sample.int(n, k)
      x[idx] <- x[src]; d[idx] <- d[src]; y[idx] <- y[src]
    } else stop("block 9: unknown corruption '", corruption, "'")
  }
  list(d = data.frame(y = y, x = x, d = d), f = as.formula("y ~ x + d"), corrupt = idx, p_true = p)
}

## The extra columns of section 1, beyond what battery_names() covers.
b9_extra_names <- function(cell) c("k_corrupt", "b.intercept", "b.x", "b.d",
                                   paste0("corrupt_in_top.", names(bt_arms(cell))),
                                   "flag.b9_error")

b9_names <- function(cell) c(battery_names(cell), b9_extra_names(cell))

## How many of the corrupted records fall in the highest-risk group, under each G arm. The grouping is the one
## EDGE uses: equal frequency by fitted risk, ties split by row order.
b9_top_group <- function(ph, idx, G) {
  n <- length(ph)
  if (!length(idx) || G < 1) return(0)
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  sum(grp[idx] == G)
}

b9_one <- function(rep, cell) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  nm   <- b9_names(cell)
  seed <- cell$seed_base + rep
  set.seed(seed)
  v <- tryCatch({
    g   <- b9_gen(cell$truth, cell$n, cell$corruption, cell$rate)
    out <- setNames(rep(NA_real_, length(nm)), nm)
    out[battery_names(cell)] <- battery_rep(list(d = g$d, f = g$f), cell)[battery_names(cell)]
    out["k_corrupt"]   <- length(g$corrupt)
    out["flag.b9_error"] <- 0
    ## the fit again, for the coefficients and the group membership; one glm beside the whole test battery
    fit <- tryCatch(suppressWarnings(glm(g$f, data = g$d, family = binomial())), error = function(e) NULL)
    if (!is.null(fit)) {
      cf <- coef(fit)
      out["b.intercept"] <- unname(cf[1]); out["b.x"] <- unname(cf[2]); out["b.d"] <- unname(cf[3])
      ph <- pmin(pmax(as.numeric(fitted(fit)), 1e-6), 1 - 1e-6)
      for (a in names(bt_arms(cell)))
        out[paste0("corrupt_in_top.", a)] <- b9_top_group(ph, g$corrupt, bt_arms(cell)[[a]])
    }
    out
  }, error = function(e) {
    z <- setNames(rep(NA_real_, length(nm)), nm); z["flag.b9_error"] <- 1; z
  })
  c(rep = rep, seed = seed, v[nm])
}

## ---- the 84 cells, in a fixed order that fixes cell_id and therefore every seed ---------------------------
## 3 truths x 3 corruptions (C1, C2, C3) x 4 rates x 2 n = 72, plus C4 at 0.005 only (3 x 2 = 6), plus clean
## (3 x 2 = 6). The order below is the order of the declaration's own list and must not be reshuffled: cell_id
## is its position, and seed_base = 300000000 + cell_id * 10000.
b9_cells <- function() {
  rows <- list()
  for (truth in B9_TRUTHS) for (n in B9_N) {
    for (corr in c("C1", "C2", "C3")) for (rate in B9_RATES)
      rows[[length(rows) + 1]] <- data.frame(truth = truth, n = n, corruption = corr, rate = rate,
                                             stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(truth = truth, n = n, corruption = "C4", rate = 0.005,
                                           stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(truth = truth, n = n, corruption = "clean", rate = 0,
                                           stringsAsFactors = FALSE)
  }
  C <- do.call(rbind, rows)
  C$cell_id   <- seq_len(nrow(C))
  C$seed_base <- B9_SEED0 + C$cell_id * 10000
  C$cell <- ifelse(C$corruption == "clean",
                   sprintf("%s_clean_n%d", C$truth, C$n),
                   sprintf("%s_%s_r%s_n%d", C$truth, C$corruption,
                           sub("^0\\.", "", format(C$rate, trailing = FALSE, scientific = FALSE)), C$n))
  C$k       <- ifelse(C$corruption == "clean", 0L, as.integer(round(C$rate * C$n)))
  C$block   <- "9"
  C$B       <- 1000L
  C$type    <- "full"
  C$ao      <- FALSE
  C$G_extra <- ""
  C$generator <- "b9"
  C$formula <- "y ~ x + d"
  ## the clean cell each corrupted cell is compared with (rule 2 of section 2)
  C$clean_cell <- sprintf("%s_clean_n%d", C$truth, C$n)
  C[, c("block", "cell", "cell_id", "truth", "n", "corruption", "rate", "k", "B", "seed_base",
        "type", "ao", "G_extra", "generator", "formula", "clean_cell")]
}
