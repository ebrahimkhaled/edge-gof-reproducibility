## _blockC1b_contam.R -- the sign-error corruption of block Z Part A (declaration sha256 f5e1a7f1...).
##
## c1b_gen() is b9_gen() with one change: the corrupted records get -4x instead of 4x, so a high-risk
## patient is recorded at the bottom of the risk scale with their outcome intact -- a bad leverage point.
## Everything else is block 9's: the draw order (x, then d, then y at the ORIGINAL covariate, then the
## corruption), the truths, the clamps.
##
## The copy is checked rather than trusted: c1b_identity_check() runs both generators from the same seed
## with corruption "C1" and requires them to agree exactly. The runner calls it before any cell.
source(file.path(SIMDIR, "_block9_contam.R"))

C1B_SEED0 <- 700000000
C1B_RATES <- c(0.001, 0.002, 0.005, 0.01)
C1B_N     <- c(1000L, 5000L)

c1b_gen <- function(truth, n, corruption, rate) {
  x <- runif(n, -3, 3)
  d <- rbinom(n, 1, 0.5)
  eta <- 0.6 * x + 0.5 * d
  p   <- b9_inv(truth)(eta)                       # the true risk, at the ORIGINAL x
  y   <- rbinom(n, 1, p)

  k   <- if (identical(corruption, "clean")) 0L else b9_k(n, rate)
  idx <- integer(0)
  if (k > 0L) {
    if (corruption == "C1") {                     # block 9's own mechanism, kept for the identity gate
      idx <- sample.int(n, k)
      x[idx] <- 4 * x[idx]
    } else if (corruption == "C1b") {             # the sign error: the recorded position is reversed
      idx <- sample.int(n, k)
      x[idx] <- -4 * x[idx]
    } else stop("block C1b: unknown corruption '", corruption, "'")
  }
  list(d = data.frame(y = y, x = x, d = d), f = as.formula("y ~ x + d"), corrupt = idx, p_true = p)
}

## the identity gate: with corruption "C1" the copy must reproduce block 9's generator exactly
c1b_identity_check <- function(seeds = 1:5) {
  for (s in seeds) {
    set.seed(s); a <- b9_gen("logit", 1000L, "C1", 0.01)
    set.seed(s); b <- c1b_gen("logit", 1000L, "C1", 0.01)
    if (!isTRUE(all.equal(a$d, b$d)) || !isTRUE(all.equal(a$corrupt, b$corrupt)) ||
        !isTRUE(all.equal(a$p_true, b$p_true)))
      stop("block C1b: the generator copy differs from b9_gen() at seed ", s)
  }
  TRUE
}

c1b_cells <- function() {
  rows <- list()
  for (n in C1B_N) {
    for (rate in C1B_RATES)
      rows[[length(rows) + 1]] <- data.frame(truth = "logit", n = n, corruption = "C1b", rate = rate,
                                             stringsAsFactors = FALSE)
    rows[[length(rows) + 1]] <- data.frame(truth = "logit", n = n, corruption = "clean", rate = 0,
                                           stringsAsFactors = FALSE)
  }
  C <- do.call(rbind, rows)
  C$cell_id   <- seq_len(nrow(C))
  C$seed_base <- C1B_SEED0 + C$cell_id * 10000
  C$cell <- ifelse(C$corruption == "clean", sprintf("c1b_clean_n%d", C$n),
                   sprintf("c1b_r%s_n%d", sub("^0\\.", "", format(C$rate, trailing = FALSE,
                                                                  scientific = FALSE)), C$n))
  C$k <- ifelse(C$corruption == "clean", 0L, as.integer(round(C$rate * C$n)))
  ## the fields battery_rep() reads from a cell
  C$type <- "full"; C$ao <- FALSE; C$G_extra <- ""; C$generator <- "c1b"; C$block <- "C1b"
  C$role <- ifelse(C$corruption == "clean", "null", "alternative")
  C$null_cell <- sprintf("c1b_clean_n%d", C$n); C$null_block <- "C1b"
  if (nrow(C) != 10L) stop("block C1b: expected 10 cells, found ", nrow(C))
  C
}

## one replicate: the battery's own tests on a sign-corrupted data set
c1b_one <- function(rep, cell) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") RNGkind("L'Ecuyer-CMRG")
  nm <- battery_names(cell)
  set.seed(cell$seed_base + rep)
  v <- tryCatch({
    g <- c1b_gen(cell$truth, cell$n, cell$corruption, cell$rate)
    out <- setNames(rep(NA_real_, length(nm)), nm)
    out[nm] <- battery_rep(list(d = g$d, f = g$f), cell)[nm]
    c(rep = rep, k_corrupt = length(g$corrupt), out)
  }, error = function(e) c(rep = rep, k_corrupt = NA_real_, setNames(rep(NA_real_, length(nm)), nm)))
  v
}
