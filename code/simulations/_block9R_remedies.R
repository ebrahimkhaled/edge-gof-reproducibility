## _block9R_remedies.R -- the six EDGE variants of block 9R (declaration sha256 7d79b85a...), all EDGE-poly3 unit form.
## Each takes the fitted quantities of bt_fit() and differs only in which records enter the statistic and how they are
## grouped; the estimation adjustment always keeps the full-sample information X'WX of the fit that made the predictions.
##   V0            equal frequency at the rule G                                    (the published test)
##   FR1, FR2      EDGE-FR: drop the floor(alpha*n) lowest and highest fitted risks, then the rule G on what is left
##   DROP1         group at the rule G, then drop the ceil(alpha*G) outer groups at each end
##   POOL1         the outer floor(alpha*n) at each end form one group each; the middle in groups of about 25
##   G10           equal frequency at G = 10

R9R_ALPHA <- c(FR1 = 0.01, FR2 = 0.02, DROP1 = 0.01, POOL1 = 0.01)
R9R_VARIANTS <- c("V0", "FR1", "FR2", "DROP1", "POOL1", "G10")

## the grouped quantities of bt_groups(), for an arbitrary partition of an arbitrary subset of the records
fr_groups <- function(fq, keep, grp) {
  V <- fq$ph * (1 - fq$ph)
  idx <- split(which(keep), grp)
  idx <- idx[vapply(idx, length, 1L) > 0L]
  if (length(idx) < 5L) return(NULL)
  og   <- vapply(idx, function(I) sum(fq$y[I]),  numeric(1))
  eg   <- vapply(idx, function(I) sum(fq$ph[I]), numeric(1))
  Vg   <- vapply(idx, function(I) sum(V[I]),     numeric(1))
  if (any(!is.finite(Vg)) || any(Vg <= 0)) return(NULL)
  pbar <- vapply(idx, function(I) mean(fq$ph[I]), numeric(1))
  U <- t(vapply(idx, function(I) colSums(fq$dmu[I] * fq$X[I, , drop = FALSE]), numeric(ncol(fq$X)))) / sqrt(Vg)
  list(G = length(idx), ng = lengths(idx), og = og, eg = eg, Vg = Vg, pbar = pbar,
       r = (og - eg) / sqrt(Vg), U = U, Omega = diag(length(idx)) - U %*% solve(fq$A) %*% t(U))
}
fr_rule_G <- function(n) max(10L, as.integer(ceiling(n / 25)))
fr_equal <- function(rank_in, G) pmin(ceiling(rank_in / (length(rank_in) / G)), G)   # bt_groups' equal-frequency rule

## the partition of one variant: a logical keep vector and the group labels of the kept records
fr_partition <- function(fq, variant) {
  n <- fq$n; rk <- rank(fq$ph, ties.method = "first")
  if (variant == "V0")  { keep <- rep(TRUE, n); grp <- fr_equal(rk, fr_rule_G(n)) }
  else if (variant == "G10") { keep <- rep(TRUE, n); grp <- fr_equal(rk, 10L) }
  else if (variant %in% c("FR1", "FR2")) {
    a <- R9R_ALPHA[[variant]]; t <- floor(a * n)
    keep <- rk > t & rk <= n - t
    if (sum(keep) < 50L) return(NULL)
    grp <- fr_equal(rank(fq$ph[keep], ties.method = "first"), fr_rule_G(sum(keep)))
  } else if (variant == "DROP1") {
    G <- fr_rule_G(n); j <- ceiling(R9R_ALPHA[["DROP1"]] * G)
    g <- fr_equal(rk, G)
    keep <- g > j & g <= G - j
    if (sum(keep) < 50L) return(NULL)
    grp <- g[keep]
  } else if (variant == "POOL1") {
    t <- floor(R9R_ALPHA[["POOL1"]] * n)
    if (t < 1L) { keep <- rep(TRUE, n); grp <- fr_equal(rk, fr_rule_G(n)) }
    else {
      mid <- rk > t & rk <= n - t
      Gm <- max(8L, as.integer(round(sum(mid) / 25)))
      grp <- integer(n)
      grp[rk <= t] <- 1L
      grp[mid] <- 1L + fr_equal(rank(fq$ph[mid], ties.method = "first"), Gm)
      grp[rk > n - t] <- Gm + 2L
      keep <- rep(TRUE, n)
    }
  } else stop("unknown variant ", variant)
  list(keep = keep, grp = grp[keep])
}

## one replicate: the p-value of every variant on the same fitted model
fr_pvalues <- function(fq) {
  out <- setNames(rep(NA_real_, length(R9R_VARIANTS)), R9R_VARIANTS)
  if (is.null(fq)) return(out)
  for (v in R9R_VARIANTS) {
    p <- tryCatch({
      part <- fr_partition(fq, v)
      if (is.null(part)) NA_real_ else {
        gs <- fr_groups(fq, part$keep, part$grp)
        Z <- if (is.null(gs)) NULL else bt_basis(gs$pbar, "poly3")
        if (is.null(Z)) NA_real_ else bt_edge_unit(gs, fq$A, Z)$p
      }
    }, error = function(e) NA_real_)
    out[v] <- if (is.finite(p) && p >= 0 && p <= 1) p else NA_real_
  }
  out
}
