## _block9R2_hybrid.R -- the hybrid partition of block 9R2 (declaration sha256 9943e551...), all EDGE-poly3 unit form.
## Coarse tails, fine middle: the outer f of the risk scale at each end forms one group, the rest keeps the rule of about
## twenty-five records a group. Everything else -- the estimation adjustment, the basis, the null -- is block 9R's.
##   V0                     equal frequency at the rule G                 (the published test, and the identity gate)
##   HYB05, HYB10, HYB20    one group of floor(f*n) at each end, the middle in groups of about 25
##   G10                    equal frequency at G = 10                     (the anchor; must equal block 9R's G10)
source(file.path(SIMDIR, "_block9R_remedies.R"))          # fr_groups(), fr_equal(), fr_rule_G()

R9R2_TAIL <- c(HYB05 = 0.05, HYB10 = 0.10, HYB20 = 0.20)
R9R2_VARIANTS <- c("V0", "HYB05", "HYB10", "HYB20", "G10")

## the partition of one variant: a logical keep vector and the group labels of the kept records
fr2_partition <- function(fq, variant) {
  n <- fq$n; rk <- rank(fq$ph, ties.method = "first")
  if (variant == "V0")       { keep <- rep(TRUE, n); grp <- fr_equal(rk, fr_rule_G(n)) }
  else if (variant == "G10") { keep <- rep(TRUE, n); grp <- fr_equal(rk, 10L) }
  else if (variant %in% names(R9R2_TAIL)) {
    t <- as.integer(floor(R9R2_TAIL[[variant]] * n))
    if (t < 1L || 2L * t >= n - 50L) return(NULL)
    mid <- rk > t & rk <= n - t
    Gm <- max(4L, as.integer(round(sum(mid) / 25)))
    grp <- integer(n)
    grp[rk <= t]     <- 1L
    grp[mid]         <- 1L + fr_equal(rank(fq$ph[mid], ties.method = "first"), Gm)
    grp[rk > n - t]  <- Gm + 2L
    keep <- rep(TRUE, n)
  } else stop("unknown variant ", variant)
  list(keep = keep, grp = grp[keep])
}

## one replicate: the p-value of every variant on the same fitted model
fr2_pvalues <- function(fq) {
  out <- setNames(rep(NA_real_, length(R9R2_VARIANTS)), R9R2_VARIANTS)
  if (is.null(fq)) return(out)
  for (v in R9R2_VARIANTS) {
    p <- tryCatch({
      part <- fr2_partition(fq, v)
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

## the group sizes of one partition, for the self-test and the record of what was actually run
fr2_sizes <- function(fq, variant) {
  part <- fr2_partition(fq, variant)
  if (is.null(part)) return(NULL)
  as.integer(table(part$grp))
}
