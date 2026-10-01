## map_pkg_proto.R -- shared by map_pkg_values.R and map_pkg_nullcheck.R (MAP_package.md, 2026-09-13).
## READ-ONLY with respect to the package: sources the ebrahim.gof dev tree into an environment `pk` and defines
## PROTOTYPES of the proposed changes -- (a) Stukel joint score u'I^-1 u with a one-df fallback, (b) EDGE basis
## "sym" and weights = "score". Nothing here edits, installs or rebuilds the package.
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
PKG <- Sys.getenv("EBRAHIM_GOF_SRC")
pk <- new.env()
for (f in list.files(file.path(PKG, "R"), pattern = "[.]R$", full.names = TRUE)) {
  r <- try(sys.source(f, envir = pk), silent = TRUE)
  if (inherits(r, "try-error")) cat("could not source", basename(f), "\n")
}
try(load(file.path(PKG, "R", "sysdata.rda"), envir = pk), silent = TRUE)
for (f in list.files(file.path(PKG, "data"), pattern = "[.]rda$", full.names = TRUE)) load(f, envir = pk)

## ---- (a) Stukel -------------------------------------------------------------------------------------------
## joint: u = Z'(y - p), I = Z'WZ - Z'WX (X'WX)^-1 X'WZ, T = u'I^-1 u on k df, k = number of non-zero half-columns.
## marginal: the CURRENT package statistic, called verbatim (sum of the two squared marginal z's on 2 df).
proto_stukel <- function(model, form = c("joint", "marginal")) {
  form <- match.arg(form)
  if (form == "marginal") return(c(pk$gof_stukel(pk$.gof_context(model)), list(fallback = NA, rho = NA_real_)))
  y   <- as.numeric(model$y)
  ph  <- pmin(pmax(as.numeric(stats::fitted(model)), 1e-6), 1 - 1e-6)
  X   <- stats::model.matrix(model)
  eta <- as.numeric(stats::predict(model, type = "link"))
  Z   <- cbind(upper = 0.5 * eta^2 * (ph >= 0.5), lower = -0.5 * eta^2 * (ph < 0.5))
  keep <- colSums(Z != 0) > 0
  if (!any(keep))
    return(list(Statistic = NA_real_, df = NA_real_, p_value = NA_real_,
                Note = "both Stukel columns are zero", fallback = NA, rho = NA_real_))
  Z <- Z[, keep, drop = FALSE]
  W <- ph * (1 - ph)
  u <- drop(crossprod(Z, y - ph))
  ZWX <- crossprod(Z, W * X)
  I <- crossprod(Z, W * Z) - ZWX %*% solve(crossprod(X, W * X), t(ZWX))
  Tst <- tryCatch(drop(crossprod(u, solve(I, u))), error = function(e) NA_real_)
  k <- ncol(Z)
  note <- if (k == 2) "" else if (keep[["upper"]])
    "one-df fallback: no fitted risk below 0.5" else "one-df fallback: no fitted risk at or above 0.5"
  list(Statistic = Tst, df = k,
       p_value = if (is.finite(Tst)) stats::pchisq(Tst, k, lower.tail = FALSE) else NA_real_,
       Note = note, fallback = (k == 1), rho = if (k == 2) stats::cov2cor(I)[1, 2] else NA_real_)
}

## ---- (b) EDGE ---------------------------------------------------------------------------------------------
proto_basis <- function(pbar, basis) {
  e <- stats::qlogis(pbar)
  if (basis == "sym") return(cbind(sym = e * abs(e)))                        # Stukel's symmetric carrier
  if (basis == "stukel2") return(cbind(e^2 * (e >= 0), -e^2 * (e < 0)))      # harness EDGE.stk.sc span (no linear column)
  pk$.def_basis(pbar, basis)
}

proto_def <- function(object, predicted_probs = NULL, X = NULL, G = 10,
                      basis = c("poly3", "poly2", "stukel", "sym", "stukel2"),
                      method = c("satterthwaite", "imhof"), weights = c("unit", "score"),
                      diag = FALSE) {
  basis <- match.arg(basis); method <- match.arg(method); weights <- match.arg(weights)
  if (inherits(object, "glm")) {
    y <- as.numeric(object$y); ph <- pmin(pmax(as.numeric(stats::fitted(object)), 1e-6), 1 - 1e-6)
    eta <- as.numeric(stats::predict(object, type = "link")); dmu <- object$family$mu.eta(eta)
    X <- stats::model.matrix(object); naive <- FALSE
  } else {
    y <- as.numeric(object); ph <- pmin(pmax(as.numeric(predicted_probs), 1e-6), 1 - 1e-6)
    dmu <- ph * (1 - ph); naive <- is.null(X)
  }
  n <- length(y); V <- ph * (1 - ph); w <- dmu^2 / V
  grp <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  idx <- split(seq_len(n), grp); Gn <- length(idx)
  og <- vapply(idx, function(I) sum(y[I]), numeric(1))
  eg <- vapply(idx, function(I) sum(ph[I]), numeric(1))
  Vg <- vapply(idx, function(I) sum(V[I]), numeric(1))
  pbar <- vapply(idx, function(I) mean(ph[I]), numeric(1))
  r <- (og - eg) / sqrt(Vg)
  if (naive) Omega <- diag(Gn) else {
    U <- t(vapply(idx, function(I) colSums(dmu[I] * X[I, , drop = FALSE]), numeric(ncol(X)))) / sqrt(Vg)
    Omega <- diag(Gn) - U %*% solve(crossprod(X, w * X)) %*% t(U)
  }
  Z <- proto_basis(pbar, basis)
  Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]
  if (weights == "unit") {
    ZtZ <- crossprod(Z); Zr <- crossprod(Z, r)
    S <- as.numeric(t(Zr) %*% solve(ZtZ) %*% Zr)
    lam <- Re(eigen(solve(ZtZ) %*% (t(Z) %*% Omega %*% Z), only.values = TRUE)$values)
    lam <- lam[lam > 1e-9]
    out <- data.frame(Test = "Directed Ebrahim-Farrington", Basis = basis, Test_Statistic = S,
                      df = sum(lam)^2 / sum(lam^2), Method = method,
                      p_value = pk$.def_pvalue(S, lam, method), stringsAsFactors = FALSE)
  } else {
    Zs <- Z * sqrt(Vg)                               # score form: column times sqrt(V_g)
    u  <- drop(crossprod(Zs, r))                     # = sum_g z_g (O_g - E_g)
    Iz <- crossprod(Zs, Omega %*% Zs)                # post-fit information of the weighted columns
    ev <- eigen((Iz + t(Iz)) / 2, symmetric = TRUE)
    ok <- ev$values > 1e-9 * max(ev$values)
    Tst <- sum(drop(crossprod(ev$vectors[, ok, drop = FALSE], u))^2 / ev$values[ok])
    k <- sum(ok)
    out <- data.frame(Test = "Directed Ebrahim-Farrington", Basis = basis, Test_Statistic = Tst,
                      df = k, Method = "score", p_value = stats::pchisq(Tst, k, lower.tail = FALSE),
                      stringsAsFactors = FALSE)
    if (diag) attr(out, "eig_rel") <- ev$values / max(ev$values)
  }
  out
}
