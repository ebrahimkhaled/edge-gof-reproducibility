## pkg280_verify_corr_guard.R -- verification of the no-fit check (2e0f6d1) and the information guard (7bc4769).
## Part A. Samples with no event or no non-event: def.gof (glm, (y, ph, X) and (y, ph) paths, every basis and form),
##         edge.gof, def.ensemble.gof, the Stukel row in every form and the fast run.all.gof battery. Controls with one
##         event and one non-event must not be called degenerate.
## Part B. How often the guard fires on 5000 null samples each of gen_sparse_link(n, -4.9) at n = 300 and 500 and the base
##         design (dgp_null("link", 1000)), plus n = 100 sparse as an extra. The guard is found two ways: by comparing the
##         installed build with the code one commit before the guard (2e0f6d1, sourced into its own environment, with the
##         ratio diag(I) / colSums(Zs^2) recorded), and by recomputing the Stukel ratio independently. For every firing
##         the p-value that the guard removed is kept, and the unguarded statistic is recomputed on a warm-restarted fit
##         (glm epsilon 1e-13) to see whether it was a number or rounding noise. The battery harness (bt_fit, bt_groups,
##         bt_basis, bt_edge_score, bt_stukel) is compared with the installed build on the same samples (E8.2).
## Seeds are new (74e6 ... 79e6 + r). Verification of the package, not a paper result.
## Run: Rscript pkg280_verify_corr_guard.R <clone dir> <workers> > ../paper_EDGE/theory/pkg280_verify_corr_guard.log 2>&1
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

args <- commandArgs(trailingOnly = TRUE)
PKG <- args[1]; NW <- as.integer(args[2])
SIM <- edge_path("code/simulations")
OUT <- edge_path("declarations")
SCR <- file.path(tempdir(), "vcorr")
options(width = 230, warn = 1)
library(parallel)
t0 <- Sys.time()
lap <- function(label) cat(sprintf("  [%s done at %.1f min]\n", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
git_show <- function(sha, file) system2("git", c("-C", shQuote(PKG), "show", paste0(sha, ":", file)), stdout = TRUE)
SRC <- list(def = git_show("2e0f6d1", "R/def_gof.R"), run = git_show("2e0f6d1", "R/run_all_gof.R"), SIM = SIM)

init <- function(SRC) {
  suppressPackageStartupMessages(library(ebrahim.gof))
  ns <- asNamespace("ebrahim.gof")
  txt <- SRC$def
  i <- grep("I  <- crossprod(Zs, Omega %*% Zs)", txt, fixed = TRUE)
  stopifnot(length(i) == 1, length(grep("ok <- d > 0$", txt)) == 1)
  txt <- append(txt, "    assign('rec', list(dI = diag(I), zz = colSums(Zs^2)), envir = V_REC)", after = i)
  PREG <- new.env(parent = ns)
  for (ex in parse(text = txt, keep.source = FALSE))
    if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]])) eval(ex, PREG)
  PRES <- new.env(parent = ns)
  for (ex in parse(text = SRC$run, keep.source = FALSE))
    if (is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]) && as.character(ex[[2]]) == "gof_stukel") eval(ex, PRES)
  assign("PREG", PREG, envir = globalenv()); assign("PRES", PRES, envir = globalenv())
  assign("V_REC", new.env(), envir = globalenv())
  source(file.path(SRC$SIM, "_battery_tests.R"))           # harness functions and _dgp_library.R, into the global env
  NULL
}

cap <- function(expr) {
  cls <- character(0); msg <- character(0)
  val <- tryCatch(withCallingHandlers(expr, warning = function(w) {
           cls <<- c(cls, class(w)[1]); msg <<- c(msg, conditionMessage(w)); invokeRestart("muffleWarning") }),
         error = function(e) structure(list(msg = conditionMessage(e)), class = "v_err"))
  list(val = val, cls = cls, msg = msg)
}
is_err <- function(x) inherits(x$val, "v_err")

## ======================================== Part A ========================================
partA <- function() {
  ns <- asNamespace("ebrahim.gof")
  res <- list()
  rec <- function(sample, what, ok, detail = "") res[[length(res) + 1]] <<- data.frame(sample = sample, check = what, ok = ok, detail = detail)
  mk_sample <- function(n, des, ny) {                 # ny = number of events
    if (des == "sparse") { x <- as.numeric(scale(rchisq(n, 4))); d <- data.frame(x = x); fo <- y ~ x }
    else { d <- data.frame(x = runif(n, -3, 3), dd = rbinom(n, 1, 0.5)); fo <- y ~ x + dd }
    y <- rep(0, n); if (ny > 0) y[order(d$x, decreasing = TRUE)[seq_len(ny)]] <- 1
    d$y <- y
    list(d = d, fit = suppressWarnings(glm(fo, family = binomial(), data = d)))
  }
  set.seed(74000000)
  for (n in c(30L, 200L, 1000L)) for (des in c("sparse", "base")) for (ny in c(0L, n)) {
    s <- mk_sample(n, des, ny); fit <- s$fit; y <- fit$y
    ph <- pmin(pmax(fitted(fit), 1e-6), 1 - 1e-6); X <- model.matrix(fit)
    lab <- sprintf("%s n=%d %s", des, n, if (ny == 0) "all 0" else "all 1")
    one_def <- function(what, expr, extra_ok = character(0)) {
      z <- cap(expr)
      dg <- sum(z$cls == "def_degenerate"); other <- setdiff(z$cls, c("def_degenerate", extra_ok))
      ok <- !is_err(z) && is.data.frame(z$val) && all(is.na(z$val$p_value)) && dg == 1 && length(other) == 0 &&
            (is.null(z$val$Test_Statistic) || all(is.na(z$val$Test_Statistic)))
      rec(lab, what, ok, if (ok) "" else paste(c(if (is_err(z)) z$val$msg, z$cls, if (!is_err(z)) format(z$val$p_value)), collapse = " | "))
    }
    for (b in c("poly2", "poly3", "stukel", "sym")) for (fm in c("sat", "imh", "sco")) {
      argsl <- list(G = 10, basis = b)
      if (fm == "imh") argsl$method <- "imhof"; if (fm == "sco") argsl$weights <- "score"
      one_def(paste("def.gof glm", b, fm), do.call(def.gof, c(list(fit), argsl)))
      one_def(paste("def.gof (y,ph,X)", b, fm), do.call(def.gof, c(list(y, ph, X = X), argsl)))
      one_def(paste("def.gof (y,ph) naive", b, fm), do.call(def.gof, c(list(y, ph), argsl)), extra_ok = "simpleWarning")
      one_def(paste("edge.gof glm", b, fm), do.call(edge.gof, c(list(fit), argsl)))
      one_def(paste("def.gof glm G auto", b, fm), do.call(def.gof, c(list(fit), modifyList(argsl, list(G = "auto")))))
    }
    one_def("def.gof glm ensemble unit", def.gof(fit, basis = "ensemble"))
    one_def("def.gof glm ensemble score", def.gof(fit, basis = "ensemble", weights = "score"))
    one_def("def.ensemble.gof glm", def.ensemble.gof(fit))
    one_def("def.ensemble.gof glm add_ef", def.ensemble.gof(fit, add_ef = TRUE))
    one_def("def.ensemble.gof (y,ph,X)", def.ensemble.gof(y, ph, X = X))
    for (fm in c("joint", "lr", "marginal")) {
      z <- cap(run.all.gof(fit, tests = "Stukel", control = list(Stukel = list(form = fm))))
      ok <- !is_err(z) && is.na(z$val$p_value[1]) && grepl("no maximum-likelihood fit", z$val$Note[1]) && length(z$cls) == 0
      rec(lab, paste("Stukel row", fm), ok, if (ok) z$val$Note[1] else paste(c(z$cls, z$val$Note[1], z$val$p_value[1]), collapse = " | "))
    }
    for (wt in c("unit", "score")) {
      ctl <- if (wt == "score") setNames(rep(list(list(weights = "score")), 4), c("DEF.poly2", "DEF.poly3", "DEF.stukel", "DEF.sym")) else list()
      z <- cap(run.all.gof(fit, include_slow = FALSE, control = ctl))
      if (is_err(z)) { rec(lab, paste("fast battery", wt), FALSE, z$val$msg); next }
      tb <- z$val
      dr <- tb[tb$Test %in% c("DEF.poly2", "DEF.poly3", "DEF.stukel", "DEF.sym"), ]
      er <- tb[grepl("^Ensemble", tb$Test), ]
      sr <- tb[tb$Test == "Stukel", ]
      rec(lab, paste("fast battery", wt, "DEF rows NA with no-fit Note"),
          nrow(dr) == 4 && all(is.na(dr$p_value)) && all(grepl("no maximum-likelihood fit", dr$Note)), paste(dr$Note, collapse = " / "))
      rec(lab, paste("fast battery", wt, "Stukel row NA with note"), nrow(sr) == 1 && is.na(sr$p_value) &&
          grepl("no maximum-likelihood fit", sr$Note), sr$Note)
      rec(lab, paste("fast battery", wt, "ensemble rows NA"), nrow(er) >= 1 && all(is.na(er$p_value)),
          paste(er$Test, er$p_value, collapse = " / "))
      esc <- unique(z$msg)
      rec(lab, paste("fast battery", wt, "warnings escaping"), !any(grepl("^def.gof", esc)),
          if (length(esc)) paste(substr(esc, 1, 90), collapse = " / ") else "none")
    }
  }
  ## controls: one event / one non-event is not degenerate
  for (n in c(30L, 200L)) for (des in c("sparse", "base")) for (ny in c(1L, n - 1L)) {
    s <- mk_sample(n, des, ny); fit <- s$fit
    lab <- sprintf("control %s n=%d events=%d", des, n, ny)
    z1 <- cap(def.gof(fit, basis = "stukel", weights = "score")); z2 <- cap(def.gof(fit, basis = "sym"))
    z3 <- cap(run.all.gof(fit, tests = "Stukel"))
    ok <- !any(c(z1$cls, z2$cls) == "def_degenerate") && !grepl("no maximum-likelihood", z3$val$Note[1])
    rec(lab, "not called degenerate", ok, paste(c(unique(c(z1$cls, z2$cls)), if (is_err(z1)) z1$val$msg, if (is_err(z2)) z2$val$msg,
                                                  z3$val$Note[1]), collapse = " | "))
  }
  do.call(rbind, res)
}

## ======================================== Part B ========================================
v_gen <- function(design, r) {
  switch(design,
    sparse100  = { set.seed(79000000 + r); gen_sparse_link(100, -4.9) },
    sparse300  = { set.seed(75000000 + r); gen_sparse_link(300, -4.9) },
    sparse500  = { set.seed(76000000 + r); gen_sparse_link(500, -4.9) },
    base1000   = { set.seed(77000000 + r); dgp_null("link", 1000) })
}
v_fit <- function(g, tight = FALSE) {
  if (!tight) return(suppressWarnings(glm(g$f, family = binomial(), data = g$d)))
  ctl <- glm.control(epsilon = 1e-13, maxit = 100)
  f0 <- suppressWarnings(glm(g$f, family = binomial(), data = g$d, control = ctl))
  suppressWarnings(glm(g$f, family = binomial(), data = g$d, control = ctl, start = unname(coef(f0))))
}
v_stk_ratio <- function(fit) {                           # independent: Z'WZ - Z'WX (X'WX)^-1 X'WZ on raw risks, over Z'WZ
  ph <- pmin(pmax(fitted(fit), 1e-6), 1 - 1e-6); pr <- as.numeric(fitted(fit)); W <- pr * (1 - pr)
  eta <- fit$linear.predictors
  Z <- cbind(0.5 * eta^2 * (ph >= 0.5), -0.5 * eta^2 * (ph < 0.5)); Z <- Z[, colSums(Z != 0) > 0, drop = FALSE]
  X <- model.matrix(fit)[, !is.na(coef(fit)), drop = FALSE]
  r <- tryCatch({ ZWX <- crossprod(Z, W * X); I <- crossprod(Z, W * Z) - ZWX %*% solve(crossprod(X, W * X), t(ZWX))
                  diag(I) / colSums(W * Z^2) }, error = function(e) NA_real_)
  min(r)
}
v_row <- function(x) if (is.null(x) || is_err(x) || !is.data.frame(x$val)) c(S = NA, df = NA, p = NA) else
  c(S = x$val$Test_Statistic[1], df = x$val$df[1], p = x$val$p_value[1])
guard_rep <- function(r, design) {
  ns <- asNamespace("ebrahim.gof")
  g <- v_gen(design, r); fit <- v_fit(g); n <- length(fit$y); ev <- sum(fit$y)
  info <- data.frame(design = design, r = r, n = n, events = ev, slope = unname(coef(fit)[2]), sd_eta = sd(fit$linear.predictors),
                     maxph = max(fitted(fit)))
  rows <- list()
  if (ev == 0 || ev == n) {
    z1 <- cap(def.gof(fit, basis = "stukel", weights = "score")); z0 <- cap(PREG$def.gof(fit, basis = "stukel", weights = "score"))
    zs <- cap(run.all.gof(fit, tests = "Stukel"))
    return(cbind(info, test = "degenerate", ratio = NA, fire_pred = NA, fire_obs = NA, pre_S = v_row(z0)["S"], pre_df = v_row(z0)["df"],
                 pre_p = v_row(z0)["p"], new_S = v_row(z1)["S"], new_df = v_row(z1)["df"], new_p = v_row(z1)["p"],
                 new_noinfo = as.numeric(any(z1$cls == "def_degenerate")),
                 h_p = NA, h_guard = NA, warm_S = NA, note = zs$val$Note[1], row.names = NULL))
  }
  Gr <- max(10, round(n / 25))
  fq <- bt_fit(g)
  for (G in unique(c(10, Gr))) {
    gs <- bt_groups(fq, G)
    for (b in c("poly2", "poly3", "stukel", "sym")) {
      rm(list = ls(V_REC), envir = V_REC)
      z0 <- cap(PREG$def.gof(fit, G = G, basis = b, weights = "score"))
      rc <- if (exists("rec", envir = V_REC)) get("rec", envir = V_REC) else NULL
      z1 <- cap(def.gof(fit, G = G, basis = b, weights = "score"))
      ratio <- if (is.null(rc)) NA_real_ else { k <- rc$dI > 0; if (any(k)) min(rc$dI[k] / rc$zz[k]) else Inf }
      fire_pred <- if (is.null(rc)) NA else any(rc$dI > 0 & rc$dI <= 1e-10 * rc$zz)
      a <- v_row(z0); bb <- v_row(z1)
      same <- identical(is.na(a), is.na(bb)) && all(a[!is.na(a)] == bb[!is.na(bb)])
      Zh <- bt_basis(gs$pbar, if (b == "stukel") "stk" else b)
      hs <- if (is.null(Zh)) list(p = NA_real_, guard = NA) else bt_edge_score(gs, fq$A, Zh)
      ws <- NA_real_
      if (isTRUE(ratio < 1e-6)) { zw <- cap(PREG$def.gof(v_fit(g, TRUE), G = G, basis = b, weights = "score")); ws <- v_row(zw)["S"] }
      rows[[length(rows) + 1]] <- cbind(info, test = sprintf("score.%s.G%d", b, G), ratio = ratio, fire_pred = fire_pred, fire_obs = !same,
        pre_S = a["S"], pre_df = a["df"], pre_p = a["p"], new_S = bb["S"], new_df = bb["df"], new_p = bb["p"],
        new_noinfo = as.numeric(any(z1$cls == "def_no_information")), h_p = hs$p, h_guard = as.numeric(isTRUE(hs$guard)), warm_S = ws,
        note = if (is_err(z1)) paste("error:", z1$val$msg) else "", row.names = NULL)
    }
  }
  ctx <- ns$.gof_context(fit)
  s1 <- ns$gof_stukel(ctx); s0 <- PRES$gof_stukel(ctx)
  ratio <- v_stk_ratio(fit)
  a <- c(S = as.numeric(s0$Statistic), df = as.numeric(s0$df), p = as.numeric(s0$p_value))
  bb <- c(S = as.numeric(s1$Statistic), df = as.numeric(s1$df), p = as.numeric(s1$p_value))
  same <- identical(is.na(a), is.na(bb)) && all(a[!is.na(a)] == bb[!is.na(bb)])
  hb <- bt_stukel(fq)
  ws <- NA_real_
  if (isTRUE(ratio < 1e-6)) ws <- as.numeric(PRES$gof_stukel(ns$.gof_context(v_fit(g, TRUE)))$Statistic)
  rows[[length(rows) + 1]] <- cbind(info, test = "stukel.joint", ratio = ratio, fire_pred = isTRUE(ratio <= 1e-10), fire_obs = !same,
    pre_S = a["S"], pre_df = a["df"], pre_p = a["p"], new_S = bb["S"], new_df = bb["df"], new_p = bb["p"],
    new_noinfo = as.numeric(grepl("no information|no Stukel direction", s1$Note)), h_p = unname(hb["Stk.joint"]),
    h_guard = unname(hb["flag.info_guard"]), warm_S = ws, note = s1$Note, row.names = NULL)
  do.call(rbind, rows)
}

## ---- run ----
cat("pkg280_verify_corr_guard | clone", system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE),
    "| workers", NW, "|", format(t0), "\n")
init(SRC)
cat("installed ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "| pre-guard code = 2e0f6d1 (def.gof instrumented, gof_stukel)\n")

cat("\n==== A. samples with no event or no non-event ====\n")
PA <- partA()
agg <- aggregate(ok ~ check, data = transform(PA, check = sub("^(def.gof|edge.gof) (glm G auto|glm|\\(y,ph,X\\)|\\(y,ph\\) naive) .*$", "\\1 \\2 (all bases, forms)", check)),
                 FUN = function(v) sprintf("%d/%d", sum(v), length(v)))
print(agg, row.names = FALSE)
cat(sprintf("\nTOTAL part A: %d checks, %d failed\n", nrow(PA), sum(!PA$ok)))
if (any(!PA$ok)) { cat("failed checks:\n"); print(PA[!PA$ok, ], row.names = FALSE) }
cat("\nnotes and escaping warnings seen (first sample of each kind):\n")
print(unique(PA[grepl("Stukel row|escaping|control", PA$check), c("check", "detail")]), row.names = FALSE)
write.csv(PA, file.path(OUT, "pkg280_verify_corr_degenerate_checks.csv"), row.names = FALSE)
lap("part A")

cl <- makeCluster(NW)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
invisible(clusterCall(cl, init, SRC))
clusterExport(cl, c("cap", "is_err", "v_gen", "v_fit", "v_stk_ratio", "v_row", "guard_rep"))
PB <- list()
for (des in c("sparse300", "sparse500", "base1000", "sparse100")) {
  PB[[des]] <- do.call(rbind, parLapplyLB(cl, 1:5000, guard_rep, design = des)); lap(des)
  saveRDS(PB, file.path(SCR, "guard_partB.rds"))
}
stopCluster(cl)
B <- do.call(rbind, PB); rownames(B) <- NULL
write.csv(B[B$test == "degenerate" | B$fire_obs %in% TRUE | B$fire_pred %in% TRUE | (!is.na(B$ratio) & B$ratio < 1e-6), ],
          file.path(OUT, "pkg280_verify_corr_guard_cases.csv"), row.names = FALSE)

cat("\n==== B. information guard on null samples (5000 each) ====\n")
cat("fire_obs = installed result differs from 2e0f6d1 (pre-guard); fire_pred = ratio <= 1e-10 on a column with positive information\n")
D <- B[B$test != "degenerate", ]
tabB <- do.call(rbind, lapply(split(D, interaction(D$design, D$test, drop = TRUE, lex.order = TRUE)), function(d) {
  f <- d$fire_obs %in% TRUE; nf <- !f
  hp_ok <- (is.na(d$h_p) & is.na(d$new_p)) | (!is.na(d$h_p) & !is.na(d$new_p) & abs(d$h_p - d$new_p) <= 1e-8)
  data.frame(design = d$design[1], test = d$test[1], samples = nrow(d), fires = sum(f), pred_mismatch = sum(f != (d$fire_pred %in% TRUE)),
             to_NA = sum(f & is.na(d$new_p)), partial = sum(f & !is.na(d$new_p)),
             min_removed_p = if (any(f & !is.na(d$pre_p))) min(d$pre_p[f & !is.na(d$pre_p)]) else NA,
             removed_p_lt_0.2 = sum(f & !is.na(d$pre_p) & d$pre_p < 0.2),
             max_events_fired = if (any(f)) max(d$events[f]) else NA, max_sdeta_fired = if (any(f)) max(d$sd_eta[f]) else NA,
             min_ratio_kept = suppressWarnings(min(d$ratio[nf & is.finite(d$ratio)])),
             kept_ratio_1e10_1e8 = sum(nf & is.finite(d$ratio) & d$ratio > 1e-10 & d$ratio <= 1e-8),
             nofire_not_identical = sum(nf & d$fire_pred %in% TRUE),
             harness_p_mismatch = sum(!hp_ok), harness_guard_mismatch = sum((d$h_guard %in% 1) != f, na.rm = TRUE),
             new_errors = sum(grepl("^error", d$note)))
}))
rownames(tabB) <- NULL
print(tabB, digits = 3)
cat(sprintf("\nTOTAL fires: %d over %d statistic-sample pairs; removed p-values below 0.2: %d; fire predicted by ratio but not observed or vice versa: %d\n",
            sum(D$fire_obs %in% TRUE), nrow(D), sum(D$fire_obs %in% TRUE & !is.na(D$pre_p) & D$pre_p < 0.2),
            sum((D$fire_obs %in% TRUE) != (D$fire_pred %in% TRUE))))
cat(sprintf("harness (bt_*) against installed p-value, |diff| > 1e-8 or NA mismatch: %d; guard flag mismatch: %d\n",
            sum(!((is.na(D$h_p) & is.na(D$new_p)) | (!is.na(D$h_p) & !is.na(D$new_p) & abs(D$h_p - D$new_p) <= 1e-8))),
            sum((D$h_guard %in% 1) != (D$fire_obs %in% TRUE))))

cat("\nfirings (every one):\n")
F <- D[D$fire_obs %in% TRUE, ]
F$warm_rel <- abs(F$warm_S - F$pre_S) / pmax(1e-12, abs(F$pre_S))
print(F[, c("design", "r", "events", "slope", "sd_eta", "test", "ratio", "pre_S", "pre_df", "pre_p", "new_df", "new_p", "warm_S", "warm_rel", "h_p", "h_guard", "note")],
      digits = 4, row.names = FALSE)
N <- D[!(D$fire_obs %in% TRUE) & is.finite(D$ratio) & D$ratio < 1e-6, ]
N$warm_rel <- abs(N$warm_S - N$new_S) / pmax(1e-12, abs(N$new_S))
cat(sprintf("\nkept (no firing) with ratio < 1e-6: %d; their statistic on a warm-restarted fit differs by at most %.2e (relative)\n",
            nrow(N), if (nrow(N)) max(N$warm_rel, na.rm = TRUE) else NA))
if (nrow(N)) print(head(N[order(N$ratio), c("design", "r", "events", "sd_eta", "test", "ratio", "new_S", "new_df", "new_p", "warm_S", "warm_rel")], 25),
                   digits = 4, row.names = FALSE)

cat("\ndegenerate null samples (no event):\n")
Dg <- B[B$test == "degenerate", ]
if (nrow(Dg)) print(aggregate(cbind(samples = 1, installed_NA = is.na(new_p), deg_warning = new_noinfo, pre_p_zero = pre_p %in% 0,
                                    stukel_note = grepl("no maximum-likelihood", note)) ~ design, data = Dg, FUN = sum), row.names = FALSE)
cat("\nevents among samples, by design (share with 1-3 events):\n")
ev <- unique(B[, c("design", "r", "events")])
print(aggregate(events ~ design, data = ev, FUN = function(v) c(mean = mean(v), zero = sum(v == 0), one_to_three = sum(v >= 1 & v <= 3))))
cat("\ndone", format(Sys.time()), sprintf("(%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
