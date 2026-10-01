## block9d_selftest.R -- does the block 9d code do what its pre-declaration says?
## Contract: paper_EDGE/theory/PREDECLARATION_block9d_robust.md (sha256 799a85aa...).
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
suppressMessages(library(data.table))
SIMDIR <- edge_path("code/simulations")
source(file.path(SIMDIR, "_battery_tests.R"))
source(file.path(SIMDIR, "_block9_contam.R"))
source(file.path(SIMDIR, "_block9d_robust.R"))

PASS <- 0L; FAIL <- 0L
ok  <- function(what, cond) { if (isTRUE(cond)) PASS <<- PASS + 1L else { FAIL <<- FAIL + 1L; cat("  FAILED: ", what, "\n", sep = "") } }
sec <- function(s) cat("\n== ", s, " ==\n", sep = "")

sec("section 1: cells and seeds")
K <- b9d_configs(); X <- b9d_cells()
ok("8 data configurations", nrow(K) == 8L)
ok("16 declared cells", nrow(X) == 16L)
ok("truths are logistic and probit", setequal(unique(K$truth), c("logit", "probit")))
ok("k in {0, 5, 10, 25}", setequal(unique(K$k), c(0L, 5L, 10L, 25L)))
ok("n = 1000 everywhere", all(K$n == 1000L))
ok("each configuration appears under both estimators",
   all(table(X$config_id) == 2L) && setequal(unique(X$estimator), c("ML", "robust")))
ok("seed_base = 600000000 + config_id * 10000", all(K$seed_base == 600000000 + K$config_id * 10000))
ok("seed streams are disjoint at B = 1000", min(diff(sort(K$seed_base))) > 1000L)
ok("seeds sit clear of block 9c (5e8) and below 7e8", min(K$seed_base) > 5e8 + 1e6 && max(K$seed_base) < 7e8)
ok("cell names are unique", !anyDuplicated(X$cell))

sec("the robust fit")
cfg <- as.list(K[K$truth == "logit" & K$k == 10L, ][1, ])
set.seed(cfg$seed_base + 1)
g   <- b9_gen("logit", 1000L, "C1", 0.01)
dat <- list(d = g$d, f = g$f)
fr  <- b9d_fit_robust(dat)
fm  <- bt_fit(dat)
ok("glmrob converges on a contaminated sample", !is.null(fr))
ok("the robust object carries every field bt_fit provides", setequal(names(fr), names(fm)))
ok("the fields have the same shapes",
   length(fr$ph) == length(fm$ph) && identical(dim(fr$X), dim(fm$X)) && identical(dim(fr$A), dim(fm$A)))
ok("the robust estimate differs from ML under contamination",
   max(abs(coef(fr$fit) - coef(fm$fit))) > 1e-4)
ok("the robust slope on x is closer to the truth (0.6) than ML's",
   abs(coef(fr$fit)[2] - 0.6) <= abs(coef(fm$fit)[2] - 0.6) + 0.02)
ok("the fitted link is logit, as the Stukel tests require", identical(fr$fit$family$link, "logit"))

sec("the override is always restored")
env <- environment(battery_rep)
before <- get("bt_fit", envir = env)
cell <- list(type = "full", ao = FALSE, G_extra = "", n = 1000L, generator = "b9d")
invisible(b9d_battery_on(dat, cell, fr))
ok("bt_fit is the original function after a robust call", identical(get("bt_fit", envir = env), before))
invisible(tryCatch(b9d_battery_on(dat, cell, NULL), error = function(e) NULL))
ok("bt_fit is restored even when the call fails", identical(get("bt_fit", envir = env), before))
r_ml_1 <- battery_rep(dat, cell)
invisible(b9d_battery_on(dat, cell, fr))
r_ml_2 <- battery_rep(dat, cell)
ok("an ML call after a robust call gives the identical ML result", isTRUE(all.equal(r_ml_1, r_ml_2)))

sec("one replicate: one data set, two fits, two rows")
res <- b9d_one(1L, cfg)
ok("two rows, named ML and robust", setequal(names(res), c("ML", "robust")))
ok("both rows record the same seed", res$ML[["seed"]] == res$robust[["seed"]])
ok("the ML row equals battery_rep on the same data", {
  set.seed(cfg$seed_base + 1); gg <- b9_gen("logit", 1000L, "C1", 0.01)
  rr <- battery_rep(list(d = gg$d, f = gg$f), cell)
  isTRUE(all.equal(unname(res$ML[names(rr)]), unname(rr)))
})
ok("the robust row reports convergence", res$robust[["robust_converged"]] == 1)
ok("k_corrupt = 10 in both rows", res$ML[["k_corrupt"]] == 10 && res$robust[["k_corrupt"]] == 10)
need <- c("EDGE.poly3.u.Grule", "EDGE.sym.u.Grule", "EDGE.poly3.sc.Grule", "EDGE.sym.sc.Grule",
          "Stk.joint", "Stk.sym1", "Stk.LR", "GiViTI", "Cubic.LR", "HL.Grule", "HLF.Grule",
          "HL.G10", "HLF.G10")
for (t in need) ok(paste("both rows carry", t), t %in% names(res$ML) && t %in% names(res$robust))
ok("the robust row's p-values that ran lie in [0, 1]",
   all(vapply(need, function(t) { v <- res$robust[[t]]; !is.finite(v) || (v >= 0 && v <= 1) }, logical(1))))
ok("no error flag in either row", res$ML[["flag.b9d_error"]] == 0 && res$robust[["flag.b9d_error"]] == 0)
## the failure this test was written after: a robust fit missing a field made battery_rep() stop, and
## the tryCatch turned the WHOLE robust row into NA -- which the [0, 1] check above passes trivially
core <- c("EDGE.poly3.u.Grule", "EDGE.sym.u.Grule", "Stk.joint", "Stk.sym1", "GiViTI", "Cubic.LR",
          "HL.Grule", "HLF.Grule")
ok("the robust row's core tests actually produced p-values (no silent NA row)",
   all(is.finite(res$robust[core])))
ok("the ML row's core tests produced p-values", all(is.finite(res$ML[core])))

sec("reproducibility")
res2 <- b9d_one(1L, cfg)
ok("the same replicate twice gives the same two rows",
   isTRUE(all.equal(res$ML, res2$ML)) && isTRUE(all.equal(res$robust, res2$robust)))

cat(sprintf("\n%d checks, %d failed\n", PASS + FAIL, FAIL))
if (FAIL > 0L) quit(save = "no", status = 1L)
