## block8R_selftest.R -- checks run_M_block8R.R before block 8R runs (declaration sha256 586d8f51...). Serial, one core.
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
SIMDIR <- edge_path("code/simulations")
options(block8R.source_only = TRUE)
suppressMessages(source(file.path(SIMDIR, "run_M_block8R.R")))
suppressPackageStartupMessages(library(BAGofT))

res <- list()
chk <- function(id, ok, msg) {
  ok <- isTRUE(ok)
  res[[length(res) + 1]] <<- data.frame(id = id, ok = ok, msg = msg, stringsAsFactors = FALSE)
  cat(sprintf("[%s] %-6s %s\n", if (ok) "PASS" else "FAIL", id, msg))
}

## ---- the plan (section 1) ---------------------------------------------------------------------------------------------
P <- r8_plan()
cnt <- table(P$test, P$role)
chk("R8.1", cnt["lecessie", "alternative"] == 56 && cnt["lecessie", "null"] == 11 && cnt["proj", "alternative"] == 19 &&
            cnt["proj", "null"] == 4 && cnt["bagoft", "alternative"] == 6 && cnt["bagoft", "null"] == 4,
    "plan: le Cessie 56 + 11, projection 19 + 4 at n = 500, BAGofT 6 + 4 at n = 500")
chk("R8.2", all(P$R[P$test == "bagoft"] == 200L) && all(P$R[P$test != "bagoft"] == 500L), "replicates: 500, 500, 200")
b <- P[P$test == "bagoft" & P$role == "alternative", "cell"]
chk("R8.3", setequal(b, R8_BAG_ALT) && all(grepl("_n500$", P$cell[P$test != "lecessie"])),
    sprintf("BAGofT alternatives are the three rough shapes and the median rungs: %s", paste(sort(b), collapse = ", ")))
chk("R8.4", all(P$null_cell[P$role == "alternative"] %in% P$cell[P$role == "null"]),
    "every alternative's matched null is in the plan for the same test")

## ---- the identity gate on every generator type ----------------------------------------------------------------------------
TEST8 <- TRUE; REPS8 <- 2L
gate <- c("binint_0.3_n200", "contint_0.5_n500", "quad_0.05_n500", "rough_osc4_n500", "osc4_n1000", "sawtooth_n1000",
          "null_quad_n500", "null_rough_n1000", "null_binint_n200", "null_contint_n1000")
G <- do.call(rbind, lapply(gate, function(cn) {
  ce <- as.list(P[P$test == "lecessie" & P$cell == cn, ][1, ]); cell <- r8_cell(ce); S <- rv_stored(ce, 2L)
  do.call(rbind, lapply(1:2, function(r) r8_one(r, cell, "lecessie", as.list(S[r, ]))))
}))
chk("R8.5", all(G[, "id_ok"] == 1), sprintf("identity gate passes on %d regenerated data sets of %d scenarios (every generator type)", nrow(G), length(gate)))
chk("R8.6", max(G[, "id.absdiff"], na.rm = TRUE) <= 1e-8, sprintf("largest |EDGE p regenerated - stored| = %.2e", max(G[, "id.absdiff"], na.rm = TRUE)))
chk("R8.7", all(is.finite(G[, "lecessie"])), "le Cessie gives a p-value on every one of them, one-covariate models included")

## ---- the Liu projection test on a one-covariate model ------------------------------------------------------------------
ce <- as.list(P[P$test == "lecessie" & P$cell == "quad_0.05_n200", ][1, ]); cell <- r8_cell(ce); S <- rv_stored(ce, 1L)
t0 <- proc.time()[["elapsed"]]
z <- r8_one(1L, cell, "proj", as.list(S[1, ]))
chk("R8.8", z[["id_ok"]] == 1 && is.finite(z[["proj"]]) && z[["flag.error"]] == 0,
    sprintf("the projection test runs on y ~ x (n = 200): p = %.3f, %.1f s", z[["proj"]], proc.time()[["elapsed"]] - t0))

## ---- BAGofT: the package fails on one covariate, the patch does not, and changes nothing with two ---------------------------
chk("R8.9", r8_patch_check(), "the patch changes exactly one line of parRF, adding drop = FALSE")
PARRF_FIX <- r8_make_parRF_dropfix()
g1 <- rv_regen(1L, cell)                                             # quad_0.05_n200: y ~ x
tiny <- function(dat, f, pf) tryCatch(suppressWarnings(suppressMessages(BAGofT::BAGofT(
  testModel = BAGofT::testGlmBi(formula = f, link = "logit"), data = dat, parFun = pf, nsplits = 3, nsim = 3))),
  error = function(e) e)
set.seed(11); u <- tiny(g1$dat$d, g1$dat$f, BAGofT::parRF())
set.seed(11); v <- tiny(g1$dat$d, g1$dat$f, PARRF_FIX())
chk("R8.10", inherits(u, "error"), sprintf("unpatched BAGofT fails on a one-covariate model: %s", if (inherits(u, "error")) conditionMessage(u) else "it did NOT fail"))
chk("R8.11", !inherits(v, "error") && is.finite(v$p.value), sprintf("patched BAGofT runs on it: p = %s", if (inherits(v, "error")) "error" else format(v$p.value)))
ce2 <- as.list(P[P$test == "lecessie" & P$cell == "binint_0.3_n200", ][1, ]); g2 <- rv_regen(1L, r8_cell(ce2))   # y ~ x + d
same <- TRUE
for (s in 1:3) {
  set.seed(100 + s); a <- tiny(g2$dat$d, g2$dat$f, BAGofT::parRF())
  set.seed(100 + s); b2 <- tiny(g2$dat$d, g2$dat$f, PARRF_FIX())
  same <- same && !inherits(a, "error") && identical(a$p.value, b2$p.value) && identical(a$p.value2, b2$p.value2) &&
          identical(a$p.value3, b2$p.value3) && identical(a$pmean, b2$pmean)
}
chk("R8.12", same, "with two covariates the patched and unpatched BAGofT return bit-identical results (3 seeds)")

## ---- output guards --------------------------------------------------------------------------------------------------------
why <- function(o) tryCatch({ r8_setup(o); "accepted" }, error = function(e) conditionMessage(e))
w <- c(b8L = why(list(out = "battery/8L")), b8 = why(list(out = "battery/8")), t8R = why(list(reps = "2", out = "battery/8R")))
chk("R8.13", grepl("another block", w[["b8L"]]) && grepl("another block", w[["b8"]]) && grepl("real block 8R", w[["t8R"]]),
    "refuses blocks 8 and 8L, and a test run into the real 8R folder, each for its own reason")
chk("R8.14", why(list(reps = "2")) == "accepted" && endsWith(OUT8, "/battery/8R_test"), "a test run goes to battery/8R_test")
chk("R8.15", why(list()) == "accepted" && endsWith(OUT8, "/battery/8R"), "the real run goes to battery/8R")

## ---- one BAGofT replicate at the package defaults, the declared call end to end (one core, several minutes) -------------------
TEST8 <- TRUE; REPS8 <- 1L
ce3 <- as.list(P[P$test == "bagoft" & P$cell == "quad_0.05_n500", ][1, ]); S3 <- rv_stored(ce3, 1L)
t0 <- proc.time()[["elapsed"]]
z3 <- r8_one(1L, r8_cell(ce3), "bagoft", as.list(S3[1, ]))
chk("R8.16", z3[["id_ok"]] == 1 && is.finite(z3[["BAGofT"]]) && z3[["flag.error"]] == 0,
    sprintf("BAGofT at its defaults, patched, on quad_0.05_n500 replicate 1: p = %s, %.0f s", format(z3[["BAGofT"]]), proc.time()[["elapsed"]] - t0))

R <- do.call(rbind, res)
cat(sprintf("\n%d of %d checks pass\n", sum(R$ok), nrow(R)))
data.table::fwrite(R, file.path(SIMDIR, "block8R_selftest.csv"))
if (!all(R$ok)) quit(save = "no", status = 1L)
