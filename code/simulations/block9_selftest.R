## block9_selftest.R -- does the block 9 code do what the pre-declaration says?
## Contract: paper_EDGE/theory/PREDECLARATION_block9_contamination.md (sha256 99eb4c27...).
## Every check below is a sentence of that document turned into a test. The driver is launched only if this passes.
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

PASS <- 0L; FAIL <- 0L
ok <- function(what, cond) {
  if (isTRUE(cond)) { PASS <<- PASS + 1L } else { FAIL <<- FAIL + 1L; cat("  FAILED: ", what, "\n", sep = "") }
}
sec <- function(s) cat("\n== ", s, " ==\n", sep = "")

C <- b9_cells()

sec("section 1: the cell list")
ok("84 cells", nrow(C) == 84L)
ok("3 truths x 3 corruptions x 4 rates x 2 n = 72 for C1-C3",
   sum(C$corruption %in% c("C1", "C2", "C3")) == 72L)
ok("C4 appears 6 times, at rate 0.005 only",
   sum(C$corruption == "C4") == 6L && all(C$rate[C$corruption == "C4"] == 0.005))
ok("clean appears 6 times, at k = 0", sum(C$corruption == "clean") == 6L &&
   all(C$k[C$corruption == "clean"] == 0L))
ok("the three truths are logistic, probit and cloglog", setequal(unique(C$truth), c("logit", "probit", "cloglog")))
ok("the rate grid is 0.001, 0.002, 0.005, 0.01",
   setequal(unique(C$rate[C$corruption %in% c("C1", "C2", "C3")]), c(0.001, 0.002, 0.005, 0.01)))
ok("the sample sizes are 1000 and 5000", setequal(unique(C$n), c(1000L, 5000L)))
ok("B = 1000 in every cell", all(C$B == 1000L))
ok("cell names are unique", !anyDuplicated(C$cell))
ok("k = round(rate * n)", all(C$k == ifelse(C$corruption == "clean", 0L, as.integer(round(C$rate * C$n)))))
ok("k at n = 1000 is 1, 2, 5, 10", setequal(C$k[C$n == 1000 & C$corruption == "C1"], c(1L, 2L, 5L, 10L)))
ok("k at n = 5000 is 5, 10, 25, 50", setequal(C$k[C$n == 5000 & C$corruption == "C1"], c(5L, 10L, 25L, 50L)))
ok("every corrupted cell names its clean partner", all(C$clean_cell %in% C$cell[C$corruption == "clean"]))

sec("section 1: the seeds")
ok("seed_base = 300000000 + cell_id * 10000", all(C$seed_base == 300000000 + C$cell_id * 10000))
ok("cell_id is 1..84 in the table's own order", identical(C$cell_id, seq_len(84L)))
ok("the seed blocks do not overlap at B = 1000", min(diff(sort(C$seed_base))) > 1000L)
## disjoint from the earlier runs: blocks 0-8 seed from the battery's own roots and from 400000000 (block 9b)
ok("the seeds sit below block 9b's 400000000", max(C$seed_base) + 1000L < 400000000)

sec("section 1: what each corruption does")
## the reference draw: the same seed with no corruption at all
clean_draw <- function(truth, n, seed) { set.seed(seed); b9_gen(truth, n, "clean", 0) }
same_seed  <- function(truth, n, corr, rate, seed) { set.seed(seed); b9_gen(truth, n, corr, rate) }
S <- 12345L

for (truth in c("logit", "probit", "cloglog")) {
  g0 <- clean_draw(truth, 1000L, S)
  ok(paste(truth, "clean: no record is corrupted"), length(g0$corrupt) == 0L)

  ## C1: k covariates multiplied by 4, outcome unchanged from the clean draw
  g1 <- same_seed(truth, 1000L, "C1", 0.005, S)
  ok(paste(truth, "C1: exactly k = 5 records changed"), length(g1$corrupt) == 5L)
  ok(paste(truth, "C1: the changed x are 4x the clean x"),
     isTRUE(all.equal(g1$d$x[g1$corrupt], 4 * g0$d$x[g1$corrupt])))
  ok(paste(truth, "C1: every other x is untouched"),
     isTRUE(all.equal(g1$d$x[-g1$corrupt], g0$d$x[-g1$corrupt])))
  ok(paste(truth, "C1: the outcome is the one drawn at the ORIGINAL x"), identical(g1$d$y, g0$d$y))
  ok(paste(truth, "C1: d is untouched"), identical(g1$d$d, g0$d$d))

  ## C3: k covariates divided by 10, outcome unchanged
  g3 <- same_seed(truth, 1000L, "C3", 0.005, S)
  ok(paste(truth, "C3: the changed x are the clean x / 10"),
     isTRUE(all.equal(g3$d$x[g3$corrupt], g0$d$x[g3$corrupt] / 10)))
  ok(paste(truth, "C3: the outcome is unchanged"), identical(g3$d$y, g0$d$y))
  ok(paste(truth, "C3: a corrupted record's prediction becomes milder, not extreme"),
     all(abs(g3$d$x[g3$corrupt]) <= abs(g0$d$x[g3$corrupt])))

  ## C2: the k highest TRUE risks have y set to 0; the covariates are untouched
  g2 <- same_seed(truth, 1000L, "C2", 0.005, S)
  ok(paste(truth, "C2: x and d are untouched"),
     identical(g2$d$x, g0$d$x) && identical(g2$d$d, g0$d$d))
  ok(paste(truth, "C2: the zeroed records are the k highest true risks"),
     setequal(g2$corrupt, order(g0$p_true, decreasing = TRUE)[1:5]))
  ok(paste(truth, "C2: those outcomes are 0"), all(g2$d$y[g2$corrupt] == 0))
  ok(paste(truth, "C2: no other outcome changed"),
     identical(g2$d$y[-g2$corrupt], g0$d$y[-g2$corrupt]))

  ## C4: the k highest-risk records copied over k randomly chosen rows
  g4 <- same_seed(truth, 1000L, "C4", 0.005, S)
  src <- order(g0$p_true, decreasing = TRUE)[1:5]
  ok(paste(truth, "C4: the destination rows carry the high-risk records"),
     isTRUE(all.equal(sort(g4$d$x[g4$corrupt]), sort(g0$d$x[src]))))
  ok(paste(truth, "C4: the row count is unchanged"), nrow(g4$d) == 1000L)
}

sec("section 1: the truth is the link it names")
for (truth in c("logit", "probit", "cloglog")) {
  g <- clean_draw(truth, 2000L, 99L)
  eta <- 0.6 * g$d$x + 0.5 * g$d$d
  ok(paste(truth, ": p_true is the inverse link at 0.6x + 0.5d"),
     isTRUE(all.equal(g$p_true, BT_LINKINV[[truth]](eta))))
}

sec("reproducibility")
r1 <- b9_one(7L, as.list(C[C$cell == "logit_C1_r005_n1000", ][1, ]))
r2 <- b9_one(7L, as.list(C[C$cell == "logit_C1_r005_n1000", ][1, ]))
ok("the same replicate twice gives the same row", isTRUE(all.equal(r1, r2)))
ok("the row records its own seed", r1[["seed"]] == 300030000 + 7L)
r8 <- b9_one(8L, as.list(C[C$cell == "logit_C1_r005_n1000", ][1, ]))
ok("a different replicate gives a different draw", !isTRUE(all.equal(r1[["b.x"]], r8[["b.x"]])))

sec("section 1: the declared tests are all present")
ce <- as.list(C[C$cell == "logit_C1_r005_n1000", ][1, ])
nm <- names(r1)
need <- c("Stk.joint", "Stk.sym1", "Stk.LR", "GiViTI", "Cubic.LR",
          "EDGE.poly3.u.G10", "EDGE.poly3.sc.G10", "EDGE.poly3.u.Grule", "EDGE.poly3.sc.Grule",
          "EDGE.sym.u.G10", "EDGE.sym.sc.G10", "EDGE.sym.u.Grule", "EDGE.sym.sc.Grule",
          "HL.G10", "HL.Grule", "HLF.G10", "HLF.Grule")
for (t in need) ok(paste("the row carries", t), t %in% nm)
ok("the fitted coefficients are stored", all(c("b.intercept", "b.x", "b.d") %in% nm))
ok("the corrupted-in-top-group count is stored per G arm",
   all(paste0("corrupt_in_top.", c("G10", "Grule")) %in% nm))
ok("k_corrupt matches the cell", r1[["k_corrupt"]] == 5L)
ok("the rule G arm at n = 1000 is G = 40", bt_arms(ce)[["Grule"]] == 40L)
ok("the p-values that ran are in [0, 1]",
   all(vapply(need, function(t) { v <- r1[[t]]; !is.finite(v) || (v >= 0 && v <= 1) }, logical(1))))
ok("the corrupted-in-top count cannot exceed k",
   r1[["corrupt_in_top.G10"]] <= 5L && r1[["corrupt_in_top.Grule"]] <= 5L)

sec("section 2 rule 1: under the logistic truth the bulk model is correct")
## a quick sanity run: with no corruption the false-alarm rate should sit near 0.05, not far from it.
set.seed(1)
cl_cell <- as.list(C[C$cell == "logit_clean_n1000", ][1, ])
pv <- vapply(1:60, function(r) b9_one(r, cl_cell)[["EDGE.poly3.u.Grule"]], numeric(1))
ok("the clean logistic cell does not reject nearly everything (EDGE-poly3 unit)",
   mean(is.finite(pv) & pv <= 0.05) < 0.25)
ok("the clean logistic cell produces p-values at all", mean(is.finite(pv)) > 0.95)

cat(sprintf("\n%d checks, %d failed\n", PASS + FAIL, FAIL))
if (FAIL > 0L) quit(save = "no", status = 1L)
