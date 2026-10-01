## block8L_selftest.R -- checks run_M_block8L.R before block 8L runs (declaration sha256 0c1a1731...). Serial, a few
## seconds; it computes nothing that block 8L keeps.
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
options(block8L.source_only = TRUE)
suppressMessages(source(file.path(SIMDIR, "run_M_block8L.R")))

res <- list()
chk <- function(id, ok, msg) {
  ok <- isTRUE(ok)
  res[[length(res) + 1]] <<- data.frame(id = id, ok = ok, msg = msg, stringsAsFactors = FALSE)
  cat(sprintf("[%s] %-6s %s\n", if (ok) "PASS" else "FAIL", id, msg))
}

## ---- the implementation (section 1) -------------------------------------------------------------------------------
v <- utils::packageVersion("ebrahim.gof")
chk("L8.1", v >= L8_MINVER, sprintf("ebrahim.gof %s installed, section 1 needs >= %s", v, L8_MINVER))

tc <- l8_transpose_check()
chk("L8.2", abs(tc[["package"]] - tc[["declared"]]) < 1e-10,
    sprintf("installed le Cessie = (I-H)'R(I-H) from its definition: %.12f vs %.12f", tc[["package"]], tc[["declared"]]))
chk("L8.3", abs(tc[["package"]] - tc[["vendor"]]) > 1e-6,
    sprintf("and not the vendor's (I-H)R(I-H): %.12f vs %.12f (the check can tell them apart)", tc[["package"]], tc[["vendor"]]))

## the public route is the declared function, and it draws no random number
set.seed(1); x <- runif(400, -3, 3); d <- rbinom(400, 1, 0.5); y <- rbinom(400, 1, plogis(0.6 * x + 0.5 * d))
fit <- glm(y ~ x + d, family = binomial())
seed_before <- .Random.seed
lc <- l8_lecessie(fit)
chk("L8.4", identical(seed_before, .Random.seed), "le Cessie draws no random number (.Random.seed unchanged)")
int <- ebrahim.gof:::gof_lecessie(ebrahim.gof:::.gof_context(fit))
chk("L8.5", isTRUE(all.equal(unname(lc), c(int$p_value, int$Statistic, int$df), tolerance = 0)),
    sprintf("run.all.gof(tests = 'le-Cessie') is the package's gof_lecessie: p %.10f, stat %.4f, df %.4f", lc[["p"]], lc[["stat"]], lc[["df"]]))
chk("L8.6", identical(l8_lecessie(fit), lc), "deterministic: the same fit gives the same p-value twice")

## ---- the cells (section 2) --------------------------------------------------------------------------------------------
C <- l8_cells()
B <- C[C$part == "B", ]
chk("L8.7", nrow(C) == 35L && sum(C$part == "A") == 30L && all(C$R[C$part == "A"] == 500L) && all(B$R == 1000L),
    "30 Part A cells at replicates 1-500, 5 Part B cells at 1-1000")
chk("L8.8", identical(sort(B$cell), sort(c("logit_clean_n1000", sprintf("logit_C1_r%s_n1000", c("001", "002", "005", "01"))))) ||
            identical(sort(B$cell), sort(c("logit_clean_n1000", sprintf("logit_C1_r%s_n1000", c("001", "002", "005", "010"))))),
    sprintf("Part B = the logistic clean cell and C1 at k = 1, 2, 5, 10 (n = 1000): %s", paste(B$cell, collapse = ", ")))
chk("L8.9", identical(as.integer(B$k[order(B$k)]), c(0L, 1L, 2L, 5L, 10L)), "Part B k = 0, 1, 2, 5, 10")
chk("L8.10", all(Cells8$cell %in% C$cell[C$part == "A"]) && all(C$null_cell[C$part == "A" & !is.na(C$null_cell)] %in% C$cell),
    "Part A = block 8's own cell list, every alternative's matched null inside it")

## ---- the identity gate (section 2, "Data") ----------------------------------------------------------------------------
TEST_RUN <- TRUE; REPS_OVERRIDE <- 3L
rows <- list()
gate_cells <- c("null_link_n500", "stk_long_n1000", "null_crossover_n1000", "crossover_n1000", "loglog_auc_n380",
                "logit_clean_n1000", "logit_C1_r010_n1000")
for (cn in gate_cells) {
  ce <- as.list(C[C$cell == cn, ]); cell <- l8_cell(ce); S <- rv_stored(ce, 3L)
  for (r in 1:3) rows[[length(rows) + 1]] <- c(cell = cn, as.list(l8_one(r, cell, as.list(S[r, ]))))
}
G <- rbindlist(rows)
chk("L8.11", all(G$id_ok == 1), sprintf("identity gate passes on %d regenerated data sets (%d cells x 3), both parts", nrow(G), length(gate_cells)))
chk("L8.12", max(G$id.absdiff, na.rm = TRUE) <= 1e-8, sprintf("largest |EDGE p regenerated - stored| = %.2e", max(G$id.absdiff, na.rm = TRUE)))
chk("L8.13", all(is.finite(G$lecessie)) && all(G$flag.error == 0), "every kept data set has a le Cessie p-value, no error flag")
chk("L8.14", all(G$id_ok != 1 | G$flag.degenerate == 1 | is.finite(G$lecessie) | G$flag.error == 1),
    "no silent NA row: a kept, non-degenerate row has a p-value or an error flag")
g10 <- l8_regen_b9(1L, l8_cell(as.list(C[C$cell == "logit_C1_r010_n1000", ])))
chk("L8.15", g10$k_corrupt == 10L && sum(abs(g10$dat$d$x) > 3) >= 1,
    sprintf("Part B regeneration corrupts k = %d records (%d with |x| > 3)", g10$k_corrupt, sum(abs(g10$dat$d$x) > 3)))

## the gate refuses a data set that does not match, and then does not run le Cessie
ce <- as.list(C[C$cell == "stk_long_n1000", ]); cell <- l8_cell(ce); S <- rv_stored(ce, 1L)
bad_p <- as.list(S[1, ]); bad_p$p <- bad_p$p + 1e-6
z1 <- l8_one(1L, cell, bad_p)
chk("L8.16", z1[["id_ok"]] == 0 && !is.finite(z1[["lecessie"]]), "a stored p-value off by 1e-6 fails the gate; le Cessie not run")
bad_s <- as.list(S[1, ]); bad_s$seed <- bad_s$seed + 1
z2 <- l8_one(1L, cell, bad_s)
chk("L8.17", z2[["id_ok"]] == 0 && !is.finite(z2[["lecessie"]]), "a stored seed that differs fails the gate")
bad_e <- as.list(S[1, ]); bad_e$events <- bad_e$events + 1
z3 <- l8_one(1L, cell, bad_e)
chk("L8.18", z3[["id_ok"]] == 0, "stored events that differ fail the gate")

## ---- rules (section 3) ------------------------------------------------------------------------------------------------
b5 <- l8_band(500); b10 <- l8_band(1000)
chk("L8.19", identical(sprintf("%.3f", b5), c("0.021", "0.079")) && identical(sprintf("%.3f", b10), c("0.029", "0.071")),
    sprintf("size bands as declared: B = 500 %.4f-%.4f, B = 1000 %.4f-%.4f", b5[1], b5[2], b10[1], b10[2]))
chk("L8.20", identical(l8_rej(c(0.05, 0.0500001, NA, 0.01)), c(TRUE, FALSE, FALSE, TRUE)), "rejection at p <= 0.05; no p-value = no rejection")
pn <- c(seq(0.01, 1, by = 0.01), NA)
chk("L8.21", isTRUE(all.equal(l8_crit(pn), stats::quantile(c(seq(0.01, 1, by = 0.01), 1), 0.05, type = 1)[[1]])),
    sprintf("size-adjusting critical value: 5%% quantile, type 1, no p-value counted as 1 (%.2f)", l8_crit(pn)))

## ---- output guards ----------------------------------------------------------------------------------------------------
## each refusal must be for its own reason (a refusal for any other reason once hid a bug here), and a valid setup must pass
why <- function(o) tryCatch({ l8_setup(o); "accepted" }, error = function(e) conditionMessage(e))
w <- c(b9 = why(list(out = "battery/9")), b8 = why(list(out = "battery/8")), b9b = why(list(out = "battery/9b")),
       t8L = why(list(reps = "3", out = "battery/8L")), outside = why(list(out = file.path(tempdir(), "x"))))
chk("L8.22", grepl("folder of another block", w[["b9"]]) && grepl("folder of another block", w[["b8"]]) &&
             grepl("folder of another block", w[["b9b"]]) && grepl("cannot write to the real block 8L", w[["t8L"]]) &&
             grepl("must lie inside", w[["outside"]]),
    sprintf("refuses blocks 8, 9, 9b, a test run into 8L, and outside battery/, each for its own reason: %s",
            paste(names(w), substr(w, 1, 28), sep = " -> ", collapse = " | ")))
ok_test <- why(list(reps = "3"))
chk("L8.23", ok_test == "accepted" && endsWith(OUT, "/battery/8L_test") && TEST_RUN && REPS_OVERRIDE == 3L,
    sprintf("a test run is accepted and goes to %s", OUT))
ok_real <- why(list())
chk("L8.24", ok_real == "accepted" && endsWith(OUT, "/battery/8L") && !TEST_RUN, sprintf("the real run is accepted and goes to %s", OUT))

## ---- cost (L4 reports it) ---------------------------------------------------------------------------------------------
tt <- G$sec[grepl("n1000$", G$cell)]
cat(sprintf("\nle Cessie at n = 1000: median %.2f s per data set (%d data sets)\n", stats::median(tt), length(tt)))

R <- rbindlist(res)
cat(sprintf("\n%d of %d checks pass\n", sum(R$ok), nrow(R)))
fwrite(R, file.path(SIMDIR, "block8L_selftest.csv"))
if (!all(R$ok)) quit(save = "no", status = 1L)
