## analysis_verify2_lib.R -- part A of the independent verification of the E13 fixes: the cell lists typed from the text of
## the pre-declaration (E3, E7, E12.4, E12.5), and a builder of synthetic battery roots. Nothing here reads a real battery
## result: every number comes from the cell table (battery/cells.csv, names, roles, n and B only) and from rnorm/runif.
## Used by analysis_verify2_sets.R and analysis_verify2_edge.R. Outputs go to battery/_review/analysis/verify2.
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
V2 <- edge_battery("_review", "analysis", "verify2")
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)                                            # the battery is using the machine
dir.create(V2, recursive = TRUE, showWarnings = FALSE)

CT <- fread(edge_battery("cells.csv"),
            colClasses = list(character = c("block", "cell", "role", "null_block", "null_cell", "param", "design", "link",
                                            "family", "generator", "G_list")))
CT[is.na(null_block), null_block := ""]; CT[is.na(null_cell), null_cell := ""]
ckey <- function(d) paste(d$block, d$cell)

## ---- the cell lists, typed from the pre-declaration ------------------------------------------------------------------------
N3 <- c(200, 500, 1000, 2000, 5000)
nm <- function(stem, ns) sprintf("%s_n%d", stem, as.integer(ns))
grd <- function(stem, sev, ns) as.vector(outer(paste0(stem, "_", sev), ns, nm))
E7 <- data.table(design = c(rep("base", 7), rep("auc", 8), rep("e12", 8)),
                 link = c("probit", "probit", "cauchit", "cauchit", "t4", "loglog", "loglog",
                          "probit", "probit", "cauchit", "cauchit", "t4", "t4", "loglog", "loglog",
                          "probit", "probit", "cauchit", "cauchit", "t4", "t4", "loglog", "loglog"),
                 n = c(9400, 16000, 2400, 4000, 20000, 610, 1000,
                       3500, 5900, 460, 760, 4700, 7800, 380, 640,
                       4300, 7200, 6400, 11000, 8700, 15000, 1000, 1700))
E7[, cell := sprintf("%s_%s_n%d", link, design, as.integer(n))]
FAMS <- rbind(
  data.table(family = 1L, block = "2", cell = E7[link %in% c("probit", "cauchit", "t4")]$cell),
  data.table(family = 1L, block = "3", cell = c(nm("link_probit", N3), nm("cauchit", N3), nm("t4", N3))),
  data.table(family = 2L, block = "2", cell = E7[link == "loglog"]$cell),
  data.table(family = 2L, block = "3", cell = c(nm("link_cloglog", N3), nm("loglog", N3), nm("stk_asym", N3))),
  data.table(family = 3L, block = "3", cell = c(nm("stk_long", N3), nm("stk_short", N3))),
  data.table(family = 4L, block = "3", cell = c(grd("quad", c("0.01", "0.02", "0.03", "0.05", "0.1", "0.2", "0.4"), N3),
                                                grd("binint", c("0.1", "0.2", "0.3", "0.5", "0.7"), N3),
                                                grd("contint", c("0.1", "0.3", "0.5", "0.7"), N3))),
  data.table(family = 5L, block = "3", cell = c(nm("rough_osc2", c(500, 1000, 2000)), nm("rough_osc4", c(500, 1000, 2000)),
                                                nm("rough_sawtooth", c(500, 1000, 2000)))),
  data.table(family = 5L, block = "4", cell = c(nm("osc4", c(1000, 2000)), nm("sawtooth", c(1000, 2000)))),
  data.table(family = 6L, block = "4", cell = c("crossover_n1000", "crossover_n2000")))
FAMS[, design := "base"]
FAMS[block == "2", design := E7$design[match(cell, E7$cell)]]
FAMS[, n := as.integer(sub("^.*_n([0-9]+)$", "\\1", cell))]
CENSUS <- nm(paste0("census_", c("logx", "int_binbin", "skew", "corr", "omit_x2x3", "omit_2int", "joint", "omit_2cov")), 1000)
H4CELLS <- rbind(FAMS[family %in% 2:4 & block == "3" & n == 1000, .(block, cell)], data.table(block = "3", cell = CENSUS))
H5B <- data.table(block = "2", cell = c("probit_base_n9400", "probit_base_n16000", "cauchit_base_n2400", "cauchit_base_n4000"))
H5S <- data.table(block = "4", cell = c("probit_skew_n9400", "probit_skew_n16000", "cauchit_skew_n2400", "cauchit_skew_n4000"))
## the matched null of a cell, from the cell table
v2_null <- function(block, cell) {
  m <- match(paste(block, cell), ckey(CT))
  data.table(null_block = ifelse(is.na(m), "", CT$null_block[m]), null_cell = ifelse(is.na(m), "", CT$null_cell[m]))
}
v2_B <- function(block, cell) CT$B[match(paste(block, cell), ckey(CT))]

## ---- the tests of a replicate ------------------------------------------------------------------------------------------------
V2_BASES4 <- c("poly3", "poly2", "stk", "sym")
V2_EDGE <- as.vector(outer(as.vector(outer(paste0("EDGE.", V2_BASES4), c("u", "sc"), paste, sep = ".")), c("G10", "Grule"), paste, sep = "."))
V2_OTHER <- c("HL.G10", "HL.Grule", "HLF.G10", "HLF.Grule", "Stk.joint", "Stk.LR", "Stk.sym1", "Stk.marg", "GiViTI", "GiViTI.t50",
              "Cubic.LR", "HL_w", "PH", "Tsiatis", "Xie", "PR")
V2_RIV <- c("HL_w", "PH", "Tsiatis", "Xie", "PR")
V2_PCOLS <- c(paste0("EDGE.", rep(c("poly3", "sym"), each = 4), ".", c("u", "sc"), ".", rep(c("G10", "Grule"), each = 2)),
              "GiViTI", "Stk.joint", "Stk.sym1", "HLF.G10", "HL.G10", "HLF.Grule", "HL.Grule", "HL_w", "PH", "Tsiatis", "Xie")
## the size with z exactly against the null standard error, and against the realised one
v2_size_null_z <- function(z, a, B) a + z * sqrt(a * (1 - a) / B)
v2_size_real_z <- function(z, a, B) { k <- z^2 / B; (2 * a + k + sqrt((2 * a + k)^2 - 4 * (1 + k) * a^2)) / (2 * (1 + k)) }

## ---- a synthetic root ----------------------------------------------------------------------------------------------------------
## prm: noise, rho, k(tests, cell) size multiplier of a null, off(tests, cell) power offset, shift (4 bases x 6 families, score
## minus unit), size (exact null sizes), declined (exact declined rates of a null cell-test, every alpha row), blank (cells and
## tests whose size-adjusted power is removed, as run_M_battery.R writes an alternative whose null has not run), post(S) a free
## hand on the summary, nofile / dropcol (per-replicate files that are missing or lack a column), pfile (cells that get one)
v2_make <- function(dir, seed, prm) {
  set.seed(seed)
  unlink(dir, recursive = TRUE); dir.create(dir, recursive = TRUE)
  file.copy(edge_battery("cells.csv"), file.path(dir, "cells.csv"))
  blocks <- if (is.null(prm[["blocks"]])) c("1a", "1b", "2", "3", "4", "5", "6", "7") else prm[["blocks"]]
  C <- merge(CT[block %in% blocks & role %in% c("null", "alternative")],
             FAMS[, .(block, cell, fam = family, fdesign = design)], by = c("block", "cell"), all.x = TRUE, sort = FALSE)
  tests0 <- c(V2_EDGE, V2_OTHER)
  L <- vector("list", nrow(C))
  for (i in seq_len(nrow(C))) {
    ce <- C[i]; Bc <- ce$B
    tests <- if (ce$generator == "external") setdiff(tests0, grep("\\.G(10|rule)$", tests0, value = TRUE)) else tests0
    notrun <- ce$n >= 10000 & tests %in% V2_RIV
    if (ce$role == "null") {
      k <- prm[["k"]](tests, ce)
      R <- rbindlist(lapply(c(0.01, 0.05, 0.10), function(a) {
        r <- rbinom(length(tests), Bc, pmin(1, a * k)) / Bc
        r[notrun] <- 0
        data.table(test = tests, alpha = a, rejection = r, size_adj_power = NA_real_, null_size = NA_real_,
                   status = ifelse(notrun, "not run (n >= 10,000)", ""))
      }))
    } else {
      pw <- runif(1, 0.05, 0.9) + rnorm(length(tests), 0, prm[["noise"]]) + prm[["off"]](tests, ce)
      e <- startsWith(tests, "EDGE.") & grepl("\\.sc\\.", tests)
      if (!is.na(ce$fam) && !is.null(prm[["shift"]]))
        pw[e] <- pw[e] + prm[["shift"]][cbind(match(sub("^EDGE\\.([a-z0-9]+)\\..*$", "\\1", tests[e]), V2_BASES4), ce$fam)]
      pw <- pmin(0.995, pmax(0.005, pw))
      R <- rbindlist(lapply(c(0.01, 0.05, 0.10), function(a) {
        s <- pmin(1, pw * c(0.55, 1, 1.25)[match(a, c(0.01, 0.05, 0.10))]); s[notrun] <- 0
        data.table(test = tests, alpha = a, rejection = pmin(1, s + 0.004), size_adj_power = s, null_size = a * 0.97,
                   status = ifelse(notrun, "not run (n >= 10,000)", ""))
      }))
    }
    R[, `:=`(block = ce$block, cell = ce$cell, role = ce$role, null_cell = ce$null_cell, n = ce$n, subset = "all", B = Bc)]
    R[, mcse := sqrt(rejection * (1 - rejection) / B)]
    R[, declined := round(runif(.N, 0, 0.03), 4)]
    R[, rejection_given_p := pmin(1, rejection / (1 - declined))]
    if (i %% 4L == 0L) {                      # decoy subsets: far-off numbers, so reading them instead of subset = all would show
      D1 <- copy(R)[, `:=`(subset = "events<Grule", B = as.integer(round(Bc / 3)), rejection = 0.5, mcse = 0.001, size_adj_power = 0.001)]
      D2 <- copy(R)[, `:=`(subset = "events>=Grule", B = Bc - as.integer(round(Bc / 3)), rejection = 0.5, mcse = 0.001, size_adj_power = 0.999)]
      R <- rbind(R, D1, D2)
    }
    FL <- data.table(test = c("flag.evlt.G10", "flag.evlt.Grule", "flag.degenerate"), alpha = NA_real_,
                     rejection = c(if (i %% 7L == 0L) 1 else 0, if (i %% 4L == 0L) 0.3 else 0, 0),
                     size_adj_power = NA_real_, null_size = NA_real_, status = "flag rate: share of replicates with the flag = 1",
                     block = ce$block, cell = ce$cell, role = ce$role, null_cell = ce$null_cell, n = ce$n, subset = "all", B = Bc,
                     mcse = NA_real_, declined = 0, rejection_given_p = NA_real_)
    L[[i]] <- rbind(FL, R, use.names = TRUE)
  }
  S <- rbindlist(L, use.names = TRUE)
  hit1 <- function(z, extra = TRUE) {
    h <- which(S$block == z$block & S$cell == z$cell & S$test == z$test & S$subset == "all" & extra)
    if (!length(h)) stop("override not found: ", z$block, " ", z$cell, " ", z$test)
    h
  }
  SZ <- prm[["size"]]; DE <- prm[["declined"]]; BL <- prm[["blank"]]; DR <- prm[["dropcells"]]      # exact names: $ matches partly
  if (!is.null(SZ)) for (j in seq_len(nrow(SZ))) {                                 # an exact null size at one level
    z <- SZ[j]; h <- hit1(z, !is.na(S$alpha) & abs(S$alpha - z$alpha) < 1e-9)
    set(S, h, c("rejection", "mcse"), list(z$rejection, sqrt(z$rejection * (1 - z$rejection) / S$B[h])))
  }
  if (!is.null(DE)) for (j in seq_len(nrow(DE))) {                                 # a declined rate on every alpha row of a null test
    z <- DE[j]; h <- hit1(z, !is.na(S$alpha))
    set(S, h, "declined", z$value)
    set(S, h, "rejection_given_p", pmin(1, S$rejection[h] / (1 - z$value)))
  }
  if (!is.null(BL)) for (j in seq_len(nrow(BL))) {                                 # no size-adjusted power
    z <- BL[j]
    h <- which(S$block == z$block & grepl(z$cell, S$cell) & grepl(z$test, S$test) & !startsWith(S$test, "flag."))
    if (!length(h)) stop("blank not found: ", z$cell, " ", z$test)
    set(S, h, c("size_adj_power", "status"), list(NA_real_, "matched null not run yet"))
  }
  if (!is.null(DR)) S <- S[!paste(block, cell) %in% paste(DR$block, DR$cell)]      # whole cells with no summary row at all
  if (!is.null(prm[["post"]])) S <- prm[["post"]](S)
  setcolorder(S, c("block", "cell", "role", "null_cell", "n", "subset", "B", "test", "alpha", "rejection", "mcse", "size_adj_power",
                   "null_size", "declined", "rejection_given_p", "status"))
  for (b in unique(S$block)) { dir.create(file.path(dir, b), showWarnings = FALSE); fwrite(S[block == b], file.path(dir, b, "_summary.csv")) }
  ## per-replicate p-values for the cells of H1-H5, drawn around each test's size-adjusted power and correlated within a cell
  PC <- if (is.null(prm[["pfile"]])) unique(rbind(FAMS[family %in% 1:2, .(block, cell)], H4CELLS, H5B, H5S)) else prm[["pfile"]]
  W5 <- S[subset == "all" & !is.na(alpha) & abs(alpha - 0.05) < 1e-9]
  wk <- paste(W5$block, W5$cell, W5$test)
  Br <- if (is.null(prm[["reps"]])) 300L else as.integer(prm[["reps"]])
  for (j in seq_len(nrow(PC))) {
    if (j %in% prm[["nofile"]]) next
    bj <- PC$block[j]; cj <- PC$cell[j]
    z <- rnorm(Br)
    P <- data.table(rep = seq_len(Br), seed = 1e6 + seq_len(Br), n = CT$n[match(paste(bj, cj), ckey(CT))], events = 100L)
    for (tt in V2_PCOLS) {
      pw <- W5$size_adj_power[match(paste(bj, cj, tt), wk)]
      if (!length(pw) || !is.finite(pw)) pw <- 0.3
      pw <- min(0.995, max(0.005, pw))
      p <- pnorm(qnorm(0.05) - qnorm(pw) + prm[["rho"]] * z + sqrt(1 - prm[["rho"]]^2) * rnorm(Br))
      p[sample(Br, 2)] <- 0.05                                    # p = alpha exactly: a rejection (p <= alpha)
      p[runif(Br) < 0.03] <- NA                                   # declined = no rejection
      set(P, j = tt, value = p)
    }
    if (j %in% prm[["dropcol"]]) P[, Stk.joint := NULL]
    ## plain text under the name the analysis expects: fwrite compresses by extension, and on this machine a gzipped
    ## read costs tens of seconds, while fread reads the plain file to the same values (analysis_verify2_probe.R)
    tmp <- file.path(dir, bj, paste0(cj, "_pvalues.csv"))
    fwrite(P, tmp)
    file.rename(tmp, file.path(dir, bj, paste0(cj, "_pvalues.csv.gz")))
  }
  invisible(S)
}
v2_finished <- function(dir)
  writeLines(c("2026-09-14 08:28:48  launch from step 1 of 19", "2026-09-15 02:00:00  launch finished"), file.path(dir, "launch.log"))

## ---- placing values by hand in a summary (used by the set and edge drivers) ---------------------------------------------------
## one value (or one per cell) of one test at one level, subset = all
v2_put <- function(S, d, test, value, col = "size_adj_power", a = 0.05) {
  h <- match(paste(d$block, d$cell, test, "all", sprintf("%.3f", a)),
             paste(S$block, S$cell, S$test, S$subset, sprintf("%.3f", S$alpha)))
  if (anyNA(h)) stop("v2_put: no row for ", paste(d$cell[is.na(h)], collapse = " "), " ", test[1])
  set(S, h, col, if (length(value) == 1L) rep(value, length(h)) else value)
  S
}
## the size of one test in one null cell at 0.05, placed relative to the nominal limit 0.05 + 3 sqrt(0.05 x 0.95 / B)
v2_place_null <- function(S, nb, nc, test, over = 0.01) {
  h <- which(S$block == nb & S$cell == nc & S$test == test & S$subset == "all" & !is.na(S$alpha) & abs(S$alpha - 0.05) < 1e-9)
  if (length(h) != 1L) stop("v2_place_null: ", nb, " ", nc, " ", test, " matched ", length(h), " rows")
  s <- 0.05 + 3 * sqrt(0.05 * 0.95 / S$B[h]) + over
  set(S, h, "rejection", s); set(S, h, "mcse", sqrt(s * (1 - s) / S$B[h]))
  S
}
