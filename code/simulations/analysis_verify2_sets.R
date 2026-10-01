## analysis_verify2_sets.R -- three fresh synthetic battery roots (seeds 8810101, 8810202, 8810303) built for the E13 readings,
## with analyse_M_battery.R and the second implementation of analysis_verify2_impl.R run on each at both arms and every output
## compared. Set 1 places sizes at and beyond the nominal limit and engineers verdicts that flip when the size-failing cells are
## dropped (E13.1, E13.3); set 2 removes values and declines matched nulls (E13.2); set 3 builds the H4 census and its Holm
## family (E13.4). No real battery result is read.
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

source(file.path(edge_path("code/simulations"), "analysis_verify2_impl.R"))
options(error = function() {                       # where a failure happened, with the sink closed first
  while (sink.number() > 0) sink()
  cl <- sys.calls(); writeLines("---- call stack ----")
  for (i in seq_along(cl)) writeLines(sprintf("%2d  %s", i, substr(paste(deparse(cl[[i]]), collapse = " "), 1, 160)))
  quit(save = "no", status = 1)
})
sink(file.path(V2, "verify2_sets.log"), split = TRUE)
cat("analysis_verify2_sets.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
shift_mat <- function(v) matrix(v, nrow = 4, ncol = 6, byrow = TRUE, dimnames = list(V2_BASES4, NULL))
set.seed(8810)
most_common_null <- function(d, skip = "") {
  n <- v2_null(d$block, d$cell); k <- paste(n$null_block, n$null_cell)
  tb <- sort(table(k[nzchar(n$null_cell)]), decreasing = TRUE)
  strsplit(setdiff(names(tb), skip)[1], " ")[[1]]
}

## ---- set 1: sizes at the limit, and verdicts that flip without the size-failing cells --------------------------------------------
n500 <- CT[block == "4" & role == "null" & B == 500]$cell[1]
n2k <- CT[block == "1b" & role == "null" & B == 2000]$cell[1:2]
bound2k <- 0.05 + 3 * sqrt(0.05 * 0.95 / 2000)
cat(sprintf("\nset 1 placements: gate %s (B 500, size 0.024 at 0.01) and %s / %s (B 2000, size %.10f and +1e-6 at 0.05)\n",
            n500, n2k[1], n2k[2], bound2k))
post1 <- function(S) {
  f1 <- FAMS[family == 1]; N1 <- most_common_null(f1)
  n1 <- v2_null(f1$block, f1$cell); sh1 <- paste(n1$null_block, n1$null_cell) == paste(N1, collapse = " ")
  x1 <- 0.09 + 0.48 / sum(sh1)                                   # mean 0.105 over the 32 cells, 0.09 without the shared ones
  S <- v2_put(S, f1, "EDGE.poly3.u.Grule", 0.55)
  S <- v2_put(S, f1, "GiViTI", 0.55 - ifelse(sh1, x1, 0.09))
  S <- v2_place_null(S, N1[1], N1[2], "GiViTI", 0.01)            # GiViTI fails size there
  f2 <- FAMS[family == 2]; N2 <- most_common_null(f2, paste(N1, collapse = " "))
  n2 <- v2_null(f2$block, f2$cell); k2 <- paste(n2$null_block, n2$null_cell)
  bad2 <- k2 %in% c(paste(N1, collapse = " "), paste(N2, collapse = " "))
  y2 <- (0.06 * nrow(f2) - 0.04 * (nrow(f2) - sum(bad2))) / sum(bad2)   # mean -0.06 over 21 cells, -0.04 without the failing ones
  S <- v2_put(S, f2, "EDGE.poly3.u.Grule", 0.5)
  S <- v2_put(S, f2, "GiViTI", 0.5 + ifelse(bad2, y2, 0.04))
  S <- v2_place_null(S, N2[1], N2[2], "GiViTI", 0.01)
  nb5 <- v2_null(H5B$block, H5B$cell)
  S <- v2_put(S, H5B, "EDGE.sym.u.Grule", 0.5)
  S <- v2_put(S, H5B, "Stk.sym1", 0.5 + c(0.055, 0.015, 0.015, 0.015))  # cost 0.025 over the four, 0.015 without the first
  S <- v2_place_null(S, nb5$null_block[1], nb5$null_cell[1], "Stk.sym1", 0.01)
  cat(sprintf("  set 1: family 1 null %s/%s carries %d cells (gain %.4f there, 0.09 elsewhere); family 2 null %s/%s, %d failing cells (diff %+.4f)\n",
              N1[1], N1[2], sum(sh1), x1, N2[1], N2[2], sum(bad2), -y2))
  S
}
## the roots carry blocks 1a-4 (every cell of Section A and H1-H5) and per-replicate files for the cells listed here; a pair
## whose cell has no file is compared as "no per-replicate file" by both implementations
## Every fread costs seconds while the battery and another job hold the machine, so the roots carry per-replicate files
## for a few cells only; a pair whose cell has no file is compared as "no per-replicate file" by both implementations.
BLK <- c("1a", "1b", "2", "3", "4")
PF1 <- unique(rbind(FAMS[family == 1][20, .(block, cell)], FAMS[family == 2][1, .(block, cell)], H4CELLS[22]))
PF3 <- unique(H4CELLS[c(1, 4, 22)])          # H4: one not detectable, one saturated, one led census cell
P1 <- list(noise = 0.03, rho = 0.5, nofile = integer(0), dropcol = integer(0), blocks = BLK, reps = 120L, pfile = PF1,
           k = function(tests, ce) ifelse(startsWith(tests, "EDGE."), runif(length(tests), 0.5, 0.85), runif(length(tests), 0.55, 1.05)),
           off = function(tests, ce) 0,
           shift = shift_mat(c(0.05, 0.05, 0.05, 0.05, 0.05, -0.04, rnorm(6, 0, 0.04), 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
                               0.02, 0.02, 0.02, 0.02, 0.02, 0.02)),
           size = rbind(data.table(block = "4", cell = n500, test = "EDGE.poly3.sc.Grule", alpha = 0.01, rejection = 0.024),
                        data.table(block = "1b", cell = n2k[1], test = "EDGE.sym.u.Grule", alpha = 0.05, rejection = bound2k),
                        data.table(block = "1b", cell = n2k[2], test = "EDGE.stk.u.Grule", alpha = 0.05, rejection = bound2k + 1e-6)),
           post = post1)

## ---- set 2: values that are missing, and matched nulls that declined --------------------------------------------------------------
nd1 <- v2_null("3", "loglog_n1000")                 # a family 2 cell: its null declines GiViTI in 95% of replicates
nd2 <- v2_null("4", "crossover_n2000")              # the negative control: 0.9499 keeps the value
cat(sprintf("\nset 2 placements: declined 0.95 for GiViTI in %s/%s; 0.9499 for EDGE.poly3.sc.Grule in %s/%s\n",
            nd1$null_block, nd1$null_cell, nd2$null_block, nd2$null_cell))
P2 <- list(noise = 0.025, rho = 0.7, nofile = c(3L, 11L), dropcol = 7L, blocks = BLK, reps = 120L, pfile = PF1,
           k = function(tests, ce) runif(length(tests), 0.5, 0.85),
           off = function(tests, ce) ifelse(tests == "GiViTI", -0.09, 0),
           shift = shift_mat(c(-0.05, -0.05, -0.05, -0.05, -0.05, 0.12, rnorm(6, 0, 0.04), -0.02, -0.02, -0.02, -0.02, -0.02, -0.02,
                               0.02, 0.02, 0.02, 0.02, 0.02, 0.02)),
           blank = data.table(block = c("2", "3", "3"),
                              cell = c("^probit_auc_n3500$", "^census_corr_n1000$", "^rough_osc4_n1000$"),
                              test = c("^GiViTI$", "^Xie$", "^EDGE\\.stk\\.sc\\.Grule$")),
           declined = rbind(data.table(block = nd1$null_block, cell = nd1$null_cell, test = "GiViTI", value = 0.95),
                            data.table(block = nd2$null_block, cell = nd2$null_cell, test = "EDGE.poly3.sc.Grule", value = 0.9499)))

## ---- set 3: the H4 census, its boundary cells and its Holm family ----------------------------------------------------------------
post3 <- function(S) {
  K <- copy(H4CELLS); idx <- seq_len(nrow(K))
  nodet <- idx[1:3]; sat <- idx[4:6]; bnd <- idx[7]; led <- idx[8:26]; notled <- setdiff(idx, c(nodet, sat, bnd, led))
  vals <- list(nodet = c(0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10),       # nothing reaches 0.15: not detectable
               sat   = c(0.99, 0.99, 0.98, 0.975, 0.97, 0.50, 0.50),      # every compared test at least 0.97: saturated
               bnd   = c(0.536, 0.55, 0.50, 0.50, 0.50, 0.40, 0.40),      # exactly best - 0.014: not led (strictly)
               led   = c(0.60, 0.55, 0.50, 0.50, 0.50, 0.45, 0.45),
               nol   = c(0.40, 0.60, 0.55, 0.50, 0.50, 0.45, 0.45))
  kind <- rep("led", nrow(K)); kind[nodet] <- "nodet"; kind[sat] <- "sat"; kind[bnd] <- "bnd"; kind[notled] <- "nol"
  M <- do.call(rbind, lapply(kind, function(z) vals[[z]]))
  for (j in seq_len(7)) {
    S <- v2_put(S, K, c("EDGE.poly3.u.G10", "HLF.G10", "HL.G10", "HL_w", "PH", "Tsiatis", "Xie")[j], M[, j])
    S <- v2_put(S, K, c("EDGE.poly3.u.Grule", "HLF.Grule", "HL.Grule", "HL_w", "PH", "Tsiatis", "Xie")[j], M[, j])
  }
  cen <- grep("^census_", K$cell)
  tgt <- K[intersect(led, cen)[1]]                                       # a led census cell with a null of its own
  nt <- v2_null(tgt$block, tgt$cell)
  S <- v2_place_null(S, nt$null_block, nt$null_cell, "HLF.G10", 0.01)    # its best partition test fails size: 19 -> 18
  nb2 <- v2_null("3", "link_cloglog_n1000")
  S <- v2_place_null(S, nb2$null_block, nb2$null_cell, "GiViTI", 0)      # a matched null exactly at the nominal limit: holds
  cat(sprintf("  set 3: census %d not detectable, %d saturated, 1 at the lead boundary, %d led, %d not led; size failure in %s/%s (HLF.G10); %s/%s GiViTI exactly at the limit\n",
              length(nodet), length(sat), length(led), length(notled), nt$null_block, nt$null_cell, nb2$null_block, nb2$null_cell))
  S
}
P3 <- list(noise = 0.02, rho = 0.35, nofile = integer(0), dropcol = integer(0), blocks = BLK, reps = 120L, pfile = PF3,
           k = function(tests, ce) runif(length(tests), 0.5, 0.9),
           off = function(tests, ce) 0,
           shift = shift_mat(c(rep(0.03, 6), rnorm(6, 0, 0.03), rep(-0.01, 6), rep(0.015, 6))),
           post = post3)

## ---- build, run and compare ---------------------------------------------------------------------------------------------------------
roots <- file.path(V2, c("set1", "set2", "set3"))
seeds <- c(8810101, 8810202, 8810303)
prms <- list(P1, P2, P3)
out <- list()
for (i in 1:3) {
  cat(sprintf("\n---- set %d (seed %d) ----\n", i, seeds[i]))
  v2_make(roots[i], seeds[i], prms[[i]])
  out[[i]] <- v2_compare(roots[i], paste0("set", i), "rule")
  if (i == 2L) out[[4L]] <- v2_compare(roots[i], paste0("set", i), "10")     # the G = 10 outputs, on the set that has missing values
}

## ---- what the sets exercised (from the second implementation) ------------------------------------------------------------------------
SET <- "coverage"
cat("\n== what the three sets exercised\n")
hy <- rbindlist(lapply(1:3, function(i) cbind(set = i, out[[i]]$mine$hyp[, .(hypothesis, form, statistic, verdict, size_fail_cells,
                                                                             statistic_without, verdict_without, final_verdict)])))
print(hy[form == "unit"], row.names = FALSE)
gt <- rbindlist(lapply(1:3, function(i) cbind(set = i, out[[i]]$mine$gate[v_pass == FALSE, .(test, block, cell, r5 = round(r5, 5), r1 = round(r1, 5))])))
cat("\ngate failures under the nominal standard error:\n"); print(gt, row.names = FALSE)
check("set 1 exercises the E13.3 flip: at least one final verdict is unresolved (size)",
      any(out[[1]]$mine$hyp$final_verdict == "unresolved (size)"),
      paste(out[[1]]$mine$hyp[final_verdict == "unresolved (size)", sprintf("%s/%s", hypothesis, form)], collapse = " "))
check("set 1 exercises the size limit exactly: a size at the nominal limit passes and one 1e-6 above fails",
      out[[1]]$mine$gate[cell == n2k[1] & test == "EDGE.sym.u.Grule"]$ok5 &&
        !out[[1]]$mine$gate[cell == n2k[2] & test == "EDGE.stk.u.Grule"]$ok5)
check("set 2 exercises E13.2: values without a mean, and a declined matched null",
      nrow(out[[2]]$mine$missing) > 0 && any(out[[2]]$mine$hyp$verdict == "no value") &&
        any(grepl("declined", out[[2]]$A$status)),
      sprintf("%d missing pairs, %d rows with no value, %d powers removed by the declined rule",
              nrow(out[[2]]$mine$missing), sum(out[[2]]$mine$hyp$verdict == "no value"), sum(grepl("declined", out[[2]]$A$status))))
check("set 3 exercises E13.4: the H4 Holm family is the detectable, non-saturated cells",
      all(out[[3]]$mine$census[form == "unit" & version == "G10", sum(headline)] == 23),
      sprintf("headline cells %d of %d; led %d", out[[3]]$mine$census[form == "unit" & version == "G10", sum(headline)],
              nrow(H4CELLS), out[[3]]$mine$hyp[hypothesis == "H4" & form == "unit"]$led_or_tied))

v2_write_checks(file.path(V2, "verify2_sets_checks.csv"))
saveRDS(lapply(out, function(z) z$mine), file.path(V2, "verify2_sets_mine.rds"))
sink()
