## battery_review_rivals_regen.R -- independent review of block 8 (E9): the cell list and the data regeneration.
## No rival is computed. One R process, BLAS on one thread. Reads the battery files; writes only under battery/_review/rivals.
##   (a) the 30 cells from my own reading of E9 and E1 against run_M_rivals.R's Cells8 and battery_cells()
##   (b) replicates 1-3 of every dry-run cell of block 8, drawn with my own generators (not bt_data), fitted with glm, and
##       EDGE-poly3 unit at G = 10 from the installed ebrahim.gof::def.gof (not the harness); compared with the driver's
##       rv_regen() data and with the stored EDGE.poly3.u.G10
##   (c) the real battery files of block 8 cells that exist now: seed = seed_base + rep, n, replicates 1-500, and the
##       driver's identity gate on replicates 1, 200 and 500
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
OUTD <- edge_battery("_review", "rivals")
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) { RhpcBLASctl::blas_set_num_threads(1); RhpcBLASctl::omp_set_num_threads(1) }
suppressPackageStartupMessages({ library(data.table); library(ebrahim.gof) })
options(rivals.source_only = TRUE)
source(file.path(SIMDIR, "run_M_rivals.R"))
LOG <- file.path(OUTD, "review_regen.log")
if (file.exists(LOG)) file.remove(LOG)
say <- function(...) { line <- sprintf(...); cat(line, "\n", sep = ""); cat(line, "\n", file = LOG, append = TRUE, sep = "") }
CA <- battery_cells()
key <- function(z) paste(z$block, z$cell)

## ---- (a) the cell list -------------------------------------------------------------------------------------------------------
links3 <- c("stk_long", "stk_short", "stk_asym", "cauchit", "t4", "loglog")
alt <- data.frame(block = c(rep("3", 12), rep("2", 7), "4"),
                  cell = c(sprintf("%s_n%d", rep(links3, each = 2), rep(c(500L, 1000L), 6)),
                           "loglog_base_n610", "loglog_base_n1000", "cauchit_auc_n460", "cauchit_auc_n760", "loglog_auc_n380",
                           "loglog_auc_n640", "loglog_e12_n1000", "crossover_n1000"), stringsAsFactors = FALSE)
alt$n <- as.integer(sub("^.*_n", "", alt$cell))
alt$null <- ifelse(alt$block == "3", sprintf("1a null_link_n%d", alt$n),
            ifelse(alt$block == "2", sprintf("1b null_%s_n%d", sub("^[a-z0-9]+_([a-z0-9]+)_n.*$", "\\1", alt$cell), alt$n),
                   "1b null_crossover_n1000"))
mine <- unique(c(key(alt), alt$null))
say("(a) my E9 list: %d cells (%d alternatives, %d nulls)", length(mine), nrow(alt), length(unique(alt$null)))
say("    driver Cells8 = my list: %s", setequal(key(Cells8), mine))
b2small <- CA$cell[CA$block == "2" & CA$n <= 1000]
say("    block 2 cells with n <= 1000 in battery_cells(): %s", paste(sort(b2small), collapse = ", "))
say("    ... equal E9's seven: %s", setequal(b2small, alt$cell[alt$block == "2"]))
m <- match(key(alt), key(CA))
say("    every alternative's null_cell in battery_cells() is the one E9 names: %s",
    all(paste(CA$null_block[m], CA$null_cell[m]) == alt$null))
mn <- match(alt$null, key(CA))
say("    null has the same n and role 'null': %s", all(CA$n[mn] == alt$n) && all(CA$role[mn] == "null"))
## E1 design parameters, read from the pre-declaration, not from the table
e12_c0 <- c(logit = -2.64, probit = -2.01, t4 = -2.18, cauchit = -3.15, loglog = -1.60, cloglog = -2.76)
exp_design <- function(z) {
  if (z$generator != "design") return(TRUE)
  ds <- z$design; lk <- z$link
  s0 <- switch(ds, base = 1, auc = 2, e12 = 1); c00 <- if (ds == "e12") e12_c0[[if (lk %in% names(e12_c0)) lk else "logit"]] else 0
  isTRUE(all.equal(c(z$s, z$c0), c(s0, c00))) && z$xdist == "uniform"
}
C8 <- CA[match(mine, key(CA)), ]
ok_des <- vapply(seq_len(nrow(C8)), function(i) exp_design(as.list(C8[i, ])), logical(1))
say("    design cells: s, c0 and x ~ U(-3, 3) as E1 (e12 c0 per link, logit null -2.64): %s", all(ok_des))
tab <- C8[, c("block", "cell", "role", "n", "B", "generator", "design", "link", "s", "c0", "family", "scen", "formula",
              "null_block", "null_cell", "seed_family", "seed_base")]
fwrite(tab, file.path(OUTD, "review_cells.csv"))
say("    B range %d-%d (>= 500 needed); seed families: %s", min(C8$B), max(C8$B), paste(unique(C8$seed_family), collapse = ", "))
say("    rows of Cells8 identical to battery_cells(): %s",
    isTRUE(all.equal(Cells8, CA[match(key(Cells8), key(CA)), ], check.attributes = FALSE)))

## ---- (b) regeneration with my own generators ------------------------------------------------------------------------------
my_h <- function(eta, a1, a2) {                                    # Stukel (1988) h-family, E0.1
  h <- eta
  for (i in seq_along(eta)) {
    e <- eta[i]
    if (e >= 0) { if (a1 > 0) h[i] <- (exp(a1 * e) - 1) / a1 else if (a1 < 0) h[i] <- -log(1 - a1 * e) / a1 }
    else { ae <- -e; if (a2 > 0) h[i] <- -(exp(a2 * ae) - 1) / a2 else if (a2 < 0) h[i] <- log(1 - a2 * ae) / a2 }
  }
  h
}
my_p <- function(link, eta) switch(link, logit = plogis(eta), probit = pnorm(eta), cauchit = pcauchy(eta), t4 = pt(eta, 4),
  loglog = exp(-exp(-eta)), cloglog = 1 - exp(-exp(eta)), stk_long = plogis(my_h(eta, -1, -1)),
  stk_short = plogis(my_h(eta, 1, 1)), stk_asym = plogis(my_h(eta, -1, 1)), stop(link))
my_gen <- function(block, cell, n) {
  if (block == "4" && grepl("^crossover_", cell)) {                 # grid_edge_loses.R crossover
    x <- runif(n, -2.5, 2.5); d <- rbinom(n, 1, 0.5); y <- rbinom(n, 1, plogis(0.9 * x - 0.5 * x * d))
  } else if (grepl("^null_crossover_", cell)) {                       # E8.4: working-model limit eta = 0.62x
    x <- runif(n, -2.5, 2.5); d <- rbinom(n, 1, 0.5); y <- rbinom(n, 1, plogis(0.62 * x))
  } else if (block == "1a") {                                          # grid_null.R link null: eta = 0.6x + 0.5d
    x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); y <- rbinom(n, 1, plogis(0.6 * x + 0.5 * d))
  } else {                                                             # E1 designs
    parts <- strsplit(sub("_n[0-9]+$", "", cell), "_")[[1]]
    if (block == "3") { link <- sub("_n[0-9]+$", "", cell); ds <- "base" }
    else if (parts[1] == "null") { link <- "logit"; ds <- parts[2] }
    else { link <- parts[1]; ds <- parts[2] }
    s <- if (ds == "auc") 2 else 1; c0 <- if (ds == "e12") e12_c0[[link]] else 0
    x <- runif(n, -3, 3); d <- rbinom(n, 1, 0.5); y <- rbinom(n, 1, my_p(link, c0 + s * (0.6 * x + 0.5 * d)))
  }
  data.frame(x = x, d = d, y = y)
}
DRY <- edge_battery("dryrun")
res <- list()
for (i in seq_len(nrow(Cells8))) {
  ce <- as.list(Cells8[i, ])
  f <- file.path(DRY, ce$block, paste0(ce$cell, "_pvalues.csv.gz"))
  if (!file.exists(f)) next
  S <- as.data.frame(fread(f))
  for (r in 1:3) {
    st <- S[S$rep == r, ]
    sb <- st$seed - r                                                  # the seed base as the battery wrote it
    RNGkind("L'Ecuyer-CMRG"); set.seed(sb + r)
    D <- my_gen(ce$block, ce$cell, ce$n)
    fit <- suppressWarnings(glm(y ~ x + d, data = D, family = binomial()))
    st_after <- .Random.seed
    pk <- tryCatch(def.gof(fit, G = 10, basis = "poly3", weights = "unit"), error = function(e) NULL)
    p_pkg <- if (is.null(pk)) NA_real_ else as.numeric(if (!is.null(pk$p_value)) pk$p_value else pk$p.value)
    g <- rv_regen(r, ce)
    Xm <- model.matrix(fit); Xd <- stats::model.matrix(g$fq$fit)
    res[[length(res) + 1]] <- data.frame(block = ce$block, cell = ce$cell, rep = r, seed_base_file = sb, seed_base_table = ce$seed_base,
      y_identical = identical(as.numeric(D$y), as.numeric(g$dat$d$y)), x_maxdiff = max(abs(D$x - g$dat$d$x)),
      d_identical = identical(as.numeric(D$d), as.numeric(g$dat$d$d)), X_maxdiff = max(abs(unname(Xm) - unname(Xd))),
      y_fit_identical = identical(unname(fit$y), unname(g$fq$fit$y)), n_ok = st$n == nrow(D), events_ok = st$events == sum(D$y),
      p_stored = st$EDGE.poly3.u.G10, p_pkg = p_pkg, p_driver = g$p, absdiff_pkg = abs(p_pkg - st$EDGE.poly3.u.G10),
      absdiff_driver = abs(g$p - st$EDGE.poly3.u.G10), gate = rv_identity(g, list(seed = st$seed, n = st$n, events = st$events,
      p = st$EDGE.poly3.u.G10), sb + r), stringsAsFactors = FALSE)
  }
}
R <- do.call(rbind, res)
fwrite(R, file.path(OUTD, "review_regen_dryrun.csv"))
say("(b) dry-run cells of block 8 regenerated with my own generators: %s", paste(unique(R$cell), collapse = ", "))
say("    seed base in the file = table: %s; y, d identical and x equal: %s (max |x diff| %.1e); model matrix max diff %.1e",
    all(R$seed_base_file == R$seed_base_table), all(R$y_identical & R$d_identical & R$y_fit_identical), max(R$x_maxdiff), max(R$X_maxdiff))
say("    n and events as stored: %s; package def.gof vs stored EDGE.poly3.u.G10 max |diff| %.2e; driver vs stored %.2e; gate passes %d/%d",
    all(R$n_ok & R$events_ok), max(R$absdiff_pkg, na.rm = TRUE), max(R$absdiff_driver, na.rm = TRUE), sum(R$gate), nrow(R))

## ---- (c) the real battery files that exist now --------------------------------------------------------------------------------
REAL <- edge_battery()
rr <- list(); gg <- list()
for (i in seq_len(nrow(Cells8))) {
  ce <- as.list(Cells8[i, ])
  f <- edge_battery(ce$block, paste0(ce$cell, "_pvalues.csv.gz"))
  if (!file.exists(f)) { rr[[length(rr) + 1]] <- data.frame(block = ce$block, cell = ce$cell, exists = FALSE, stringsAsFactors = FALSE); next }
  S <- as.data.frame(fread(f, select = c("rep", "seed", "n", "events", "EDGE.poly3.u.G10", "flag.degenerate", "flag.rep_error")))
  S5 <- S[S$rep %in% 1:500, ]
  dryf <- file.path(DRY, ce$block, paste0(ce$cell, "_pvalues.csv.gz"))
  same_dry <- NA
  if (file.exists(dryf)) { Dd <- as.data.frame(fread(dryf, select = c("rep", "seed", "EDGE.poly3.u.G10")))
    mm <- match(Dd$rep, S$rep); same_dry <- all(Dd$seed == S$seed[mm]) && max(abs(Dd$EDGE.poly3.u.G10 - S$EDGE.poly3.u.G10[mm])) == 0 }
  rr[[length(rr) + 1]] <- data.frame(block = ce$block, cell = ce$cell, exists = TRUE, rows = nrow(S), B = ce$B,
    reps_1_500 = all(1:500 %in% S$rep), seed_is_base_plus_rep = all(S$seed - S$rep == ce$seed_base), n_ok = all(S$n == ce$n),
    na_edge_1_500 = sum(!is.finite(S5$EDGE.poly3.u.G10)), degenerate_1_500 = sum(S5$flag.degenerate %in% 1),
    rep_error_1_500 = sum(S5$flag.rep_error %in% 1), min_events_1_500 = min(S5$events), same_as_dryrun = same_dry, stringsAsFactors = FALSE)
  for (r in c(1L, 200L, 500L)) {
    st <- as.list(S[S$rep == r, ]); g <- rv_regen(r, ce)
    gg[[length(gg) + 1]] <- data.frame(block = ce$block, cell = ce$cell, rep = r, absdiff = abs(g$p - st$EDGE.poly3.u.G10),
      gate = rv_identity(g, list(seed = st$seed, n = st$n, events = st$events, p = st$EDGE.poly3.u.G10), ce$seed_base + r),
      stringsAsFactors = FALSE)
  }
}
RR <- rbindlist(rr, fill = TRUE); GG <- rbindlist(gg)
fwrite(RR, file.path(OUTD, "review_real_files.csv")); fwrite(GG, file.path(OUTD, "review_real_gate.csv"))
say("(c) real battery files of block 8 cells present now: %d of 30 (%s missing)", sum(RR$exists), paste(RR$cell[!RR$exists], collapse = ", "))
E <- RR[RR$exists == TRUE]
say("    replicates 1-500 present: %s; seed = seed_base + rep on every row: %s; n as the table: %s",
    all(E$reps_1_500), all(E$seed_is_base_plus_rep), all(E$n_ok))
say("    in replicates 1-500: EDGE.poly3.u.G10 missing %d, degenerate %d, replicate errors %d, fewest events %d",
    sum(E$na_edge_1_500), sum(E$degenerate_1_500), sum(E$rep_error_1_500), min(E$min_events_1_500))
say("    real replicates equal the dry-run replicates where both exist: %s", paste(na.omit(E$same_as_dryrun), collapse = ", "))
say("    identity gate on replicates 1, 200, 500 of the present cells: %d/%d pass, max |diff| %.2e", sum(GG$gate), nrow(GG), max(GG$absdiff))
say("done")
