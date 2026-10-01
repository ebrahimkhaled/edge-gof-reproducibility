## analysis_review_indep.R -- a second, separately written implementation of Section A and H1-H5 of
## paper_EDGE/theory/PREDECLARATION_restructure_battery.md (A, D, E12), run beside analyse_M_battery.R on three fresh
## synthetic summary sets built from battery/cells.csv (seeds 5550101, 5550202, 5550303), with per-replicate files for the
## paired tests of E12.6. The cell lists are typed from the text of E3, E7, E12.4 and E12.5, not matched by pattern.
## No value of the real battery is read. Everything is written under battery/_review/analysis/review/.
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
RV <- edge_battery("_review", "analysis", "review")
suppressPackageStartupMessages(library(data.table))
dir.create(RV, recursive = TRUE, showWarnings = FALSE)
sink(file.path(RV, "review_indep.log"), split = TRUE)
cat("analysis_review_indep.R ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")

CT <- fread(edge_battery("cells.csv"),
            colClasses = list(character = c("block", "cell", "role", "null_block", "null_cell", "param", "design", "link",
                                            "family", "generator", "G_list")))
CT[is.na(null_block), null_block := ""]; CT[is.na(null_cell), null_cell := ""]

RES <- data.frame()
SET <- "cells"
check <- function(label, ok, detail = "") {
  RES <<- rbind(RES, data.frame(set = SET, check = label, ok = isTRUE(ok), detail = detail, stringsAsFactors = FALSE))
  cat(sprintf("  [%s] %s  %s\n", if (isTRUE(ok)) "ok" else "DIFF", label, detail))
}
same_num <- function(a, b, tol = 1e-10) length(a) == length(b) && all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol))

## ---- the cell lists, typed from the pre-declaration ------------------------------------------------------------------------------
N3 <- c(200, 500, 1000, 2000, 5000)
nm <- function(stem, ns) sprintf("%s_n%d", stem, as.integer(ns))
grid <- function(stem, sev, ns) as.vector(outer(paste0(stem, "_", sev), ns, nm))
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
  data.table(family = 4L, block = "3", cell = c(grid("quad", c("0.01", "0.02", "0.03", "0.05", "0.1", "0.2", "0.4"), N3),
                                                grid("binint", c("0.1", "0.2", "0.3", "0.5", "0.7"), N3),
                                                grid("contint", c("0.1", "0.3", "0.5", "0.7"), N3))),
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

## ---- the builder's script, functions only ----------------------------------------------------------------------------------------
BE <- new.env()
sys.source(file.path(SIMDIR, "analyse_M_battery.R"), envir = BE)
run_builder <- function(args) { res <- NULL; txt <- capture.output(res <- BE$an_main(args)); res }

cat("\n== cell table: the typed lists against the real cell table and against the builder's patterns\n")
key <- function(d) paste(d$block, d$cell)
check("every typed family cell exists in cells.csv as an alternative",
      all(key(FAMS) %in% key(CT[role == "alternative"])), paste(setdiff(key(FAMS), key(CT)), collapse = ", "))
check("family sizes", identical(as.integer(table(FAMS$family)), c(32L, 21L, 10L, 80L, 13L, 2L)), paste(table(FAMS$family), collapse = " "))
check("no cell name appears in two blocks", !any(duplicated(CT$cell)), paste(CT$cell[duplicated(CT$cell)], collapse = ", "))
Mb <- BE$an_membership(BE$an_cell_table(edge_battery()), character(0))
check("builder membership on the real cell table = typed lists (same cells, same families)",
      setequal(paste(Mb$block, Mb$cell, Mb$family), paste(FAMS$block, FAMS$cell, FAMS$family)) && nrow(Mb) == nrow(FAMS),
      sprintf("builder %d rows, typed %d", nrow(Mb), nrow(FAMS)))
out_alt <- CT[block %in% c("2", "3", "4") & role == "alternative" & !(key(CT[block %in% c("2", "3", "4") & role == "alternative"]) %in% key(FAMS))]
cat("  alternatives of blocks 2-4 in no family:", nrow(out_alt), "\n   ", paste(out_alt$block, out_alt$cell, collapse = "; "), "\n")
check("H4 cells: 21 block 3 cells of families 2-4 at n = 1000 plus 8 census cells, all in cells.csv",
      nrow(H4CELLS) == 29 && all(key(H4CELLS) %in% key(CT)))
check("H5 cells in cells.csv, skew n = base n", all(key(rbind(H5B, H5S)) %in% key(CT)) &&
        setequal(CT$n[match(key(H5B), key(CT))], CT$n[match(key(H5S), key(CT))]))

## ---- synthetic roots ------------------------------------------------------------------------------------------------------------
RV_BASES4 <- c("poly3", "poly2", "stk", "sym")
RV_EDGE <- as.vector(outer(as.vector(outer(paste0("EDGE.", RV_BASES4), c("u", "sc"), paste, sep = ".")), c("G10", "Grule"), paste, sep = "."))
RV_OTHER <- c("HL.G10", "HL.Grule", "HLF.G10", "HLF.Grule", "Stk.joint", "Stk.LR", "Stk.sym1", "Stk.marg", "GiViTI", "GiViTI.t50",
              "Cubic.LR", "HL_w", "PH", "Tsiatis", "Xie", "PR")
RV_RIV <- c("HL_w", "PH", "Tsiatis", "Xie", "PR")
PCOLS <- c(paste0("EDGE.", rep(c("poly3", "sym"), each = 4), ".", c("u", "sc"), ".", rep(c("G10", "Grule"), each = 2)),
           "GiViTI", "Stk.joint", "Stk.sym1", "HLF.G10", "HL.G10", "HLF.Grule", "HL.Grule", "HL_w", "PH", "Tsiatis", "Xie")

rv_make <- function(dir, seed, prm) {
  set.seed(seed)
  unlink(dir, recursive = TRUE); dir.create(dir, recursive = TRUE)
  file.copy(edge_battery("cells.csv"), file.path(dir, "cells.csv"))
  C <- merge(CT[block %in% c("1a", "1b", "2", "3", "4", "5", "6", "7") & role %in% c("null", "alternative")],
             FAMS[, .(block, cell, fam = family, fdesign = design)], by = c("block", "cell"), all.x = TRUE, sort = FALSE)
  tests0 <- c(RV_EDGE, RV_OTHER)
  L <- vector("list", nrow(C))
  for (i in seq_len(nrow(C))) {
    ce <- C[i]; Bc <- ce$B
    tests <- if (ce$generator == "external") setdiff(tests0, grep("\\.G(10|rule)$", tests0, value = TRUE)) else tests0
    notrun <- ce$n >= 10000 & tests %in% RV_RIV
    if (ce$role == "null") {
      k <- prm$k(tests, ce)
      R <- rbindlist(lapply(c(0.01, 0.05, 0.10), function(a) {
        r <- rbinom(length(tests), Bc, pmin(1, a * k)) / Bc
        r[notrun] <- 0
        data.table(test = tests, alpha = a, rejection = r, size_adj_power = NA_real_, null_size = NA_real_,
                   status = ifelse(notrun, "not run (n >= 10,000)", ""))
      }))
    } else {
      pw <- runif(1, 0.05, 0.9) + rnorm(length(tests), 0, prm$noise) + prm$off(tests, ce)
      e <- startsWith(tests, "EDGE.") & grepl("\\.sc\\.", tests)
      if (!is.na(ce$fam)) pw[e] <- pw[e] + prm$shift[cbind(match(sub("^EDGE\\.([a-z0-9]+)\\..*$", "\\1", tests[e]), RV_BASES4), ce$fam)]
      pw <- pmin(1, pmax(0, pw))
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
    if (i %% 4L == 0L) {          # decoy subsets with far-off numbers: using them instead of subset = all would show
      D1 <- copy(R)[, `:=`(subset = "events<Grule", B = as.integer(round(Bc / 3)), rejection = 0.5, mcse = 0.001, size_adj_power = 0.001)]
      D2 <- copy(R)[, `:=`(subset = "events>=Grule", B = Bc - as.integer(round(Bc / 3)), rejection = 0.5, mcse = 0.001, size_adj_power = 0.999)]
      R <- rbind(R, D1, D2)
    }
    FL <- data.table(test = c("flag.evlt.G10", "flag.evlt.Grule", "flag.degenerate"), alpha = NA_real_,
                     rejection = c(0, if (i %% 4L == 0L) 0.3 else 0, 0), size_adj_power = NA_real_, null_size = NA_real_,
                     status = "flag rate: share of replicates with the flag = 1", block = ce$block, cell = ce$cell, role = ce$role,
                     null_cell = ce$null_cell, n = ce$n, subset = "all", B = Bc, mcse = NA_real_, declined = 0, rejection_given_p = NA_real_)
    L[[i]] <- rbind(FL, R, use.names = TRUE)
  }
  S <- rbindlist(L, use.names = TRUE)
  setcolorder(S, c("block", "cell", "role", "null_cell", "n", "subset", "B", "test", "alpha", "rejection", "mcse", "size_adj_power",
                   "null_size", "declined", "rejection_given_p", "status"))
  for (b in unique(S$block)) { dir.create(file.path(dir, b), showWarnings = FALSE); fwrite(S[block == b], file.path(dir, b, "_summary.csv")) }
  ## per-replicate files for the cells of H1-H5
  PC <- unique(rbind(FAMS[family %in% 1:2, .(block, cell)], H4CELLS, H5B, H5S))
  W5 <- S[subset == "all" & !is.na(alpha) & abs(alpha - 0.05) < 1e-9]
  wk <- paste(W5$block, W5$cell, W5$test)
  Br <- 300L
  for (j in seq_len(nrow(PC))) {
    if (j %in% prm$nofile) next
    bj <- PC$block[j]; cj <- PC$cell[j]
    z <- rnorm(Br)
    P <- data.table(rep = seq_len(Br), seed = 1e6 + seq_len(Br), n = CT$n[CT$block == bj & CT$cell == cj], events = 100L)
    for (tt in PCOLS) {
      pw <- W5$size_adj_power[match(paste(bj, cj, tt), wk)]
      if (!length(pw) || !is.finite(pw)) pw <- 0.3
      pw <- min(0.995, max(0.005, pw))
      p <- pnorm(qnorm(0.05) - qnorm(pw) + prm$rho * z + sqrt(1 - prm$rho^2) * rnorm(Br))
      p[sample(Br, 2)] <- 0.05; p[runif(Br) < 0.03] <- NA
      set(P, j = tt, value = p)
    }
    if (j %in% prm$dropcol) P[, Stk.joint := NULL]
    fwrite(P, file.path(dir, bj, paste0(cj, "_pvalues.csv.gz")), compress = "gzip")
  }
  invisible(S)
}

## ---- the reviewer's implementation --------------------------------------------------------------------------------------------------
rv_read <- function(root) {
  S <- rbindlist(lapply(c("0", "1a", "1b", "2", "3", "4", "5", "6", "7"), function(b) {
    f <- file.path(root, b, "_summary.csv")
    if (file.exists(f)) fread(f, colClasses = list(character = c("block", "cell", "role", "null_cell", "subset", "test", "status"))) else NULL
  }), use.names = TRUE, fill = TRUE)
  A <- S[subset == "all" & !is.na(alpha) & !startsWith(test, "flag.")]
  A[!is.na(status) & startsWith(status, "not run"), `:=`(rejection = NA_real_, mcse = NA_real_, size_adj_power = NA_real_)]
  A[, k4 := paste(block, cell, test, sprintf("%.3f", alpha))]
  A
}
gv <- function(A, bl, ce, tt, col, a = 0.05) A[[col]][match(paste(bl, ce, tt, sprintf("%.3f", a)), A$k4)]

rv_sectionA <- function(A, arm) {
  gate <- rbindlist(lapply(c("poly3", "sym", "stk"), function(b) rbindlist(lapply(c("u", "sc"), function(f) {
    tt <- sprintf("EDGE.%s.%s.%s", b, f, arm)
    nc <- unique(A[role == "null" & test == tt, .(block, cell)])
    r5 <- gv(A, nc$block, nc$cell, tt, "rejection", 0.05); m5 <- gv(A, nc$block, nc$cell, tt, "mcse", 0.05)
    r1 <- gv(A, nc$block, nc$cell, tt, "rejection", 0.01); m1 <- gv(A, nc$block, nc$cell, tt, "mcse", 0.01)
    data.table(basis = b, form = f, test = tt, block = nc$block, cell = nc$cell, r5 = r5, m5 = m5, r1 = r1, m1 = m1,
               ok5 = r5 <= 0.05 + 3 * m5, ok1 = r1 <= 0.01 + 3 * m1)
  }))))
  gate[, ok := ok5 & ok1]
  fam <- rbindlist(lapply(c("poly3", "sym", "stk"), function(b) rbindlist(lapply(c("u", "sc"), function(f) {
    p <- gv(A, FAMS$block, FAMS$cell, sprintf("EDGE.%s.%s.%s", b, f, arm), "size_adj_power")
    data.table(basis = b, form = f, family = 1:6, mean = vapply(1:6, function(k) mean(p[FAMS$family == k]), 0))
  }))))
  dec <- rbindlist(lapply(c("poly3", "sym", "stk"), function(b) {
    eu <- all(gate[basis == b & form == "u"]$ok); es <- all(gate[basis == b & form == "sc"]$ok)
    fu <- fam[basis == b & form == "u"][order(family)]$mean; fs <- fam[basis == b & form == "sc"][order(family)]$mean
    mu <- mean(fu); ms <- mean(fs); d <- ms - mu
    branch <- ""
    if (eu && es) {
      if (d >= 0.014) { ch <- if (max(fu - fs) > 0.10) "unit" else "score"; branch <- if (ch == "unit") "score higher, vetoed" else "score higher" }
      else if (-d >= 0.014) { ch <- if (max(fs - fu) > 0.10) "score" else "unit"; branch <- if (ch == "score") "unit higher, vetoed" else "unit higher" }
      else { ch <- "unit"; branch <- "tie" }
    } else if (eu) { ch <- "unit"; branch <- "only unit eligible" } else if (es) { ch <- "score"; branch <- "only score eligible" } else { ch <- "unit"; branch <- "neither eligible" }
    data.table(basis = b, elig_u = eu, elig_s = es, macro_u = mu, macro_s = ms, diff = d, loss_s = max(fu - fs), loss_u = max(fs - fu),
               choice = ch, branch = branch)
  }))
  list(gate = gate, fam = fam, dec = dec)
}

rv_census <- function(A, ver, f) {
  edge <- sprintf("EDGE.poly3.%s.%s", f, ver)
  part <- c(paste0("HLF.", ver), paste0("HL.", ver), "HL_w", "PH"); six <- c(part, "Tsiatis", "Xie")
  K <- rbindlist(lapply(seq_len(nrow(H4CELLS)), function(i) {
    g <- function(tt) gv(A, H4CELLS$block[i], H4CELLS$cell[i], tt, "size_adj_power")
    p3 <- g(edge); rp <- sapply(part, g); r6 <- sapply(six, g)
    data.table(block = H4CELLS$block[i], cell = H4CELLS$cell[i], p3 = p3, best = max(rp, na.rm = TRUE), best_name = names(rp)[which.max(rp)],
               detectable = max(c(p3, r6), na.rm = TRUE) >= 0.15, saturated = min(c(p3, rp), na.rm = TRUE) >= 0.97,
               lead = p3 > max(rp, na.rm = TRUE) - 0.014, lead6 = p3 > max(r6, na.rm = TRUE) - 0.014)
  }))
  K[, headline := detectable & !saturated]
  K
}

rv_hyp <- function(A, arm) {
  H <- list(); PR <- list(); CEN <- list()
  F1 <- FAMS[family == 1]; F2 <- FAMS[family == 2]
  sap <- function(d, tt) gv(A, d$block, d$cell, tt, "size_adj_power")
  addp <- function(hyp, f, d, ta, tb) PR[[length(PR) + 1]] <<- data.table(hypothesis = hyp, form = f, block = d$block, cell = d$cell, test_a = ta, test_b = tb)
  for (f in c("u", "sc")) {
    e3 <- sprintf("EDGE.poly3.%s.%s", f, arm); es <- sprintf("EDGE.sym.%s.%s", f, arm)
    d1 <- sap(F1, e3) - sap(F1, "GiViTI")
    dm <- vapply(c("base", "auc", "e12"), function(g) mean(d1[F1$design == g]), 0)
    H[[length(H) + 1]] <- data.table(hypothesis = "H1", form = f, statistic = mean(d1), verdict = if (mean(d1) >= 0.10 && all(dm >= -0.014)) "holds" else "against",
                                     mean_base = dm[["base"]], mean_auc = dm[["auc"]], mean_e12 = dm[["e12"]], below = paste(names(dm)[dm < -0.014], collapse = " "))
    addp("H1", f, F1, e3, "GiViTI")
    d2 <- sap(F1, es) - sap(F1, "Stk.joint")
    H[[length(H) + 1]] <- data.table(hypothesis = "H2", form = f, statistic = mean(d2), verdict = if (mean(d2) >= -0.014) "holds" else "against")
    addp("H2", f, F1, es, "Stk.joint")
    d3 <- sap(F2, e3) - sap(F2, "GiViTI")
    H[[length(H) + 1]] <- data.table(hypothesis = "H3", form = f, statistic = mean(d3), verdict = if (mean(d3) >= -0.05) "holds" else "against")
    addp("H3", f, F2, e3, "GiViTI")
    for (ver in c("G10", "Grule")) {
      K <- rv_census(A, ver, f); hyp <- if (ver == "G10") "H4" else "H4.ruleG"
      led <- sum(K$headline & K$lead)
      H[[length(H) + 1]] <- data.table(hypothesis = hyp, form = f, statistic = led, headline = sum(K$headline), led6 = sum(K$headline & K$lead6),
                                       verdict = if (ver == "Grule") "reported" else if (led >= 19) "holds" else "against")
      CEN[[length(CEN) + 1]] <- cbind(data.table(version = ver, form = f), K)
      addp(hyp, f, H4CELLS, sprintf("EDGE.poly3.%s.%s", f, ver), K$best_name)
    }
    d5 <- sap(H5B, "Stk.sym1") - sap(H5B, es)
    H[[length(H) + 1]] <- data.table(hypothesis = "H5", form = f, statistic = mean(d5), verdict = if (mean(d5) <= 0.02) "holds" else "against")
    addp("H5", f, H5B, "Stk.sym1", es)
    d6 <- sap(H5S, "Stk.sym1") - sap(H5S, es)
    H[[length(H) + 1]] <- data.table(hypothesis = "H5.skew", form = f, statistic = mean(d6), verdict = "reported")
    addp("H5.skew", f, H5S, "Stk.sym1", es)
  }
  list(hyp = rbindlist(H, fill = TRUE), pairs = rbindlist(PR), census = rbindlist(CEN))
}

holm_by_hand <- function(p) {
  m <- length(p); o <- order(p); adj <- numeric(m); run <- 0
  for (j in seq_len(m)) { run <- max(run, min(1, (m - j + 1) * p[o[j]])); adj[o[j]] <- run }
  adj
}
rv_mcnemar <- function(root, A, PR) {
  X <- merge(PR, CT[, .(block, cell, null_block, null_cell)], by = c("block", "cell"), all.x = TRUE, sort = FALSE)
  X[, `:=`(sa = gv(A, null_block, null_cell, test_a, "rejection"), ma = gv(A, null_block, null_cell, test_a, "mcse"),
           sb = gv(A, null_block, null_cell, test_b, "rejection"), mb = gv(A, null_block, null_cell, test_b, "mcse"))]
  X[, hold := !is.na(sa) & !is.na(sb) & sa <= 0.05 + 3 * ma & sb <= 0.05 + 3 * mb]
  X[, `:=`(n11 = NA_integer_, n10 = NA_integer_, n01 = NA_integer_, n00 = NA_integer_, p = NA_real_)]
  cache <- list()
  for (i in seq_len(nrow(X))) {
    f <- file.path(root, X$block[i], paste0(X$cell[i], "_pvalues.csv.gz"))
    if (!file.exists(f)) next
    if (is.null(cache[[f]])) cache[[f]] <- fread(f)
    P <- cache[[f]]
    if (!(X$test_a[i] %in% names(P)) || !(X$test_b[i] %in% names(P))) next
    ra <- !is.na(P[[X$test_a[i]]]) & P[[X$test_a[i]]] <= 0.05; rb <- !is.na(P[[X$test_b[i]]]) & P[[X$test_b[i]]] <= 0.05
    n10 <- sum(ra & !rb); n01 <- sum(!ra & rb)
    set(X, i, c("n11", "n10", "n01", "n00"), list(sum(ra & rb), n10, n01, sum(!ra & !rb)))
    if (X$hold[i]) set(X, i, "p", if (n10 + n01 == 0) 1 else binom.test(n10, n10 + n01, 0.5)$p.value)
  }
  X[, used := !is.na(p)]
  X[, p_holm := NA_real_]
  X[used == TRUE, p_holm := holm_by_hand(p), by = .(hypothesis, form)]
  X
}

## ---- compare the two implementations on one root ------------------------------------------------------------------------------------
FORMNAME <- c(u = "unit", sc = "score")
compare_root <- function(root, label) {
  SET <<- label
  cat("\n== ", label, "  ", root, "\n", sep = "")
  A <- rv_read(root)
  mine <- rv_sectionA(A, "Grule")
  forms <- setNames(ifelse(mine$dec$choice == "unit", "u", "sc"), mine$dec$basis)
  cat("  reviewer's decision:", paste(sprintf("%s %s (%s; diff %+.4f, loss score %.3f, loss unit %.3f)", mine$dec$basis, mine$dec$choice,
                                              mine$dec$branch, mine$dec$diff, mine$dec$loss_s, mine$dec$loss_u), collapse = "; "), "\n")
  cat("  gate failures (reviewer):", paste(mine$gate[ok == FALSE, sprintf("%s %s/%s", test, block, cell)], collapse = "; "), "\n")
  b <- run_builder(c("--root", root, "--out", file.path(root, "out_builder")))
  ## gate
  G <- merge(b$gate[, .(test, block, cell, pass_05, pass_01, pass, z_05)], mine$gate, by = c("test", "block", "cell"), all = TRUE)
  check("gate: same null cells per variant", !anyNA(G$pass) && !anyNA(G$ok), sprintf("%d rows", nrow(G)))
  check("gate: same pass at 0.05, at 0.01 and overall", identical(G$pass_05, G$ok5) && identical(G$pass_01, G$ok1) && identical(G$pass, G$ok),
        sprintf("%d failing rows", sum(!G$ok, na.rm = TRUE)))
  check("gate: z at 0.05 = (size - 0.05) / MCSE", same_num(G[m5 > 0]$z_05, (G[m5 > 0]$r5 - 0.05) / G[m5 > 0]$m5, 1e-9))
  ## families, decision
  Fb <- b$families[family != "macro", .(basis, form, family = as.integer(family), mean_power)]
  FF <- merge(Fb, mine$fam, by = c("basis", "form", "family"), all = TRUE)
  check("family means (3 bases x 2 forms x 6 families)", nrow(FF) == 36 && same_num(FF$mean_power, FF$mean))
  D <- merge(b$decision, mine$dec, by = "basis")
  check("eligibility", identical(D$eligible_unit, D$elig_u) && identical(D$eligible_score, D$elig_s))
  check("macro-averages, difference, worst losses", same_num(D$macro_unit, D$macro_u) && same_num(D$macro_score, D$macro_s) &&
          same_num(D$difference, D$diff) && same_num(D$worst_loss_score, D$loss_s) && same_num(D$worst_loss_unit, D$loss_u))
  check("choice per basis", identical(D$choice, D$choice.y) || identical(D$choice.x, D$choice.y),
        paste(D$basis, D[[if ("choice.x" %in% names(D)) "choice.x" else "choice"]], collapse = ", "))
  check("forms used for H1-H5", identical(unname(b$forms_used[c("poly3", "sym")]), unname(forms[c("poly3", "sym")])))
  ## hypotheses
  mh <- rv_hyp(A, "Grule")
  Hb <- copy(b$hypotheses); Hb[, form2 := ifelse(form == "unit", "u", "sc")][, form := NULL]
  HH <- merge(Hb, mh$hyp, by.x = c("hypothesis", "form2"), by.y = c("hypothesis", "form"), all = TRUE, suffixes = c(".b", ".r"))
  check("hypotheses: same rows (H1, H2, H3, H4, H4.ruleG, H5, H5.skew x 2 forms)", nrow(HH) == 14 && !anyNA(HH$verdict.b) && !anyNA(HH$verdict.r))
  check("hypotheses: statistics", same_num(HH$statistic.b, HH$statistic.r, 1e-10),
        paste(sprintf("%s/%s %.4f", HH$hypothesis, HH$form2, HH$statistic.r), collapse = "; "))
  check("hypotheses: verdicts", identical(HH$verdict.b, HH$verdict.r), paste(sprintf("%s/%s %s", HH$hypothesis, HH$form2, HH$verdict.r), collapse = "; "))
  h1 <- HH[hypothesis == "H1"]
  check("H1: design means and designs below -0.014", same_num(h1$mean_base.b, h1$mean_base.r) && same_num(h1$mean_auc.b, h1$mean_auc.r) &&
          same_num(h1$mean_e12.b, h1$mean_e12.r) && identical(h1$designs_below_margin, h1$below),
        paste(sprintf("%s: base %.4f auc %.4f e12 %.4f below '%s'", h1$form2, h1$mean_base.r, h1$mean_auc.r, h1$mean_e12.r, h1$below), collapse = "; "))
  h4 <- HH[hypothesis %in% c("H4", "H4.ruleG")]
  check("H4: headline cells, led or tied", identical(as.integer(h4$headline_cells), as.integer(h4$headline)) && identical(as.integer(h4$led_or_tied), as.integer(h4$statistic.r)),
        paste(sprintf("%s/%s %d of %d", h4$hypothesis, h4$form2, as.integer(h4$statistic.r), as.integer(h4$headline)), collapse = "; "))
  check("forms_used flags on the hypothesis rows", all(HH[hypothesis %in% c("H1", "H3", "H4", "H4.ruleG")]$forms_used == (HH[hypothesis %in% c("H1", "H3", "H4", "H4.ruleG")]$form2 == forms[["poly3"]])) &&
          all(HH[hypothesis %in% c("H2", "H5", "H5.skew")]$forms_used == (HH[hypothesis %in% c("H2", "H5", "H5.skew")]$form2 == forms[["sym"]])))
  ## census
  Kb <- copy(b$census); Kb[, form2 := ifelse(form == "unit", "u", "sc")][, form := NULL]
  KK <- merge(Kb, mh$census, by.x = c("version", "form2", "block", "cell"), by.y = c("version", "form", "block", "cell"), all = TRUE)
  check("census cell by cell: detectable, saturated, leads/ties partition, all six, headline, best partition test",
        nrow(KK) == 116 && identical(KK$detectable.x, KK$detectable.y) && identical(KK$saturated.x, KK$saturated.y) &&
          identical(KK$beats_partition, KK$lead) && identical(KK$beats_all6, KK$lead6) && identical(KK$in_headline, KK$headline) &&
          identical(KK$best_partition_name, KK$best_name), sprintf("%d rows", nrow(KK)))
  ## McNemar and Holm
  mm <- rv_mcnemar(root, A, mh$pairs)
  Mc <- copy(b$mcnemar); Mc[, form2 := ifelse(form == "unit", "u", "sc")][, form := NULL]
  MM <- merge(Mc, mm, by.x = c("hypothesis", "form2", "block", "cell"), by.y = c("hypothesis", "form", "block", "cell"), all = TRUE, suffixes = c(".b", ".r"))
  check("McNemar: same rows and test pairs", nrow(MM) == nrow(mm) && nrow(MM) == nrow(Mc) && identical(MM$test_a.b, MM$test_a.r) &&
          identical(MM$test_b.b, MM$test_b.r), sprintf("%d rows", nrow(MM)))
  check("McNemar: size held (both tests, matched null) and used", identical(MM$holds_a & MM$holds_b, MM$hold) && identical(MM$used.b, MM$used.r),
        sprintf("used %d of %d; builder used %d", sum(MM$used.r), nrow(MM), sum(MM$used.b)))
  cnt <- MM[!is.na(n11)]
  check("McNemar: 2x2 counts on the shared replicates", nrow(cnt) > 0 && identical(as.integer(cnt$both), cnt$n11) && identical(as.integer(cnt$a_only), cnt$n10) &&
          identical(as.integer(cnt$b_only), cnt$n01) && identical(as.integer(cnt$neither), cnt$n00), sprintf("%d rows with a file and both columns", nrow(cnt)))
  check("McNemar: exact p", same_num(MM$p_mcnemar, MM$p, 1e-12))
  check("Holm within each hypothesis and form", same_num(MM$p_holm.b, MM$p_holm.r, 1e-12),
        paste(MM[used.r == TRUE, .(sig = sum(p_holm.r <= 0.05), m = .N), by = .(hypothesis, form2)][, sprintf("%s/%s %d of %d", hypothesis, form2, sig, m)], collapse = "; "))
  cat("  McNemar rows not used, by reason (builder notes):\n")
  print(Mc[used == FALSE, .N, by = .(reason = sub(" \\(.*$", "", sub("^(.*) fails size.*$", "fails size", note)))][order(-N)][1:min(.N, 6)], row.names = FALSE)
  ## the G = 10 arm
  b10 <- run_builder(c("--root", root, "--out", file.path(root, "out_builder"), "--arm", "10"))
  m10 <- rv_sectionA(A, "G10")
  G10 <- merge(b10$gate[, .(test, block, cell, pass)], m10$gate[, .(test, block, cell, ok)], by = c("test", "block", "cell"), all = TRUE)
  D10 <- merge(b10$decision, m10$dec, by = "basis")
  check("--arm 10: gate, macro-averages and choice at G = 10", !anyNA(G10$pass) && !anyNA(G10$ok) && identical(G10$pass, G10$ok) &&
          same_num(D10$macro_unit, D10$macro_u) && same_num(D10$macro_score, D10$macro_s) &&
          identical(D10[[if ("choice.x" %in% names(D10)) "choice.x" else "choice"]], D10$choice.y),
        paste(D10$basis, D10$choice.y, D10$branch, collapse = ", "))
  mh10 <- rv_hyp(A, "G10")
  H10 <- copy(b10$hypotheses); H10[, form2 := ifelse(form == "unit", "u", "sc")][, form := NULL]
  HH10 <- merge(H10, mh10$hyp, by.x = c("hypothesis", "form2"), by.y = c("hypothesis", "form"), suffixes = c(".b", ".r"))
  check("--arm 10: H1-H3, H5 at G10, H4 both versions, forms from the rule-G choice",
        nrow(HH10) == 14 && same_num(HH10$statistic.b, HH10$statistic.r) && identical(HH10$verdict.b, HH10$verdict.r) &&
          identical(unname(b10$forms_used[c("poly3", "sym")]), unname(forms[c("poly3", "sym")])))
  check("--arm 10: the eight _G10 files", all(file.exists(file.path(root, "out_builder", paste0(c("rule_A_gate", "rule_A_membership", "rule_A_families", "hypotheses", "h4_census", "declined", "mcnemar"), "_G10.csv")))) &&
          file.exists(file.path(root, "out_builder", "rule_A_decision_G10.txt")))
  fwrite(MM, file.path(RV, paste0("compare_mcnemar_", label, ".csv")))
  fwrite(HH, file.path(RV, paste0("compare_hypotheses_", label, ".csv")))
  invisible(list(mine = mine, builder = b))
}

## ---- the three sets -----------------------------------------------------------------------------------------------------------------
shift_mat <- function(v) matrix(v, nrow = 4, ncol = 6, byrow = TRUE, dimnames = list(RV_BASES4, NULL))
set.seed(4040)
P1 <- list(noise = 0.03, rho = 0.5, nofile = 5L, dropcol = 10L,
           k = function(tests, ce) { k <- runif(length(tests), 0.55, 1.05)
                                     if (ce$cell %in% c("null_link_n2000", "null_quad_n5000")) k[tests == "EDGE.sym.u.Grule"] <- 1.9
                                     k },
           off = function(tests, ce) ifelse(tests == "GiViTI", -0.13, 0),
           shift = shift_mat(c(0.06, 0.06, 0.06, 0.06, 0.06, -0.15,
                               rnorm(6, 0, 0.05),
                               0.004, 0.004, 0.004, 0.004, 0.004, 0.004,
                               0.02, 0.03, 0.02, 0.03, 0.02, 0.03)))
P2 <- list(noise = 0.05, rho = 0.3, nofile = c(2L, 60L), dropcol = integer(0),
           k = function(tests, ce) { k <- ifelse(startsWith(tests, "EDGE."), runif(length(tests), 0.85, 1.2), runif(length(tests), 0.7, 1.15))
                                     if (ce$block == "1b") k[tests == "GiViTI"] <- 1.35
                                     k },
           off = function(tests, ce) 0,
           shift = NULL)
P3 <- list(noise = 0.025, rho = 0.7, nofile = integer(0), dropcol = 3L,
           k = function(tests, ce) { k <- runif(length(tests), 0.6, 1.0)
                                     if (ce$cell == "null_base_n2400") k[tests == "EDGE.stk.u.Grule"] <- 1.9
                                     if (ce$cell == "null_link_n500") k[tests == "EDGE.sym.sc.G10"] <- 2.2
                                     k },
           off = function(tests, ce) ifelse(tests == "GiViTI", ifelse(identical(ce$fdesign, "auc"), 0.03, -0.11), ifelse(tests == "Stk.sym1", 0.02, 0)),
           shift = shift_mat(c(-0.06, -0.06, -0.06, -0.06, -0.06, 0.15,
                               rnorm(6, 0, 0.05),
                               -0.03, -0.03, -0.03, -0.03, -0.03, -0.03,
                               0.02, 0.02, 0.02, 0.02, 0.02, 0.02)))
P2$shift <- shift_mat(rnorm(24, 0, 0.05))
roots <- file.path(RV, c("set1", "set2", "set3"))
rv_make(roots[1], 5550101, P1)
rv_make(roots[2], 5550202, P2)
rv_make(roots[3], 5550303, P3)
for (i in 1:3) compare_root(roots[i], paste0("set", i))

## ---- guard: path spellings of the real root (no data read) --------------------------------------------------------------------------
SET <- "guard"
cat("\n== guard\n")
check("real root recognised in other spellings (backslashes, trailing slash, dot segments, case)",
      BE$an_is_real(edge_battery()) &&
        BE$an_is_real("battery/") && BE$an_is_real("battery/./") && BE$an_is_real("battery/_review/../") &&
        BE$an_is_real(edge_battery()))

fwrite(RES, file.path(RV, "review_indep_checks.csv"))
cat(sprintf("\n%d checks, %d differ\n", nrow(RES), sum(!RES$ok)))
if (any(!RES$ok)) print(RES[!RES$ok, ], row.names = FALSE)
sink()
