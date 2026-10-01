## analysis_verify2_impl.R -- part B: a second implementation of Section A, H1-H5, H6 and the paired tests in the reading of
## E12 AND E13 (nominal standard error; no mean over fewer cells; the declined-null rule; the size-failure columns and the
## final verdict; the H4 Holm family), written from the text of the pre-declaration, plus a comparator that runs
## analyse_M_battery.R on the same root and compares every output, value by value. No real battery result is read.
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

source(file.path(edge_path("code/simulations"), "analysis_verify2_lib.R"))

BE <- new.env()
sys.source(file.path(SIMDIR, "analyse_M_battery.R"), envir = BE)

EPS2 <- 1e-9
RES <- data.frame()
SET <- ""
check <- function(label, ok, detail = "") {
  RES <<- rbind(RES, data.frame(set = SET, check = label, ok = isTRUE(ok), detail = detail, stringsAsFactors = FALSE))
  cat(sprintf("  [%s] %s  %s\n", if (isTRUE(ok)) "ok" else "DIFF", label, detail))
}
eqn <- function(a, b, tol = 1e-10) length(a) == length(b) && all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= tol))
eqc <- function(a, b) length(a) == length(b) && all((is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b))
mean_all <- function(v) if (length(v) && all(is.finite(v))) mean(v) else NA_real_      # E13.2: no mean over fewer cells
cellstr <- function(block, cell) paste(sprintf("%s/%s", block, cell), collapse = " ")

## ---- the numbers a run uses (E12.1) and the declined-null rule (E13.2) ---------------------------------------------------------
v2_values <- function(root) {
  S <- rbindlist(lapply(c("0", "1a", "1b", "2", "3", "4", "5", "6", "7"), function(b) {
    f <- file.path(root, b, "_summary.csv")
    if (file.exists(f)) fread(f, colClasses = list(character = c("block", "cell", "role", "null_cell", "subset", "test", "status"))) else NULL
  }), use.names = TRUE, fill = TRUE)
  S[is.na(status), status := ""]
  A <- S[subset == "all" & !is.na(alpha) & !startsWith(test, "flag.")]
  A[startsWith(status, "not run"), `:=`(rejection = NA_real_, mcse = NA_real_, size_adj_power = NA_real_)]
  A[, k4 := paste(block, cell, test, sprintf("%.3f", alpha))]
  nl <- v2_null(A$block, A$cell)
  A[, `:=`(nb = nl$null_block, nc = nl$null_cell)]
  A[, nd := A$declined[match(paste(nb, nc, test, sprintf("%.3f", alpha)), A$k4)]]
  A[role == "null" | !nzchar(nc), nd := NA_real_]
  A[is.finite(size_adj_power) & is.finite(nd) & nd >= 1 - alpha - EPS2,
    `:=`(size_adj_power = NA_real_, status = "no value: matched null declined (E13.2)")]
  list(S = S, A = A)
}
g <- function(A, bl, ce, tt, col, a = 0.05) A[[col]][match(paste(bl, ce, tt, sprintf("%.3f", a)), A$k4)]
## C1 / E12.6 / E13.1: size at 0.05 in the matched null at most 0.05 + 3 sqrt(0.05 x 0.95 / B_null)
v2_holds <- function(A, block, cell, tt) {
  n0 <- v2_null(block, cell)
  s <- g(A, n0$null_block, n0$null_cell, tt, "rejection")
  B <- g(A, n0$null_block, n0$null_cell, tt, "B")
  ok <- nzchar(n0$null_cell) & is.finite(s) & is.finite(B) & s <= 0.05 + 3 * sqrt(0.05 * 0.95 / B) + EPS2
  ok[is.na(ok)] <- FALSE
  ok
}

## ---- Section A ----------------------------------------------------------------------------------------------------------------
v2_gate <- function(A, arm) {
  rbindlist(lapply(c("poly3", "sym", "stk"), function(b) rbindlist(lapply(c("u", "sc"), function(f) {
    tt <- sprintf("EDGE.%s.%s.%s", b, f, arm)
    nc <- unique(A[role == "null" & test == tt & is.finite(rejection) & alpha %in% c(0.01, 0.05), .(block, cell)])
    if (!nrow(nc)) return(NULL)
    one <- function(a) {
      r <- g(A, nc$block, nc$cell, tt, "rejection", a); B <- g(A, nc$block, nc$cell, tt, "B", a); m <- g(A, nc$block, nc$cell, tt, "mcse", a)
      list(r = r, mn = sqrt(a * (1 - a) / B), mr = m,
           okn = r <= a + 3 * sqrt(a * (1 - a) / B) + EPS2, okr = is.finite(m) & r <= a + 3 * m + EPS2)
    }
    a5 <- one(0.05); a1 <- one(0.01)
    data.table(basis = b, form = f, test = tt, block = nc$block, cell = nc$cell,
               r5 = a5$r, r1 = a1$r, mcse_null_05 = a5$mn, mcse_null_01 = a1$mn,
               z_null_05 = (a5$r - 0.05) / a5$mn, z_null_01 = (a1$r - 0.01) / a1$mn,
               z_real_05 = (a5$r - 0.05) / a5$mr, z_real_01 = (a1$r - 0.01) / a1$mr,
               ok5 = a5$okn, ok1 = a1$okn, ok_real5 = a5$okr, ok_real1 = a1$okr)
  }))))
}
v2_gate_pass <- function(G) {
  G <- copy(G)
  G[, `:=`(pass = !(ok5 %in% FALSE) & !(ok1 %in% FALSE), pass_real = !(ok_real5 %in% FALSE) & !(ok_real1 %in% FALSE))]
  G
}
v2_fam <- function(A, arm) {
  rbindlist(lapply(c("poly3", "sym", "stk"), function(b) rbindlist(lapply(c("u", "sc"), function(f) {
    p <- g(A, FAMS$block, FAMS$cell, sprintf("EDGE.%s.%s.%s", b, f, arm), "size_adj_power")
    fm <- rbindlist(lapply(1:6, function(k) data.table(family = k, cells = sum(FAMS$family == k), with_value = sum(is.finite(p[FAMS$family == k])),
                                                       mean = mean_all(p[FAMS$family == k]),
                                                       without = cellstr(FAMS$block[FAMS$family == k & !is.finite(p)], FAMS$cell[FAMS$family == k & !is.finite(p)]))))
    rbind(fm, data.table(family = 0L, cells = nrow(FAMS), with_value = sum(is.finite(p)), mean = mean_all(fm$mean),
                         without = cellstr(FAMS$block[!is.finite(p)], FAMS$cell[!is.finite(p)])))[, `:=`(basis = b, form = f)][]
  }))))
}
v2_dec <- function(G, Fm) {
  rbindlist(lapply(c("poly3", "sym", "stk"), function(b) {
    eu <- all(G[basis == b & form == "u"]$pass); es <- all(G[basis == b & form == "sc"]$pass)
    if (!nrow(G[basis == b & form == "u"])) eu <- FALSE
    if (!nrow(G[basis == b & form == "sc"])) es <- FALSE
    fu <- Fm[basis == b & form == "u" & family %in% 1:6][order(family)]$mean
    fs <- Fm[basis == b & form == "sc" & family %in% 1:6][order(family)]$mean
    mu <- mean_all(fu); ms <- mean_all(fs); d <- ms - mu
    ls <- if (all(is.finite(fu - fs))) max(fu - fs) else NA_real_       # what score loses against unit, per family
    lu <- if (all(is.finite(fs - fu))) max(fs - fu) else NA_real_
    if (eu && es) {
      if (!is.finite(d)) { ch <- "no decision" }
      else if (d >= 0.014 - EPS2) ch <- if (is.finite(ls) && ls > 0.10 + EPS2) "unit" else "score"
      else if (-d >= 0.014 - EPS2) ch <- if (is.finite(lu) && lu > 0.10 + EPS2) "score" else "unit"
      else ch <- "unit"
    } else if (eu) ch <- "unit" else if (es) ch <- "score" else ch <- "unit"
    data.table(basis = b, elig_u = eu, elig_s = es, macro_u = mu, macro_s = ms, diff = d, loss_s = ls, loss_u = lu, choice = ch)
  }))
}

## ---- H4 census (E12.5, the rule of make_headline_recount.R) ----------------------------------------------------------------------
v2_census <- function(A, ver, f) {
  edge <- sprintf("EDGE.poly3.%s.%s", f, ver)
  part <- c(paste0("HLF.", ver), paste0("HL.", ver), "HL_w", "PH"); six <- c(part, "Tsiatis", "Xie")
  rbindlist(lapply(seq_len(nrow(H4CELLS)), function(i) {
    gg <- function(tt) g(A, H4CELLS$block[i], H4CELLS$cell[i], tt, "size_adj_power")
    p3 <- gg(edge); rp <- vapply(part, gg, 0); r6 <- vapply(six, gg, 0)
    big <- function(v) if (any(!is.na(v))) max(v, na.rm = TRUE) else NA_real_
    low <- function(v) if (any(!is.na(v))) min(v, na.rm = TRUE) else NA_real_
    data.table(block = H4CELLS$block[i], cell = H4CELLS$cell[i], edge = edge, p3 = p3, best = big(rp),
               best_name = if (any(!is.na(rp))) names(which.max(rp)) else NA_character_,
               best6 = big(r6), best6_name = if (any(!is.na(r6))) names(which.max(r6)) else NA_character_,
               detect = big(c(p3, r6)) >= 0.15, sat = low(c(p3, rp)) >= 0.97,
               lead = p3 > big(rp) - 0.014, lead6 = p3 > big(r6) - 0.014, incomplete = anyNA(c(p3, r6)))
  }))[, headline := detect & !sat][]
}

## ---- H1-H5 with E13.2 and E13.3 ---------------------------------------------------------------------------------------------------
v2_hyp <- function(A, arm, forms) {
  H <- list(); PR <- list(); CEN <- list()
  F1 <- FAMS[family == 1]; F2 <- FAMS[family == 2]
  sap <- function(d, tt) g(A, d$block, d$cell, tt, "size_adj_power")
  vd <- function(st, ok) if (!is.finite(st)) "no value" else if (ok(st)) "holds" else "against"
  fin <- function(v, vw) if (is.na(vw) || v %in% c("no value", "reported") || identical(v, vw)) v else "unresolved (size)"
  dsgn <- function(d) ifelse(d$block == "2", sub("^[a-z0-9]+_(base|auc|e12)_n[0-9]+$", "\\1", d$cell), "base")
  mrow <- function(hyp, fn, used, d, ta, tb, dif, ok) {
    sf <- !(v2_holds(A, d$block, d$cell, ta) & v2_holds(A, d$block, d$cell, tb))
    st <- mean_all(dif); sw <- mean_all(dif[!sf]); v <- vd(st, ok); vw <- vd(sw, ok)
    data.table(hypothesis = hyp, arm = arm, form = fn, forms_used = used, cells = length(dif), cells_with_value = sum(is.finite(dif)),
               statistic = st, verdict = v, size_fail_cells = sum(sf), size_fail_cell_names = cellstr(d$block[sf], d$cell[sf]),
               statistic_without = sw, verdict_without = vw, final_verdict = fin(v, vw),
               cells_without_value = cellstr(d$block[!is.finite(dif)], d$cell[!is.finite(dif)]))
  }
  addp <- function(hyp, a, fn, used, d, ta, tb, pa, pb, inh = NA)
    PR[[length(PR) + 1L]] <<- data.table(hypothesis = hyp, arm = a, form = fn, forms_used = used, block = d$block, cell = d$cell,
                                         test_a = ta, test_b = tb, power_a = pa, power_b = pb, diff = pa - pb, in_headline = inh)
  for (f in c("u", "sc")) {
    fn <- c(u = "unit", sc = "score")[[f]]
    e3 <- sprintf("EDGE.poly3.%s.%s", f, arm); es <- sprintf("EDGE.sym.%s.%s", f, arm)
    u3 <- identical(unname(forms[["poly3"]]), f); us <- identical(unname(forms[["sym"]]), f)
    ## H1
    pa <- sap(F1, e3); pb <- sap(F1, "GiViTI"); d1 <- pa - pb; dg <- dsgn(F1)
    sf <- !(v2_holds(A, F1$block, F1$cell, e3) & v2_holds(A, F1$block, F1$cell, "GiViTI"))
    h1v <- function(keep) {
      st <- mean_all(d1[keep]); dm <- vapply(c("base", "auc", "e12"), function(z) mean_all(d1[keep & dg == z]), 0)
      below <- names(dm)[is.finite(dm) & dm < -0.014 - EPS2]
      list(st = st, dm = dm, below = below, v = if (!is.finite(st)) "no value" else if (st >= 0.10 - EPS2 && !length(below)) "holds" else "against")
    }
    aa <- h1v(rep(TRUE, nrow(F1))); ww <- h1v(!sf)
    H[[length(H) + 1L]] <- data.table(hypothesis = "H1", arm = arm, form = fn, forms_used = u3, cells = nrow(F1),
                                      cells_with_value = sum(is.finite(d1)), statistic = aa$st, verdict = aa$v,
                                      size_fail_cells = sum(sf), size_fail_cell_names = cellstr(F1$block[sf], F1$cell[sf]),
                                      statistic_without = ww$st, verdict_without = ww$v, final_verdict = fin(aa$v, ww$v),
                                      cells_without_value = cellstr(F1$block[!is.finite(d1)], F1$cell[!is.finite(d1)]),
                                      mean_base = aa$dm[["base"]], mean_auc = aa$dm[["auc"]], mean_e12 = aa$dm[["e12"]],
                                      below = paste(aa$below, collapse = " "))
    addp("H1", arm, fn, u3, F1, e3, "GiViTI", pa, pb)
    ## H2
    pa <- sap(F1, es); pb <- sap(F1, "Stk.joint")
    H[[length(H) + 1L]] <- mrow("H2", fn, us, F1, es, "Stk.joint", pa - pb, function(s) s >= -0.014 - EPS2)
    addp("H2", arm, fn, us, F1, es, "Stk.joint", pa, pb)
    ## H3
    pa <- sap(F2, e3); pb <- sap(F2, "GiViTI")
    H[[length(H) + 1L]] <- mrow("H3", fn, u3, F2, e3, "GiViTI", pa - pb, function(s) s >= -0.05 - EPS2)
    addp("H3", arm, fn, u3, F2, e3, "GiViTI", pa, pb)
    ## H4 at G = 10 and the rule-G version beside it
    for (ver in c("G10", "Grule")) {
      K <- v2_census(A, ver, f); hyp <- if (ver == "G10") "H4" else "H4.ruleG"
      edge <- sprintf("EDGE.poly3.%s.%s", f, ver)
      sf <- !(v2_holds(A, K$block, K$cell, rep(edge, nrow(K))) & v2_holds(A, K$block, K$cell, K$best_name))
      nov <- any(K$incomplete); nh <- sum(K$headline %in% TRUE)
      led <- sum(K$headline %in% TRUE & K$lead %in% TRUE); ledw <- sum(K$headline %in% TRUE & K$lead %in% TRUE & !sf)
      v  <- if (nov) "no value" else if (ver == "G10") (if (led >= 19) "holds" else "against") else "reported"
      vw <- if (nov) "no value" else if (ver == "G10") (if (ledw >= 19) "holds" else "against") else "reported"
      H[[length(H) + 1L]] <- data.table(hypothesis = hyp, arm = ver, form = fn, forms_used = u3, cells = nrow(K),
                                        cells_with_value = sum(!K$incomplete), statistic = if (nov) NA_real_ else led, verdict = v,
                                        size_fail_cells = sum(sf), size_fail_cell_names = cellstr(K$block[sf], K$cell[sf]),
                                        statistic_without = if (nov) NA_real_ else ledw, verdict_without = vw, final_verdict = fin(v, vw),
                                        cells_without_value = cellstr(K$block[K$incomplete], K$cell[K$incomplete]),
                                        headline_cells = nh, led_or_tied = led)
      CEN[[length(CEN) + 1L]] <- cbind(data.table(version = ver, form = fn), K)
      addp(hyp, ver, fn, u3, K, rep(edge, nrow(K)), K$best_name, K$p3, K$best, K$headline)
    }
    ## H5 and the skew reading
    pa <- sap(H5B, "Stk.sym1"); pb <- sap(H5B, es)
    H[[length(H) + 1L]] <- mrow("H5", fn, us, H5B, "Stk.sym1", es, pa - pb, function(s) s <= 0.02 + EPS2)
    addp("H5", arm, fn, us, H5B, "Stk.sym1", es, pa, pb)
    qa <- sap(H5S, "Stk.sym1"); qb <- sap(H5S, es); dq <- qa - qb
    sfq <- !(v2_holds(A, H5S$block, H5S$cell, "Stk.sym1") & v2_holds(A, H5S$block, H5S$cell, es))
    stq <- mean_all(dq); swq <- mean_all(dq[!sfq])
    H[[length(H) + 1L]] <- data.table(hypothesis = "H5.skew", arm = arm, form = fn, forms_used = us, cells = nrow(H5S),
                                      cells_with_value = sum(is.finite(dq)), statistic = stq,
                                      verdict = if (is.finite(stq)) "reported" else "no value", size_fail_cells = sum(sfq),
                                      size_fail_cell_names = cellstr(H5S$block[sfq], H5S$cell[sfq]), statistic_without = swq,
                                      verdict_without = if (is.finite(swq)) "reported" else "no value",
                                      final_verdict = if (is.finite(stq)) "reported" else "no value",
                                      cells_without_value = cellstr(H5S$block[!is.finite(dq)], H5S$cell[!is.finite(dq)]))
    addp("H5.skew", arm, fn, us, H5S, "Stk.sym1", es, qa, qb)
  }
  list(hyp = rbindlist(H, fill = TRUE), pairs = rbindlist(PR, fill = TRUE), census = rbindlist(CEN))
}

## ---- the pairs a run needs (the refusal list of E13.2) -------------------------------------------------------------------------------
v2_missing <- function(A, arm) {
  ef <- function(b, gg) sprintf("EDGE.%s.%s.%s", b, c("u", "sc"), gg)
  one <- function(set, d, tests) data.table(set = set, block = rep(d$block, each = length(tests)), cell = rep(d$cell, each = length(tests)),
                                            test = rep(tests, nrow(d)))
  six <- c("HLF.G10", "HL.G10", "HL_w", "PH", "Tsiatis", "Xie", "HLF.Grule", "HL.Grule")
  N <- rbind(one("Section A families", FAMS, unlist(lapply(unique(c(arm, "Grule")), function(gg) unlist(lapply(c("poly3", "sym", "stk"), ef, gg = gg))))),
             one("H1", FAMS[family == 1], c(ef("poly3", arm), "GiViTI")),
             one("H2", FAMS[family == 1], c(ef("sym", arm), "Stk.joint")),
             one("H3", FAMS[family == 2], c(ef("poly3", arm), "GiViTI")),
             one("H4", H4CELLS, unique(c(ef("poly3", "G10"), ef("poly3", "Grule"), six))),
             one("H5", H5B, c("Stk.sym1", ef("sym", arm))),
             one("H5.skew", H5S, c("Stk.sym1", ef("sym", arm))))
  N[!is.finite(g(A, block, cell, test, "size_adj_power"))]
}

## ---- H6: the declined rate of every test in every cell, overall and by subset ----------------------------------------------------------
v2_h6 <- function(S) {
  D <- S[!startsWith(test, "flag.") & !is.na(alpha) & abs(alpha - 0.05) < 1e-9,
         .(block, cell, test, subset, B, declined, rejection_05 = rejection, source = "summary")]
  Fl <- S[startsWith(test, "flag.evlt.") & subset == "all" & rejection %in% c(0, 1),
          .(block, cell, G = sub("flag.evlt.", "", test, fixed = TRUE), rate = rejection)]
  have <- unique(D[subset != "all", .(block, cell, subset)])
  base <- D[subset == "all"][, subset := NULL]
  add <- function(sub_of, src, blank) {
    z <- Fl[, .(block, cell, subset = sub_of(rate, G))][!have, on = c("block", "cell", "subset")]
    if (!nrow(z)) return(NULL)
    y <- base[z, on = c("block", "cell"), nomatch = NULL, allow.cartesian = TRUE][, source := src]
    if (blank) y[, `:=`(B = 0, declined = NA_real_, rejection_05 = NA_real_)]
    y
  }
  a1 <- add(function(r, G) paste0(ifelse(r == 1, "events<", "events>="), G), "all samples are in this subset (flag rate 0 or 1)", FALSE)
  a0 <- add(function(r, G) paste0(ifelse(r == 1, "events>=", "events<"), G), "no sample in this subset", TRUE)
  D <- rbind(D, a1, a0, use.names = TRUE)
  D[, arm := ifelse(grepl("\\.G([0-9]+|rule)$", test), sub("^.*\\.(G([0-9]+|rule))$", "\\1", test),
                    ifelse(test %in% c("HL_w", "PH", "Tsiatis", "Xie", "PR"), "G10", "none"))]
  D[, sub_G := ifelse(subset == "all", "", sub("^events(<|>=)", "", subset))]
  D[, matches := subset == "all" | arm == "none" | sub_G == arm]
  D[, key := paste(block, cell, test, subset)][]
}

## ---- the paired tests (C3, E12.6, E13.4) -------------------------------------------------------------------------------------------
holm_hand <- function(p) { m <- length(p); o <- order(p); run <- 0; adj <- numeric(m)
  for (j in seq_len(m)) { run <- max(run, min(1, (m - j + 1) * p[o[j]])); adj[o[j]] <- run }; adj }
mcn_hand <- function(b, c) { m <- b + c; if (m == 0) 1 else min(1, 2 * stats::pbinom(min(b, c), m, 0.5)) }
v2_mcnemar <- function(root, A, PR) {
  X <- copy(PR)
  n0 <- v2_null(X$block, X$cell); X[, `:=`(null_block = n0$null_block, null_cell = n0$null_cell)]
  gs <- function(tt, col) g(A, X$null_block, X$null_cell, tt, col)
  X[, `:=`(size_a = gs(test_a, "rejection"), Ba = gs(test_a, "B"), mcse_a = gs(test_a, "mcse"),
           size_b = gs(test_b, "rejection"), Bb = gs(test_b, "B"), mcse_b = gs(test_b, "mcse"))]
  X[, `:=`(lim_a = 0.05 + 3 * sqrt(0.05 * 0.95 / Ba), lim_b = 0.05 + 3 * sqrt(0.05 * 0.95 / Bb))]
  X[, `:=`(hold_a = is.finite(size_a) & is.finite(lim_a) & size_a <= lim_a + EPS2,
           hold_b = is.finite(size_b) & is.finite(lim_b) & size_b <= lim_b + EPS2,
           z_nom_a = (size_a - 0.05) / sqrt(0.05 * 0.95 / Ba), z_nom_b = (size_b - 0.05) / sqrt(0.05 * 0.95 / Bb),
           z_real_a = (size_a - 0.05) / mcse_a, z_real_b = (size_b - 0.05) / mcse_b)]
  X[, in_family := !startsWith(hypothesis, "H4") | in_headline %in% TRUE]      # E13.4
  X[, `:=`(both = NA_integer_, a_only = NA_integer_, b_only = NA_integer_, neither = NA_integer_, p = NA_real_, why = "")]
  cache <- list()
  for (i in seq_len(nrow(X))) {
    f <- file.path(root, X$block[i], paste0(X$cell[i], "_pvalues.csv.gz"))
    why <- character(0)
    if (!nzchar(X$null_cell[i])) why <- c(why, "no matched null")
    else { if (!X$hold_a[i]) why <- c(why, "a")
           if (is.na(X$test_b[i]) || !X$hold_b[i]) why <- c(why, "b") }
    if (!X$in_family[i]) why <- c(why, "E13.4")
    if (!file.exists(f)) why <- c(why, "no file")
    else {
      if (is.null(cache[[f]])) cache[[f]] <- fread(f)
      P <- cache[[f]]
      if (is.na(X$test_b[i]) || !all(c(X$test_a[i], X$test_b[i]) %in% names(P))) why <- c(why, "no column")
      else {
        ra <- !is.na(P[[X$test_a[i]]]) & P[[X$test_a[i]]] <= 0.05
        rb <- !is.na(P[[X$test_b[i]]]) & P[[X$test_b[i]]] <= 0.05
        set(X, i, c("both", "a_only", "b_only", "neither"), list(sum(ra & rb), sum(ra & !rb), sum(!ra & rb), sum(!ra & !rb)))
        if (!length(why)) set(X, i, "p", mcn_hand(sum(ra & !rb), sum(!ra & rb)))
      }
    }
    set(X, i, "why", paste(why, collapse = ";"))
  }
  X[, used := is.finite(p)]
  X[, p_holm := NA_real_]
  X[used == TRUE, p_holm := holm_hand(p), by = .(hypothesis, arm, form)]
  X[, key := paste(hypothesis, arm, form, block, cell)][]
}

## ---- run analyse_M_battery.R on the same root and compare every output ----------------------------------------------------------------
v2_run_builder <- function(args, bat = BE$BAT) {
  res <- NULL
  msg <- tryCatch({ txt <- capture.output(res <- BE$an_main(args, bat = bat)); "ran" }, error = function(e) conditionMessage(e))
  list(msg = msg, res = res)
}
v2_compare <- function(root, label, arm = "rule", bat = BE$BAT, extra = character(0)) {
  SET <<- paste(label, arm)
  cat("\n== ", label, ", arm ", arm, "  ", root, "\n", sep = "")
  armk <- if (arm == "rule") "Grule" else "G10"
  suf <- if (arm == "rule") "" else "_G10"
  out <- file.path(root, paste0("out_", arm))
  V <- v2_values(root); A <- V$A
  G <- v2_gate_pass(v2_gate(A, armk)); Fm <- v2_fam(A, armk)
  Gr <- if (armk == "Grule") G else v2_gate_pass(v2_gate(A, "Grule"))
  Fr <- if (armk == "Grule") Fm else v2_fam(A, "Grule")
  Dc <- v2_dec(G, Fm); Dr <- v2_dec(Gr, Fr)                       # the forms come from the rule-G decision (E12.2)
  forms <- setNames(c(unit = "u", score = "sc", "no decision" = NA_character_)[Dr$choice], Dr$basis)
  cat("  second implementation:", paste(sprintf("%s %s (diff %s)", Dc$basis, Dc$choice, ifelse(is.finite(Dc$diff), sprintf("%+.4f", Dc$diff), "NA")), collapse = "; "), "\n")
  B <- v2_run_builder(c("--root", root, "--out", out, if (arm == "10") c("--arm", "10"), extra), bat = bat)
  check("the builder runs", B$msg == "ran", B$msg)
  if (B$msg != "ran") return(invisible(list(msg = B$msg)))
  b <- B$res
  ## gate
  Gv <- copy(G); kk <- setdiff(names(Gv), c("test", "block", "cell")); setnames(Gv, kk, paste0("v_", kk))
  GG <- merge(b$gate[, .(test, block, cell, pass, pass_null, pass_realised, pass_05, pass_01, z_05, z_01, z_null_05, z_null_01,
                         z_realised_05, z_realised_01, mcse_null_05, mcse_null_01)],
              Gv, by = c("test", "block", "cell"), all = TRUE)
  check("gate: the same variant-cells", nrow(GG) == nrow(Gv) && nrow(GG) == nrow(b$gate) && !anyNA(GG$pass) && !anyNA(GG$v_pass),
        sprintf("%d rows, %d fail the nominal gate", nrow(GG), sum(!GG$v_pass)))
  check("gate: pass at 0.05, at 0.01 and overall, under the nominal standard error (E13.1)",
        eqc(GG$pass_05, GG$v_ok5) && eqc(GG$pass_01, GG$v_ok1) && eqc(GG$pass, GG$v_pass) && eqc(GG$pass_null, GG$v_pass))
  check("gate: the realised reading is reported and does not decide", eqc(GG$pass_realised, GG$v_pass_real),
        sprintf("nominal fails %d, realised fails %d", sum(!GG$v_pass), sum(!GG$v_pass_real)))
  check("gate: z and MCSE columns", eqn(GG$z_05, GG$v_z_null_05, 1e-9) && eqn(GG$z_null_05, GG$v_z_null_05, 1e-9) &&
          eqn(GG$z_null_01, GG$v_z_null_01, 1e-9) && eqn(GG$z_realised_05, GG$v_z_real_05, 1e-9) &&
          eqn(GG$z_realised_01, GG$v_z_real_01, 1e-9) && eqn(GG$mcse_null_05, GG$v_mcse_null_05, 1e-12))
  ## membership
  Mb <- b$M
  check("membership: the same cells and families as the typed lists",
        setequal(paste(Mb$block, Mb$cell, Mb$family), paste(FAMS$block, FAMS$cell, FAMS$family)) && nrow(Mb) == nrow(FAMS) &&
          identical(as.integer(table(FAMS$family)), c(32L, 21L, 10L, 80L, 13L, 2L)))
  ## families and the decision
  Fv <- copy(Fm); setnames(Fv, c("cells", "with_value", "mean", "without"), c("v_cells", "v_with_value", "v_mean", "v_without"))
  FB <- merge(b$families[, .(basis, form, family = ifelse(family == "macro", 0L, suppressWarnings(as.integer(family))), cells, with_value,
                             mean_power, mean_power_available, cells_without_value)],
              Fv, by = c("basis", "form", "family"), all = TRUE)
  check("family means and the macro-average: no mean over fewer cells (E13.2)", nrow(FB) == 42 && eqn(FB$mean_power, FB$v_mean) &&
          eqc(as.integer(FB$cells), as.integer(FB$v_cells)) && eqc(as.integer(FB$with_value), as.integer(FB$v_with_value)),
        sprintf("%d of %d means have no value", sum(is.na(FB$v_mean)), nrow(FB)))
  check("families: the cells without a value are named", eqc(FB$cells_without_value, FB$v_without))
  ## E13.2 says the analysis never averages over fewer cells: the available columns are reported beside the mean
  av <- b$families[is.na(mean_power) & is.finite(mean_power_available)]
  check("NOTE (finding 2): rule_A_families.csv still carries a mean over the cells with a value", TRUE,
        sprintf("%d rows have no mean but a mean_power_available", nrow(av)))
  DB <- merge(b$decision, Dc, by = "basis")
  check("eligibility, macro-averages, difference and the family losses",
        eqc(DB$eligible_unit, DB$elig_u) && eqc(DB$eligible_score, DB$elig_s) && eqn(DB$macro_unit, DB$macro_u) &&
          eqn(DB$macro_score, DB$macro_s) && eqn(DB$difference, DB$diff) && eqn(DB$worst_loss_score, DB$loss_s) &&
          eqn(DB$worst_loss_unit, DB$loss_u))
  check("the choice per basis", eqc(DB$choice.x, DB$choice.y), paste(DB$basis, DB$choice.y, collapse = ", "))
  check("the forms used for H1-H5 come from the rule-G decision",
        identical(unname(b$forms_used[c("poly3", "sym")]), unname(forms[c("poly3", "sym")])),
        paste(names(forms), ifelse(is.na(forms), "none", forms), collapse = " "))
  ## hypotheses
  mh <- v2_hyp(A, armk, forms)
  HB <- merge(b$hypotheses, mh$hyp, by = c("hypothesis", "arm", "form"), all = TRUE, suffixes = c(".b", ".v"))
  check("hypotheses: the same 14 rows", nrow(HB) == 14 && !anyNA(HB$verdict.b) && !anyNA(HB$verdict.v))
  check("hypotheses: statistic, verdict, cells and cells with a value (E13.2)",
        eqn(HB$statistic.b, HB$statistic.v) && eqc(HB$verdict.b, HB$verdict.v) && eqc(as.integer(HB$cells.b), as.integer(HB$cells.v)) &&
          eqc(as.integer(HB$cells_with_value.b), as.integer(HB$cells_with_value.v)),
        paste(sprintf("%s/%s %s", HB$hypothesis, substr(HB$form, 1, 1), HB$verdict.v), collapse = "; "))
  check("hypotheses: size failures, the statistic and verdict without them, the final verdict (E13.3)",
        eqc(as.integer(HB$size_fail_cells.b), as.integer(HB$size_fail_cells.v)) &&
          eqc(HB$size_fail_cell_names.b, HB$size_fail_cell_names.v) &&
          eqn(HB$statistic_without_size_failures, HB$statistic_without) &&
          eqc(HB$verdict_without_size_failures, HB$verdict_without) && eqc(HB$final_verdict.b, HB$final_verdict.v),
        sprintf("rows with a size failure %d; unresolved (size) %d", sum(HB$size_fail_cells.v > 0), sum(HB$final_verdict.v == "unresolved (size)")))
  check("hypotheses: the cells without a value are named", eqc(HB$cells_without_value.b, HB$cells_without_value.v))
  h1 <- HB[hypothesis == "H1"]
  check("H1: the three design means and the designs below -0.014 (no design mean over fewer cells)",
        eqn(h1$mean_base.b, h1$mean_base.v) && eqn(h1$mean_auc.b, h1$mean_auc.v) && eqn(h1$mean_e12.b, h1$mean_e12.v) &&
          eqc(h1$designs_below_margin, h1$below),
        paste(sprintf("%s base %s auc %s e12 %s", substr(h1$form, 1, 1), format(h1$mean_base.v, digits = 3), format(h1$mean_auc.v, digits = 3),
                      format(h1$mean_e12.v, digits = 3)), collapse = "; "))
  h4 <- HB[hypothesis %in% c("H4", "H4.ruleG")]
  check("H4: headline cells and cells led or tied", eqc(as.integer(h4$headline_cells.b), as.integer(h4$headline_cells.v)) &&
          eqc(as.integer(h4$led_or_tied.b), as.integer(h4$led_or_tied.v)),
        paste(sprintf("%s/%s %d of %d", h4$hypothesis, substr(h4$form, 1, 1), as.integer(h4$led_or_tied.v), as.integer(h4$headline_cells.v)), collapse = "; "))
  ## census
  KB <- merge(b$census[, .(version, form, block, cell, edge_poly3, best_partition, best_partition_name, best_all6, best_all6_name,
                           detectable, saturated, beats_partition, beats_all6, in_headline)],
              mh$census, by.x = c("version", "form", "block", "cell"), by.y = c("version", "form", "block", "cell"), all = TRUE)
  check("census cell by cell: detectable, saturated, leads or ties, best partition and best of six, headline",
        nrow(KB) == 4 * nrow(H4CELLS) && eqn(KB$edge_poly3, KB$p3) && eqn(KB$best_partition, KB$best) &&
          eqc(KB$best_partition_name, KB$best_name) && eqn(KB$best_all6, KB$best6) && eqc(KB$best_all6_name, KB$best6_name) &&
          eqc(KB$detectable, KB$detect) && eqc(KB$saturated, KB$sat) && eqc(KB$beats_partition, KB$lead) &&
          eqc(KB$beats_all6, KB$lead6) && eqc(KB$in_headline, KB$headline), sprintf("%d rows", nrow(KB)))
  ## the refusal list
  MV <- v2_missing(A, armk)
  check("missing_values: the same (set, cell, test) pairs without a size-adjusted power",
        setequal(paste(b$missing_values$set, b$missing_values$block, b$missing_values$cell, b$missing_values$test),
                 paste(MV$set, MV$block, MV$cell, MV$test)) && nrow(b$missing_values) == nrow(MV),
        sprintf("builder %d, second implementation %d", nrow(b$missing_values), nrow(MV)))
  ## H6
  D6 <- v2_h6(V$S)
  BD <- copy(b$declined)[, key := paste(block, cell, test, subset)]
  MD <- merge(BD, D6, by = "key", all = TRUE, suffixes = c(".b", ".v"))
  check("H6 declined: the same rows, B, declined rate and subset arithmetic",
        nrow(MD) == nrow(D6) && nrow(MD) == nrow(BD) && eqn(as.numeric(MD$B.b), as.numeric(MD$B.v)) &&
          eqn(MD$declined.b, MD$declined.v) && eqn(MD$rejection_05.b, MD$rejection_05.v) &&
          eqc(MD$test_arm, MD$arm) && eqc(MD$subset_matches_arm, MD$matches) && eqc(MD$source.b, MD$source.v),
        sprintf("%d rows, %d from a flag rate of 0 or 1", nrow(MD), sum(MD$source.v != "summary")))
  ## paired tests
  mm <- v2_mcnemar(root, A, mh$pairs)
  MB <- copy(b$mcnemar)[, key := paste(hypothesis, arm, form, block, cell)]
  MC <- merge(MB, mm, by = "key", all = TRUE, suffixes = c(".b", ".v"))
  check("paired tests: the same rows and the same test pairs", nrow(MC) == nrow(mm) && nrow(MC) == nrow(MB) &&
          eqc(MC$test_a.b, MC$test_a.v) && eqc(MC$test_b.b, MC$test_b.v), sprintf("%d rows", nrow(MC)))
  check("paired tests: size held under the nominal standard error (E13.1), the family of E13.4, and which rows are compared",
        eqc(MC$holds_a, MC$hold_a) && eqc(MC$holds_b, MC$hold_b) && eqc(MC$in_family.b, MC$in_family.v) && eqc(MC$used.b, MC$used.v),
        sprintf("in family %d of %d, compared %d", sum(MC$in_family.v), nrow(MC), sum(MC$used.v)))
  check("paired tests: the nominal and realised z, reported beside the check",
        eqn(MC$z_nominal_a, MC$z_nom_a, 1e-9) && eqn(MC$z_nominal_b, MC$z_nom_b, 1e-9) && eqn(MC$z_realised_a, MC$z_real_a, 1e-9) &&
          eqn(MC$z_realised_b, MC$z_real_b, 1e-9) && eqn(MC$limit_a, MC$lim_a, 1e-12) && eqn(MC$limit_b, MC$lim_b, 1e-12))
  cnt <- MC[!is.na(both.b) | !is.na(both.v)]
  check("paired tests: the 2 x 2 counts on the shared replicates", eqn(as.numeric(cnt$both.b), as.numeric(cnt$both.v)) &&
          eqn(as.numeric(cnt$a_only.b), as.numeric(cnt$a_only.v)) && eqn(as.numeric(cnt$b_only.b), as.numeric(cnt$b_only.v)) &&
          eqn(as.numeric(cnt$neither.b), as.numeric(cnt$neither.v)), sprintf("%d rows with a file and both columns", nrow(cnt)))
  check("paired tests: the exact McNemar p and Holm within each hypothesis, arm and form (E13.4)",
        eqn(MC$p_mcnemar, MC$p, 1e-12) && eqn(MC$p_holm.b, MC$p_holm.v, 1e-12),
        paste(MC[used.v == TRUE, .(m = .N, sig = sum(p_holm.v <= 0.05)), by = .(hypothesis.v, arm.v, form.v)][, sprintf("%s/%s %d of %d", hypothesis.v, substr(form.v, 1, 1), sig, m)], collapse = "; "))
  ## the files on disk
  fs <- c(paste0(c("rule_A_gate", "rule_A_membership", "rule_A_families", "hypotheses", "h4_census", "declined", "mcnemar",
                   "missing_values"), suf, ".csv"), paste0("rule_A_decision", suf, ".txt"))
  check("the nine outputs are written", all(file.exists(file.path(out, fs))), paste(fs[!file.exists(file.path(out, fs))], collapse = " "))
  fh <- fread(file.path(out, paste0("hypotheses", suf, ".csv")))
  fk <- fread(file.path(out, paste0("mcnemar", suf, ".csv")))
  check("the written files hold what the run returned", nrow(fh) == nrow(b$hypotheses) && nrow(fk) == nrow(b$mcnemar) &&
          eqn(fh$statistic, b$hypotheses$statistic, 1e-9) && eqc(fh$final_verdict, b$hypotheses$final_verdict) &&
          eqn(fk$p_holm, b$mcnemar$p_holm, 1e-9))
  dl <- readLines(file.path(out, paste0("rule_A_decision", suf, ".txt")))
  check("rule_A_decision names the nominal reading and ends with the three choices",
        any(grepl("^MCSE of the size gate: null", dl)) && all(sprintf("choice %s %s", Dc$basis, Dc$choice) %in% dl))
  invisible(list(A = A, mine = list(gate = G, fam = Fm, dec = Dc, hyp = mh$hyp, census = mh$census, mcnemar = mm, missing = MV, h6 = D6),
                 builder = b, merged = list(hyp = HB, mcnemar = MC, fam = FB)))
}
v2_write_checks <- function(file) {
  fwrite(RES, file)
  cat(sprintf("\n%d checks, %d differ\n", nrow(RES), sum(!RES$ok)))
  if (any(!RES$ok)) print(RES[!RES$ok, ], row.names = FALSE)
}
