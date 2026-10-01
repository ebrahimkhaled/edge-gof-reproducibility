## _battery_cells.R -- the cell table of the EDGE restructure battery (PREDECLARATION_restructure_battery.md, E3).
## battery_cells() returns one row per cell: block, cell name, scenario, design, link, n, G list, B, seed root, cell id,
## seed base, generator, fitted formula, the artifact it reproduces, and its matched null.
##
## Seeds (E0.4): cells of the July grids (blocks 1a, 3, 4) keep their old seed roots and cell ids, with the ids rebuilt
## from each old script's own enumeration code below (copied, not retyped), so seed_base = root + cell_id * 1e5 as before.
## New cells use seed_base = 100000000 + block * 1e7 + cell_id * 1e4, with block 1b counted as 1. In blocks that also hold
## old cells (3, 4) the new ids start at 201 so an id never names two cells of one block.
## verify_old_ids() checks the rebuilt ids against the seeds stored in the July p-value files; seed_overlaps() lists
## every pair of cells whose seed ranges [seed_base + 1, seed_base + B] intersect.
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

BC_SIMDIR <- edge_path("code/simulations")
BC_BLOCKS <- c("0", "1a", "1b", "2", "3", "4", "5", "6", "7")
BC_BLOCKNUM <- c("0" = 0, "1a" = 1, "1b" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5, "6" = 6, "7" = 7)

## ---- old enumerations, copied from the July scripts ---------------------------------------------------------------
old_cells_grid_null <- function() {                                # grid_null.R
  FAMILIES <- c("quad", "binint", "contint", "link", "rough")
  N_MAIN   <- c(200L, 500L, 1000L, 2000L, 5000L)
  main_grid <- expand.grid(family = FAMILIES, n = N_MAIN, G = 10L, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  gsens_grid <- expand.grid(family = c("link", "quad"), n = 1000L, G = c(6L, 8L, 12L, 14L, 20L),
                            stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  cells <- rbind(main_grid, gsens_grid)
  cells$cell_id <- seq_len(nrow(cells))
  cells$B <- 10000L; cells$seed_root <- 20260705L
  cells
}

old_cells_grid_power_broad <- function() {                         # grid_power_broad.R
  G_DEFAULT <- 10L; GSENS_GVALS <- c(6L, 8L, 10L, 12L, 14L, 20L)
  N_LADDER <- c(200L, 500L, 1000L, 2000L, 5000L)
  build_family_cells <- function(family, params, n_grid) {
    expand.grid(family = family, param = as.character(params), n = n_grid, G = G_DEFAULT,
                stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  }
  cells_quad    <- build_family_cells("quad",    c(0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.40), N_LADDER)
  cells_binint  <- build_family_cells("binint",  c(0.1, 0.2, 0.3, 0.5, 0.7),                   N_LADDER)
  cells_contint <- build_family_cells("contint", c(0.1, 0.3, 0.5, 0.7),                        N_LADDER)
  cells_link    <- build_family_cells("link",    c("probit", "cloglog", "stukel_heavy", "stukel_light", "stukel_asym"), N_LADDER)
  cells_rough   <- build_family_cells("rough",   c("osc2", "osc4", "sawtooth", "bump"),         c(500L, 1000L, 2000L))
  cells <- rbind(cells_quad, cells_binint, cells_contint, cells_link, cells_rough)
  gsens_cells <- rbind(
    build_family_cells("link", "cloglog", 1000L)[rep(1, length(GSENS_GVALS)), , drop = FALSE],
    build_family_cells("quad", "0.05",    1000L)[rep(1, length(GSENS_GVALS)), , drop = FALSE])
  gsens_cells$G <- rep(GSENS_GVALS, times = 2)
  gsens_cells$is_gsens <- TRUE
  cells$is_gsens <- FALSE
  all_cells <- rbind(cells, gsens_cells)
  all_cells$cell_id <- seq_len(nrow(all_cells))
  all_cells$B <- ifelse(all_cells$is_gsens, 5000L, ifelse(all_cells$n == 200L, 10000L, 5000L))
  all_cells$seed_root <- 20260708L
  all_cells
}

old_cells_edge_loses <- function() {                               # grid_edge_loses.R (A_delta column left out)
  cells <- rbind(
    data.frame(scenario = "crossover", shape = NA_character_, n = c(1000L, 2000L), G = 10L, B = 5000L, seed_root = 20260710L, stringsAsFactors = FALSE),
    data.frame(scenario = "osc4", shape = "osc4", n = c(1000L, 2000L), G = 10L, B = 5000L, seed_root = 20260711L, stringsAsFactors = FALSE),
    data.frame(scenario = "sawtooth", shape = "sawtooth", n = c(1000L, 2000L), G = 10L, B = 5000L, seed_root = 20260711L, stringsAsFactors = FALSE))
  cells$cell_id <- seq_len(nrow(cells))
  cells
}

old_cells_proj <- function() {                                     # grid_proj_power.R (loop order: cells x NS)
  cells <- list(c("link", "cloglog"), c("link", "probit"), c("link", "stukel_heavy"), c("link", "stukel_light"),
                c("link", "stukel_asym"), c("quad", "0.02"), c("binint", "0.3"), c("contint", "0.5"),
                c("omit_2cov", "off-index"), c("null_quad", "null"), c("null_link", "null"), c("null_contint", "null"))
  NS <- c(500L, 1000L)
  out <- list(); cell_id <- 0L
  for (cell0 in cells) for (n in NS) {
    cell_id <- cell_id + 1L
    out[[cell_id]] <- data.frame(scenario = cell0[1], param = cell0[2], n = n, cell_id = cell_id, B = 500L,
                                 seed_root = 20260720L, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

## ---- one row ---------------------------------------------------------------------------------------------------------------
bc_cell <- function(block, cell, generator, n, B, formula, reproduces, scenario = cell, design = NA, link = NA,
                    G_extra = "", seed_root = NA, cell_id = NA, family = NA, param = NA, s = NA, c0 = NA, xdist = NA,
                    intercept = NA, slope = NA, scen = NA, pstar = NA, kind = NA, dist_s = NA, dataset = NA,
                    type = "full", ao = FALSE, null_block = NA, null_cell = NA, source = "", pair_proj = FALSE,
                    stored = NA, stored_key = NA, orders = NA_integer_, notes = "") {
  data.frame(block = block, cell = cell, scenario = scenario, design = design, link = link, n = n, G_extra = G_extra,
             B = B, orders = orders, seed_root = seed_root, cell_id = cell_id, generator = generator, family = family, param = param,
             s = s, c0 = c0, xdist = xdist, intercept = intercept, slope = slope, scen = scen, pstar = pstar, kind = kind,
             dist_s = dist_s, dataset = dataset, type = type, ao = ao, formula = formula, reproduces = reproduces,
             null_block = null_block, null_cell = null_cell, source = source, pair_proj = pair_proj, stored = stored,
             stored_key = stored_key, notes = notes, stringsAsFactors = FALSE)
}
PRE_LAUNCH <- "pre-launch addition after review (not in E3)"

bc_read_nplan <- function(path = edge_battery("nplan.csv")) {
  if (!file.exists(path)) stop("battery/nplan.csv is missing: run battery_nplan.R first")
  np <- utils::read.csv(path, stringsAsFactors = FALSE)
  np <- np[np$in_block2, ]
  np$ns <- lapply(np$n_run, function(z) as.integer(strsplit(z, ",")[[1]]))
  np
}

PLATEAU_NAME <- c(stukel_heavy = "plateau_upper", stukel_light = "plateau_lower", stukel_asym = "plateau_both")
FAM_FORMULA  <- c(quad = "y ~ x", binint = "y ~ x + d", contint = "y ~ x + z", link = "y ~ x + d", rough = "y ~ x")
CENSUS <- c("logx", "int_binbin", "skew", "corr", "omit_x2x3", "omit_2int", "joint", "omit_2cov")
CENSUS_FORMULA <- c(logx = "y ~ x", int_binbin = "y ~ d1 + d2", skew = "y ~ x (x ~ chi-square(4))", corr = "y ~ x1 + x2",
                    omit_x2x3 = "y ~ x", omit_2int = "y ~ x + d + z", joint = "y ~ x", omit_2cov = "y ~ x")
PROJ_FORMULA <- c(link = "y ~ x + d", quad = "y ~ x", binint = "y ~ x + d", contint = "y ~ x + z", omit_2cov = "y ~ x",
                  null_quad = "y ~ x", null_link = "y ~ x + d", null_contint = "y ~ x + z")

## ---- the table ---------------------------------------------------------------------------------------------------------------
battery_cells <- function(nplan = bc_read_nplan()) {
  R <- list(); add <- function(df) R[[length(R) + 1]] <<- df

  gn <- old_cells_grid_null(); gp <- old_cells_grid_power_broad(); ge <- old_cells_edge_loses(); gj <- old_cells_proj()
  gn_row <- function(fam, n) gn[gn$family == fam & gn$n == n & gn$G == 10L, ]

  ## ---- block 0: validation --------------------------------------------------------------------------------------------------
  z <- gn_row("quad", 1000L)
  add(bc_cell("0", "id_null_quad_n1000", "dgp_null", 1000L, 200L, "y ~ x", "block 0 (i): identity with sim_null_pvalues.csv",
              family = "quad", seed_root = z$seed_root, cell_id = z$cell_id, source = "grid_null.R", stored = "sim_null"))
  z <- gp[gp$family == "link" & gp$param == "cloglog" & gp$n == 1000L & !gp$is_gsens, ]
  add(bc_cell("0", "id_link_cloglog_n1000", "dgp_alt", 1000L, 200L, "y ~ x + d", "block 0 (i): identity with sim_power_broad_pvalues.csv",
              family = "link", param = "cloglog", link = "cloglog", seed_root = z$seed_root, cell_id = z$cell_id,
              source = "grid_power_broad.R", stored = "sim_power_broad"))
  z <- ge[ge$scenario == "crossover" & ge$n == 1000L, ]
  add(bc_cell("0", "id_crossover_n1000", "crossover", 1000L, 200L, "y ~ x + d", "block 0 (i): identity with sim_edge_loses_pvalues.csv",
              seed_root = z$seed_root, cell_id = z$cell_id, source = "grid_edge_loses.R", stored = "sim_edge_loses"))
  add(bc_cell("0", "v_sparse49_n300", "sparse", 300L, 200L, "y ~ x", "block 0 (ii): harness vs package, one-df Stukel fallback and events < G",
              link = "logit", intercept = -4.9, slope = 1, cell_id = 1L))
  add(bc_cell("0", "v_probit_auc_n2000", "design", 2000L, 200L, "y ~ x + d", "block 0 (ii): harness vs package at high AUC (score form)",
              design = "auc", link = "probit", s = 2, c0 = 0, xdist = "uniform", cell_id = 2L))
  z <- gn_row("quad", 5000L)
  add(bc_cell("0", "w_null_quad_n5000", "dgp_null", 5000L, 200L, "y ~ x", "block 0 (iii): worker timing (E0.8)",
              family = "quad", seed_root = z$seed_root, cell_id = z$cell_id, source = "grid_null.R"))

  ## ---- block 1a: size grid, old seeds -------------------------------------------------------------------------------------------
  for (i in which(gn$G == 10L)) {
    z <- gn[i, ]
    rep_txt <- "Fig 7 size; Sec 6.1 size text; Sec 6.5 min-p size"
    if (z$family == "quad" && z$n == 1000L) rep_txt <- paste(rep_txt, "; Table 4 null row")
    if (z$family == "link" && z$n <= 1000L) rep_txt <- paste(rep_txt, "; tab:nulllevel (base-design null)")
    add(bc_cell("1a", sprintf("null_%s_n%d", z$family, z$n), "dgp_null", z$n, 10000L, FAM_FORMULA[[z$family]], rep_txt,
                scenario = paste0("null_", z$family), family = z$family, seed_root = z$seed_root, cell_id = z$cell_id,
                G_extra = if (z$family %in% c("quad", "link") && z$n == 1000L) "6,8,12,14,20" else "",
                source = "grid_null.R", stored = "sim_null"))
  }

  ## ---- block 1b: new nulls ------------------------------------------------------------------------------------------------------
  dz_par <- list(base = c(1, 0), auc = c(2, 0), e12 = c(1, -2.64))
  for (dz in c("base", "auc", "e12")) {
    ns <- sort(unique(unlist(nplan$ns[nplan$design == dz])))
    for (n in ns) add(bc_cell("1b", sprintf("null_%s_n%d", dz, n), "design", n, 2000L, "y ~ x + d",
                              sprintf("matched logistic null of block 2 (%s)", dz), scenario = "null", design = dz, link = "logit",
                              s = dz_par[[dz]][1], c0 = dz_par[[dz]][2], xdist = "uniform"))
  }
  ns_skew <- sort(unique(unlist(nplan$ns[nplan$design == "base" & nplan$link %in% c("probit", "cauchit")])))
  for (n in ns_skew) add(bc_cell("1b", sprintf("null_skew_n%d", n), "design", n, 2000L, "y ~ x + d",
                                 "matched null of the block 4 skew cells (H5)", scenario = "null", design = "skew", link = "logit",
                                 s = 1, c0 = 0, xdist = "skewed"))
  for (ic in c(-4.9, -1)) for (n in c(100L, 150L, 200L, 300L, 500L, 1000L, 2000L, 5000L)) {
    tag <- if (ic == -4.9) "sparse49" else "moderate1"
    add(bc_cell("1b", sprintf("%s_n%d", tag, n), "sparse", n, 5000L, "y ~ x",
                sprintf("Figure 8 (%s); size on sparse data%s", if (ic == -4.9) "sparse, 1.5%" else "moderate, 29%",
                        if (ic == -4.9 && n %in% c(200L, 300L, 500L)) "; matched null of block 4 sparse cloglog" else ""),
                scenario = sprintf("gen_sparse_link(%g)", ic), link = "logit", intercept = ic, slope = 1))
  }
  for (n in c(500L, 2000L, 10000L)) add(bc_cell("1b", sprintf("runG_null_sparse_n%d", n), "sparse", n, 2000L, "y ~ x",
                                                "size: run G sparse null (eta = -4 + 0.9x)", scenario = "runG null_sparse",
                                                link = "logit", intercept = -4, slope = 0.9))
  add(bc_cell("1b", "runA_null_p20_n1000", "pstar", 1000L, 2000L, "y ~ x1 + d + z1..z17", "size at p* = 20 (run A null)",
              scen = "null", pstar = 20L, G_extra = "20"))
  add(bc_cell("1b", "null_quad_n100", "dgp_null", 100L, 10000L, "y ~ x", "Sec 4 size at n = 100 with G = 5 and 10",
              scenario = "null_quad", family = "quad", G_extra = "5"))
  ## matched nulls that E3 did not list (review F1, F2, F4): logistic truth at pi* of the design, same covariates
  for (n in c(1000L, 2000L)) add(bc_cell("1b", sprintf("null_crossover_n%d", n), "design_null", n, 5000L, "y ~ x + d",
      "matched null of block 4 crossover (Section A family 6)", scenario = "null crossover", scen = "crossover", notes = PRE_LAUNCH))
  for (sc in CENSUS) add(bc_cell("1b", sprintf("null_census_%s_n1000", sc), "design_null", 1000L, 5000L, CENSUS_FORMULA[[sc]],
      sprintf("matched null of census_%s_n1000%s", sc, if (sc == "omit_2cov") " and of block 4 omit_2cov_n1000" else ""),
      scenario = paste("null census", sc), scen = sc, notes = PRE_LAUNCH))
  add(bc_cell("1b", "null_census_omit_2cov_n500", "design_null", 500L, 5000L, "y ~ x", "matched null of block 4 omit_2cov_n500",
              scenario = "null census omit_2cov", scen = "omit_2cov", notes = PRE_LAUNCH))
  add(bc_cell("1b", "null_rough_n10000", "dgp_null", 10000L, 1000L, "y ~ x", "matched null of run G osc4 at n = 10,000",
              scenario = "null_rough", family = "rough", notes = PRE_LAUNCH))

  ## ---- block 2: symmetric tails ---------------------------------------------------------------------------------------------------
  for (i in seq_len(nrow(nplan))) {
    np <- nplan[i, ]
    for (n in np$ns[[1]]) add(bc_cell("2", sprintf("%s_%s_n%d", np$link, np$design, n), "design", n, 2000L, "y ~ x + d",
      if (isTRUE(np$control)) "symmetric-tail headline: undetectable control (n_env > 40,000)" else "symmetric-tail headline (H1, H2; weighting rule A)",
      scenario = np$link, design = np$design, link = np$link, s = np$s, c0 = np$c0, xdist = "uniform",
      null_block = "1b", null_cell = sprintf("null_%s_n%d", np$design, n)))
  }

  ## ---- block 3: Table 4 battery ------------------------------------------------------------------------------------------------------
  for (i in which(!gp$is_gsens)) {
    z <- gp[i, ]
    nm_param <- if (z$family == "link" && z$param %in% names(PLATEAU_NAME)) PLATEAU_NAME[[z$param]] else z$param
    rep_txt <- "Sec 6.2 n ladder; SI"
    if (z$n == 1000L) {
      rep_txt <- "census (SI Table S2)"
      if (z$family %in% c("quad", "binint", "contint")) rep_txt <- paste(rep_txt, "; Fig 2 severity")
      if ((z$family == "quad" && z$param == "0.02") || (z$family == "binint" && z$param == "0.3") ||
          (z$family == "contint" && z$param == "0.5") || z$family == "link") rep_txt <- paste("Table 4 row;", rep_txt)
      if (z$family == "link" && z$param %in% names(PLATEAU_NAME)) rep_txt <- paste(rep_txt, "; risk-plateau group (E0.2)")
      if (z$family == "link" && z$param == "cloglog") rep_txt <- paste(rep_txt, "; Sec 6.6 text")
    }
    nullfam <- z$family
    add(bc_cell("3", sprintf("%s_%s_n%d", z$family, nm_param, z$n), "dgp_alt", z$n, z$B, FAM_FORMULA[[z$family]], rep_txt,
                scenario = paste(z$family, nm_param), link = if (z$family == "link") nm_param else NA,
                family = z$family, param = z$param, seed_root = z$seed_root, cell_id = z$cell_id,
                null_block = "1a", null_cell = sprintf("null_%s_n%d", nullfam, z$n), source = "grid_power_broad.R",
                stored = "sim_power_broad"))
  }
  id3 <- 200L
  for (lk in c("stk_long", "stk_short", "stk_asym", "cauchit", "t4", "loglog")) for (n in c(200L, 500L, 1000L, 2000L, 5000L)) {
    id3 <- id3 + 1L
    add(bc_cell("3", sprintf("%s_n%d", lk, n), "design", n, if (n == 200L) 10000L else 5000L, "y ~ x + d",
                if (grepl("^stk_", lk)) "Table 4 exact Stukel rows (E0.1); SI n ladder" else "Table 4 new link rows; SI n ladder",
                scenario = lk, design = "base", link = lk, s = 1, c0 = 0, xdist = "uniform", cell_id = id3,
                null_block = "1a", null_cell = sprintf("null_link_n%d", n)))
  }
  for (sc in CENSUS) {
    id3 <- id3 + 1L
    add(bc_cell("3", sprintf("census_%s_n1000", sc), "census", 1000L, 5000L, CENSUS_FORMULA[[sc]], "SI Table S2 census row (June _master_add.R)",
                scenario = sc, scen = sc, cell_id = id3, null_block = "1b", null_cell = sprintf("null_census_%s_n1000", sc),
                notes = "matched null added before launch (review F2)"))
  }

  ## ---- block 4: boundary and grouping ----------------------------------------------------------------------------------------------
  for (i in seq_len(nrow(ge))) {
    z <- ge[i, ]
    rough <- z$scenario != "crossover"
    add(bc_cell("4", sprintf("%s_n%d", z$scenario, z$n), if (rough) "dgp_alt" else "crossover", z$n, z$B,
                if (rough) "y ~ x" else "y ~ x + d",
                if (z$n == 1000L) "Sec 6.3 boundary text" else "Fig 3 (plots n = 2000)",
                scenario = z$scenario, family = if (rough) "rough" else NA, param = if (rough) z$shape else NA,
                seed_root = z$seed_root, cell_id = z$cell_id, source = "grid_edge_loses.R", stored = "sim_edge_loses",
                null_block = if (rough) "1a" else "1b",
                null_cell = if (rough) sprintf("null_rough_n%d", z$n) else sprintf("null_crossover_n%d", z$n),
                notes = if (rough) "" else "matched null added before launch (review F1)"))
  }
  ## all 24 grid_proj_power.R cells with their old seeds, paired with the stored projection p-values (B "as now"; E3 listed
  ## only omit_2cov, review F7 adds the other 22). Alternatives are matched to the projection grid's own nulls where it has
  ## one (so the projection column can be size-adjusted too), binint to the block 1a null, omit_2cov to its new design null.
  for (i in seq_len(nrow(gj))) {
    z <- gj[i, ]; sc <- z$scenario; key <- paste(sc, z$param, sep = "|")
    if (sc == "omit_2cov") {
      nm <- sprintf("omit_2cov_n%d", z$n); gen <- "omit_2cov"; fam <- NA; par <- NA; lk <- NA
      nb <- "1b"; nc <- sprintf("null_census_omit_2cov_n%d", z$n); nt <- ""
    } else if (grepl("^null_", sc)) {
      fam <- sub("^null_", "", sc); nm <- sprintf("proj_%s_n%d", sc, z$n); gen <- "dgp_null"; par <- NA; lk <- NA
      nb <- NA; nc <- NA; nt <- PRE_LAUNCH
    } else {
      fam <- sc; par <- z$param; gen <- "dgp_alt"; nt <- PRE_LAUNCH
      pn <- if (sc == "link" && par %in% names(PLATEAU_NAME)) PLATEAU_NAME[[par]] else par
      lk <- if (sc == "link") pn else NA
      nm <- sprintf("proj_%s_%s_n%d", sc, pn, z$n)
      if (sc == "binint") { nb <- "1a"; nc <- sprintf("null_binint_n%d", z$n) } else { nb <- "4"; nc <- sprintf("proj_null_%s_n%d", sc, z$n) }
    }
    add(bc_cell("4", nm, gen, z$n, 500L, PROJ_FORMULA[[sc]], "Sec 6.3 projection head-to-head (paired with stored p_proj)",
                scenario = sc, family = fam, param = par, link = lk, seed_root = z$seed_root, cell_id = z$cell_id,
                source = "grid_proj_power.R", pair_proj = TRUE, stored = "proj_power_grid", stored_key = key,
                null_block = nb, null_cell = nc, notes = nt))
  }
  id4 <- 200L
  for (lk in c("probit", "cauchit")) {
    ns <- nplan$ns[[which(nplan$design == "base" & nplan$link == lk)]]
    for (n in ns) { id4 <- id4 + 1L
      add(bc_cell("4", sprintf("%s_skew_n%d", lk, n), "design", n, 2000L, "y ~ x + d", "H5 grouping cost under a skewed covariate",
                  scenario = lk, design = "skew", link = lk, s = 1, c0 = 0, xdist = "skewed", cell_id = id4,
                  null_block = "1b", null_cell = sprintf("null_skew_n%d", n))) }
  }
  for (n in c(200L, 300L, 500L)) { id4 <- id4 + 1L
    add(bc_cell("4", sprintf("sparse49_cloglog_n%d", n), "sparse", n, 2000L, "y ~ x",
                "events < groups: sparse cloglog, reported unconditionally and split by events < G",
                scenario = "sparse cloglog", link = "cloglog", intercept = -4.9, slope = 1, cell_id = id4,
                null_block = "1b", null_cell = sprintf("sparse49_n%d", n))) }

  ## ---- block 5: rest of the paper -------------------------------------------------------------------------------------------------
  for (ps in c(3L, 5L, 10L, 20L)) for (sc in c("null", "cloglog", "quad", "binint"))
    add(bc_cell("5", sprintf("runA_%s_p%d", sc, ps), "pstar", 1000L, 2000L, sprintf("y ~ x1 + d + %d noise", ps - 3L),
                "tab:pstar, tab:giviti (run A; G = 10 and 20 in one replicate)", scenario = paste("runA", sc), scen = sc,
                pstar = ps, G_extra = "20", null_block = if (sc == "null") NA else "5",
                null_cell = if (sc == "null") NA else sprintf("runA_null_p%d", ps)))
  NSF <- c(500L, 1000L, 2000L, 5000L, 20000L); MS <- c(15L, 25L, 50L, 100L, 200L)
  cf <- do.call(rbind, lapply(NSF, function(n) {                   # run_F_grule.R, then one cell per (n, G)
    G <- pmax(4L, pmin(as.integer(round(n / MS)), 400L))
    unique(data.frame(n = n, m = MS, G = G, stringsAsFactors = FALSE))
  }))
  cf <- do.call(rbind, lapply(split(cf, paste(cf$n, cf$G)), function(z) data.frame(n = z$n[1], G = z$G[1], m = paste(z$m, collapse = "/"))))
  cf <- cf[order(cf$n, -cf$G), ]
  for (sc in c("null", "cloglog", "probit")) for (i in seq_len(nrow(cf))) {
    z <- cf[i, ]
    add(bc_cell("5", sprintf("runF_%s_n%d_G%d", sc, z$n, z$G), "design", z$n, if (z$n >= 20000L) 800L else 2000L, "y ~ x + d",
                sprintf("tab:grule, Fig 9 (run F, m = %s)", z$m), scenario = paste("runF", sc), design = "base",
                link = if (sc == "null") "logit" else sc, s = 1, c0 = 0, xdist = "uniform",
                G_extra = if (z$G == 10L) "" else as.character(z$G), null_block = if (sc == "null") NA else "5",
                null_cell = if (sc == "null") NA else sprintf("runF_null_n%d_G%d", z$n, z$G)))
  }
  for (sc in c("null", "quad", "binint", "osc4", "null_sparse", "cloglog_sparse")) for (n in c(500L, 2000L, 10000L)) {
    sparse <- grepl("sparse", sc)
    add(bc_cell("5", sprintf("runG_%s_n%d", sc, n), if (sparse) "sparse" else "runG", n, if (n >= 10000L) 1000L else 2000L,
                if (sc == "osc4" || sparse) "y ~ x" else "y ~ x + d", "Sec 6.9 hold-out text (run G; both G arms in one replicate)",
                scenario = paste("runG", sc), scen = if (sparse) NA else sc,
                link = if (sc == "cloglog_sparse") "cloglog" else if (sparse) "logit" else NA,
                intercept = if (sparse) -4 else NA, slope = if (sparse) 0.9 else NA,
                ## osc4 is drawn on x alone and fitted y ~ x: matched to the rough null of the same design (review F4)
                null_block = if (sc %in% c("null", "null_sparse")) NA else if (sc == "osc4") (if (n == 10000L) "1b" else "1a") else "5",
                null_cell = if (sc %in% c("null", "null_sparse")) NA else if (sc == "osc4") sprintf("null_rough_n%d", n)
                            else sprintf("runG_%s_n%d", if (sparse) "null_sparse" else "null", n)))
  }
  k2 <- rbind(expand.grid(n = c(500L, 1000L, 2000L), scen = c("null", "cloglog", "loglog"), stringsAsFactors = FALSE),
              data.frame(n = 1000L, scen = c("probit", "binint", "osc4")), data.frame(n = 5000L, scen = c("null", "probit")))
  for (i in seq_len(nrow(k2))) {
    z <- k2[i, ]
    add(bc_cell("5", sprintf("runK2_%s_n%d", z$scen, z$n), "runK", z$n, 2000L, if (z$scen == "osc4") "y ~ x" else "y ~ x + d",
                "tab:ao, Sec 6.8 (run K2, plus EDGE-ao unit)", scenario = paste("runK2", z$scen), scen = z$scen, ao = TRUE,
                null_block = if (z$scen == "null") NA else if (z$scen == "osc4") "1a" else "5",
                null_cell = if (z$scen == "null") NA else if (z$scen == "osc4") sprintf("null_rough_n%d", z$n) else sprintf("runK2_null_n%d", z$n)))
  }
  NSH <- c(500L, 1000L, 2000L, 5000L)
  for (n in NSH) add(bc_cell("5", sprintf("runH_null_n%d", n), "external", n, 2000L, "none (frozen predictions)",
                             "tab:external null (run H)", scenario = "runH null", kind = "null", dist_s = 1, type = "external"))
  hc <- expand.grid(kind = c("temp", "asym_hi", "asym_lo", "probit_comp"), s = c(0.85, 0.7), n = NSH, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(hc))) {
    z <- hc[i, ]
    add(bc_cell("5", sprintf("runH_%s_s%g_n%d", z$kind, z$s, z$n), "external", z$n, 2000L, "none (frozen predictions)",
                "tab:external, Sec 6.10 (run H, plus EDGE-sym external)", scenario = paste("runH", z$kind), kind = z$kind,
                dist_s = z$s, type = "external", null_block = "5", null_cell = sprintf("runH_null_n%d", z$n)))
  }

  ## ---- block 6: real data -------------------------------------------------------------------------------------------------------
  for (kd in c("null", "link", "nonlin", "int"))
    add(bc_cell("6", sprintf("glow_%s_n500", kd), "glow", 500L, 2000L, "y ~ age + weight + priorfrac + momfrac",
                "tab:rdpower, Sec 7.6 (GLOW semi-synthetic)", scenario = paste("glow", kd), scen = kd,
                null_block = if (kd == "null") NA else "6", null_cell = if (kd == "null") NA else "glow_null_n500"))
  real <- list(c("beetle_logit", "y ~ dose", "tab:concord (beetle, logit fit)", "8"),
               c("beetle_cloglog", "y ~ dose, cloglog link", "tab:concord (beetle, cloglog refit; link-agnostic tests only)", "8"),
               c("lbw_additive", "y ~ age + lwt + race + smoke", "tab:concord (LBW)", ""),
               c("lbw_fix", "y ~ age + lwt + race + smoke + age:lwt + smoke:lwt", "tab:concord (LBW, with interactions)", ""),
               c("vaso_linear", "y ~ vol + rate", "tab:concord (vaso)", "6"),
               c("vaso_log", "y ~ log(vol) + log(rate)", "tab:concord (vaso, log fit)", "6"),
               c("uis_model19", "Hosmer model 19", "Sec 7.4 UIS text", ""),
               c("diabetes_dev", "13 main effects", "tab:diabetes development half (G sweep), Fig 10", "50,200,1000,2000"))
  ## E8.5: every one-off fit on its stored row order (rep 0) and on random row orders (rep k, set.seed(20260914 + k)),
  ## 1000 orders, 100 for each Diabetes-130 half; B counts the rows of the output (1 + orders)
  RO <- "stored row order plus random row orders, seed 20260914 + k (E8.5)"
  for (r in real) {
    no <- if (grepl("^diabetes", r[1])) 100L else 1000L
    add(bc_cell("6", paste0("real_", r[1]), "real", NA_integer_, 1L + no, r[2], r[3], scenario = r[1], dataset = r[1],
                G_extra = r[4], orders = no, notes = paste0(if (r[1] == "beetle_cloglog") paste("one-off fit;", PRE_LAUNCH) else "one-off fit",
                                                            "; ", RO)))
  }
  add(bc_cell("6", "real_diabetes_val", "real", NA_integer_, 101L, "frozen development-half model", "tab:diabetes validation half (G sweep)",
              scenario = "diabetes_val", dataset = "diabetes_val", type = "external", G_extra = "10,50,200,1000,2000", orders = 100L,
              notes = paste0("one-off, external mode; G sweep added before launch (review F8); ", RO)))

  ## ---- block 7: long tail ---------------------------------------------------------------------------------------------------------
  have <- do.call(rbind, R)
  probit_base_n <- unique(c(have$n[have$block == "3" & have$family %in% "link" & have$param %in% "probit"],
                            have$n[have$block == "2" & have$link %in% "probit" & have$design %in% "base"]))
  for (n in setdiff(c(1000L, 5000L, 20000L, 50000L), probit_base_n)) {
    nul_1b <- sprintf("null_base_n%d", n)
    has_null <- any(have$block == "1b" & have$cell == nul_1b)
    add(bc_cell("7", sprintf("runD_probit_n%d", n), "design", n, if (n >= 20000L) 800L else 2000L, "y ~ x + d",
                "Sec 8 probit ladder (run D1)", scenario = "runD probit", design = "base", link = "probit", s = 1, c0 = 0,
                xdist = "uniform", null_block = if (has_null) "1b" else "7",
                null_cell = if (has_null) nul_1b else sprintf("runD_null_n%d", n)))
    if (!has_null) add(bc_cell("7", sprintf("runD_null_n%d", n), "design", n, if (n >= 20000L) 800L else 2000L, "y ~ x + d",
                               "matched null of the run D probit ladder", scenario = "runD null", design = "base", link = "logit",
                               s = 1, c0 = 0, xdist = "uniform"))
  }
  for (sc in c("null", "cloglog"))
    add(bc_cell("7", sprintf("runJ_%s_n2000", sc), "design", 2000L, 2000L, "y ~ x + d",
                "Sec 6.9 weights and tiny-m text (run J, EDGE only)", scenario = paste("runJ", sc), design = "base",
                link = if (sc == "null") "logit" else "cloglog", s = 1, c0 = 0, xdist = "uniform",
                G_extra = "20,50,100,200,400,1000", type = "edgeonly",
                null_block = if (sc == "null") NA else "7", null_cell = if (sc == "null") NA else "runJ_null_n2000"))

  C <- do.call(rbind, R)
  ## new cell ids per block, then seed bases
  for (b in BC_BLOCKS) {
    i <- which(C$block == b & is.na(C$seed_root) & C$generator != "real")
    j <- i[is.na(C$cell_id[i])]
    if (length(j)) C$cell_id[j] <- seq_along(j)
  }
  C$seed_family <- ifelse(is.na(C$seed_root), ifelse(C$generator == "real", "none", "new"), paste0("old:", C$seed_root))
  C$seed_base <- ifelse(C$seed_family == "none", NA_real_,
                        ifelse(is.na(C$seed_root), 1e8 + BC_BLOCKNUM[C$block] * 1e7 + C$cell_id * 1e4, C$seed_root + C$cell_id * 1e5))
  C$G_list <- vapply(seq_len(nrow(C)), function(i) {
    ex <- if (nzchar(C$G_extra[i])) strsplit(C$G_extra[i], ",")[[1]] else character(0)
    if (C$type[i] == "external")
      return(paste(unique(c(if (is.na(C$n[i])) "max(6, n/25)" else as.character(max(6L, as.integer(round(C$n[i] / 25)))), ex)), collapse = ","))
    rule <- if (is.na(C$n[i])) "rule" else as.character(max(10, round(C$n[i] / 25)))
    paste(unique(c("10", rule, ex)), collapse = ",")
  }, character(1))
  ## role of each cell in the summary (review F6: set from the generator, not from the cell name)
  C$role <- ifelse(C$generator == "real", "real data",
            ifelse(C$generator %in% c("dgp_null", "design_null") | (C$generator %in% c("design", "sparse") & C$link %in% "logit") |
                   (C$generator %in% c("pstar", "runG", "runK", "glow") & C$scen %in% "null") |
                   (C$type == "external" & C$kind %in% "null"), "null", "alternative"))
  stopifnot(!anyDuplicated(paste(C$block, C$cell)))
  rownames(C) <- NULL
  C
}

## ---- checks ----------------------------------------------------------------------------------------------------------------------
## every pair of result cells (block 0 excluded: it re-draws cells of other blocks on purpose) whose seed ranges intersect
seed_overlaps <- function(C) {
  z <- C[!is.na(C$seed_base) & C$block != "0", ]
  lo <- z$seed_base + 1; hi <- z$seed_base + z$B
  o <- order(lo); z <- z[o, ]; lo <- lo[o]; hi <- hi[o]
  out <- list(); k <- length(lo)
  for (i in seq_len(k)) {
    j <- i + 1L
    while (j <= k && lo[j] <= hi[i]) {
      both_old <- z$seed_family[i] != "new" && z$seed_family[j] != "new"
      out[[length(out) + 1]] <- data.frame(cell_a = paste(z$block[i], z$cell[i]), cell_b = paste(z$block[j], z$cell[j]),
        family_a = z$seed_family[i], family_b = z$seed_family[j], offset = lo[j] - lo[i], shared_seeds = min(hi[i], hi[j]) - lo[j] + 1,
        kind = if (both_old) "inherited from the July grids" else "conflict", stringsAsFactors = FALSE)
      j <- j + 1L
    }
  }
  if (length(out)) do.call(rbind, out) else data.frame(cell_a = character(0), cell_b = character(0), family_a = character(0),
    family_b = character(0), offset = numeric(0), shared_seeds = numeric(0), kind = character(0))
}

## rebuilt old ids against the seeds written by the July runs: stored seed - rep must be the table's seed_base for every
## replicate, and the stored number of replicates must be the table's B
verify_old_ids <- function(C, simdir = BC_SIMDIR) {
  suppressPackageStartupMessages(library(data.table))
  base_of <- function(dt, keys) dt[, .(stored_base = unique(seed - rep)[1], n_bases = uniqueN(seed - rep), stored_B = max(rep)), by = keys]
  res <- list()
  a <- fread(file.path(simdir, "sim_null_pvalues.csv"), select = c("family", "n", "G", "rep", "seed"))
  a <- base_of(a, c("family", "n", "G"))
  z <- C[C$source == "grid_null.R" & C$block %in% c("1a", "0"), ]
  for (i in seq_len(nrow(z))) {
    s <- a[family == z$family[i] & n == z$n[i] & G == 10L]
    res[[length(res) + 1]] <- data.frame(block = z$block[i], cell = z$cell[i], file = "sim_null_pvalues.csv", table_base = z$seed_base[i],
      stored_base = s$stored_base, n_bases = s$n_bases, table_B = if (z$block[i] == "0") 10000L else z$B[i], stored_B = s$stored_B)
  }
  a <- fread(file.path(simdir, "sim_power_broad_pvalues.csv"), select = c("family", "param", "n", "G", "rep", "seed"),
             colClasses = list(character = "param"))
  a <- base_of(a, c("family", "param", "n", "G"))
  z <- C[C$source == "grid_power_broad.R", ]
  for (i in seq_len(nrow(z))) {
    s <- a[family == z$family[i] & param == z$param[i] & n == z$n[i] & G == 10L]
    res[[length(res) + 1]] <- data.frame(block = z$block[i], cell = z$cell[i], file = "sim_power_broad_pvalues.csv", table_base = z$seed_base[i],
      stored_base = s$stored_base, n_bases = s$n_bases, table_B = if (z$block[i] == "0") ifelse(z$n[i] == 200, 10000L, 5000L) else z$B[i],
      stored_B = s$stored_B)
  }
  a <- fread(file.path(simdir, "sim_edge_loses_pvalues.csv"), select = c("scenario", "n", "G", "rep", "seed"))
  a <- base_of(a, c("scenario", "n", "G"))
  z <- C[C$source == "grid_edge_loses.R", ]
  for (i in seq_len(nrow(z))) {
    s <- a[scenario == z$scenario[i] & n == z$n[i]]
    if (z$block[i] == "0") s <- a[scenario == "crossover" & n == z$n[i]]
    res[[length(res) + 1]] <- data.frame(block = z$block[i], cell = z$cell[i], file = "sim_edge_loses_pvalues.csv", table_base = z$seed_base[i],
      stored_base = s$stored_base, n_bases = s$n_bases, table_B = if (z$block[i] == "0") 5000L else z$B[i], stored_B = s$stored_B)
  }
  a <- fread(file.path(simdir, "proj_power_grid_pvalues.csv"), select = c("scenario", "param", "n", "rep", "seed"),
             colClasses = list(character = "param"))
  a <- base_of(a, c("scenario", "param", "n"))
  z <- C[C$source == "grid_proj_power.R", ]
  for (i in seq_len(nrow(z))) {
    sk <- strsplit(z$stored_key[i], "|", fixed = TRUE)[[1]]
    s <- a[scenario == sk[1] & param == sk[2] & n == z$n[i]]
    res[[length(res) + 1]] <- data.frame(block = z$block[i], cell = z$cell[i], file = "proj_power_grid_pvalues.csv", table_base = z$seed_base[i],
      stored_base = s$stored_base, n_bases = s$n_bases, table_B = z$B[i], stored_B = s$stored_B)
  }
  V <- do.call(rbind, res)
  V$ok <- !is.na(V$stored_base) & V$table_base == V$stored_base & V$n_bases == 1 & V$table_B == V$stored_B
  V
}

## counts per block and the E3 rows each block must contain
cells_expected <- function(C, nplan = bc_read_nplan()) {
  b <- function(x) C[C$block == x, ]
  n2 <- sum(lengths(nplan$ns))
  chk <- list(
    c("0: three old-seed identity cells", sum(!is.na(b("0")$stored)) == 3),
    c("0: worker timing cell at n = 5000", any(b("0")$n == 5000L & grepl("worker", b("0")$reproduces))),
    c("1a: 25 grid_null cells, ids 1-25, B 10,000", nrow(b("1a")) == 25 && setequal(b("1a")$cell_id, 1:25) && all(b("1a")$B == 10000L)),
    c("1a: quad and link n = 1000 carry G 6, 8, 12, 14, 20", sum(b("1a")$G_extra == "6,8,12,14,20") == 2),
    c("1b: logistic null for every (design, n) of block 2",
      all(vapply(seq_len(nrow(b("2"))), function(i) b("2")$null_cell[i] %in% b("1b")$cell, logical(1)))),
    c("1b: skew null at every n of the block 4 skew cells",
      all(b("4")$null_cell[b("4")$design %in% "skew"] %in% b("1b")$cell)),
    c("1b: gen_sparse_link(-4.9) and (-1) at 8 n each, B 5000", sum(b("1b")$B == 5000L & b("1b")$generator == "sparse") == 16),
    c("1b: run G sparse null at 500, 2000, 10,000", sum(grepl("^runG_null_sparse", b("1b")$cell)) == 3),
    c("1b: run A null p* = 20", any(b("1b")$cell == "runA_null_p20_n1000")),
    c("1b: quad null n = 100, G 5 and 10, B 10,000", any(b("1b")$cell == "null_quad_n100" & b("1b")$B == 10000L & b("1b")$G_extra == "5")),
    c(sprintf("2: {probit, cauchit, t4, loglog} x {base, auc, e12} x {n_lo, n_hi} = %d cells from nplan.csv", n2), nrow(b("2")) == n2),
    c("3: all 117 grid_power_broad cells, ids 1-117", sum(b("3")$source == "grid_power_broad.R") == 117 &&
        setequal(b("3")$cell_id[b("3")$source == "grid_power_broad.R"], 1:117)),
    c("3: exact Stukel long, short, asym x 5 n", sum(grepl("^stk_", b("3")$cell)) == 15),
    c("3: cauchit, t4, loglog x 5 n", sum(grepl("^(cauchit|t4|loglog)_n", b("3")$cell)) == 15),
    c("3: 8 census scenarios at n = 1000", sum(grepl("^census_", b("3")$cell)) == 8),
    c("3: plateau generators keep old seeds (15 cells)", sum(grepl("plateau", b("3")$cell) & !is.na(b("3")$seed_root)) == 15),
    c("4: 6 grid_edge_loses cells with old seeds", sum(b("4")$source == "grid_edge_loses.R") == 6),
    c("4: all 24 grid_proj_power cells with their seeds (ids 1-24), B 500, paired with the stored projection p-values",
      sum(b("4")$pair_proj) == 24 && setequal(b("4")$cell_id[b("4")$pair_proj], 1:24) && all(b("4")$B[b("4")$pair_proj] == 500L) &&
        !anyDuplicated(paste(b("4")$stored_key, b("4")$n)[b("4")$pair_proj])),
    c("1b: crossover null at n 1000, 2000 and one null per census design (+ omit_2cov n 500), B 5000",
      all(c("null_crossover_n1000", "null_crossover_n2000", sprintf("null_census_%s_n1000", CENSUS), "null_census_omit_2cov_n500") %in% b("1b")$cell) &&
        all(b("1b")$B[b("1b")$generator == "design_null"] == 5000L)),
    c("1b: rough null at n = 10,000 for run G osc4", any(b("1b")$cell == "null_rough_n10000")),
    c("every alternative outside block 0 has a matched null that exists at the same n (C2)", {
      alt <- C[C$role == "alternative" & C$block != "0", ]
      key <- paste(C$block, C$cell)
      all(!is.na(alt$null_cell)) && all(paste(alt$null_block, alt$null_cell) %in% key) &&
        all(C$role[match(paste(alt$null_block, alt$null_cell), key)] == "null") &&
        all(C$n[match(paste(alt$null_block, alt$null_cell), key)] == alt$n) }),
    c("6: Diabetes G sweep 10, 50, 200, 1000, 2000 on both halves; beetle cloglog refit",
      all(c("10", "50", "200", "1000", "2000") %in% strsplit(C$G_list[C$cell == "real_diabetes_dev"], ",")[[1]]) &&
        all(c("10", "50", "200", "1000", "2000") %in% strsplit(C$G_list[C$cell == "real_diabetes_val"], ",")[[1]]) &&
        any(b("6")$cell == "real_beetle_cloglog")),
    c("4: skew x {probit, cauchit} at the block 2 base n", sum(b("4")$design %in% "skew") ==
        sum(lengths(nplan$ns[nplan$design == "base" & nplan$link %in% c("probit", "cauchit")]))),
    c("4: sparse cloglog at 200, 300, 500", sum(grepl("^sparse49_cloglog", b("4")$cell)) == 3),
    c("5: run A 4 scenarios x 4 p* (32 original cells, G 10 and 20 inside)", sum(grepl("^runA_", b("5")$cell)) == 16),
    c("5: run F 75 cells with n = 20,000 de-duplicated to 9 (= 69)", sum(grepl("^runF_", b("5")$cell)) == 69 &&
        sum(grepl("^runF_.*_n20000_", b("5")$cell)) == 9),
    c("5: run G 6 scenarios x 3 n", sum(grepl("^runG_", b("5")$cell)) == 18),
    c("5: run K2 14 cells with EDGE-ao", sum(grepl("^runK2_", b("5")$cell) & b("5")$ao) == 14),
    c("5: run H 36 cells, external", sum(grepl("^runH_", b("5")$cell) & b("5")$type == "external") == 36),
    c("6: GLOW 4 cells n = 500", sum(grepl("^glow_", b("6")$cell)) == 4),
    c("6: real one-off fits on the stored row order and 1000 random orders, 100 for each Diabetes-130 half (E8.5)", {
      rl <- b("6")[b("6")$generator == "real", ]
      nrow(rl) == 9 && all(rl$orders == ifelse(grepl("^real_diabetes", rl$cell), 100L, 1000L)) && all(rl$B == 1L + rl$orders) }),
    c("6: beetle, LBW, vaso, UIS, Diabetes-130 one-off fits", all(c("real_beetle_logit", "real_lbw_additive", "real_vaso_linear",
        "real_uis_model19", "real_diabetes_dev", "real_diabetes_val") %in% b("6")$cell)),
    c("7: run D probit ladder to n = 50,000 without cells already in blocks 2-4", any(b("7")$cell == "runD_probit_n50000") &&
        !any(b("7")$cell %in% c("runD_probit_n1000", "runD_probit_n5000"))),
    c("7: run J, EDGE only", sum(b("7")$type == "edgeonly") == 2))
  data.frame(check = vapply(chk, `[`, "", 1), ok = as.logical(vapply(chk, `[`, "", 2)), stringsAsFactors = FALSE)
}

## ---- the E4 run order ------------------------------------------------------------------------------------------------------------
## "0; 2 with its 1b nulls; 1a and the rest of 1b; 3 (the n = 1000 slice first); 4; 5; 6; 7", then a summary of every block.
## One row per Rscript call of run_M_battery.R; an empty 'cells' means the whole block (finished cells are skipped).
launch_plan <- function(C) {
  b2null <- unique(C$null_cell[C$block == "2" & C$null_block %in% "1b"])
  n1000 <- C$cell[C$block == "3" & C$n %in% 1000L]
  S <- rbind(c("0", "", "block 0: identity, harness against package, worker timing (E0.8); writes RULE_weighting.txt"),
             c("1b", paste(b2null, collapse = ","), "the matched nulls of block 2"),
             c("2", "", "block 2, symmetric tails"),
             c("1a", "", "block 1a, size grid"),
             c("1b", "", "the rest of block 1b"),
             c("3", paste(n1000, collapse = ","), "block 3, the n = 1000 slice"),
             c("3", "", "the rest of block 3"),
             c("4", "", "block 4"), c("5", "", "block 5"), c("6", "", "block 6"), c("7", "", "block 7"))
  P <- data.frame(step = seq_len(nrow(S)), action = "run", block = S[, 1], cells = S[, 2], note = S[, 3], stringsAsFactors = FALSE)
  rbind(P, data.frame(step = nrow(S) + seq_along(BC_BLOCKS[-1]), action = "summary", block = BC_BLOCKS[-1], cells = "",
                      note = "summary of every finished cell of the block", stringsAsFactors = FALSE))
}
