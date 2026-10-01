## pkg280_reference.R -- outputs of ebrahim.gof that the 2.8.0 change (joint Stukel, basis "sym",
## weights = "score", G = "auto", few-events warning) must NOT move.
##
## Run before the change:  Rscript pkg280_reference.R before devtree
## Run after the change:   Rscript pkg280_reference.R after installed     (or: after devtree)
##
## "devtree" loads the package with pkgload::load_all() from the dev tree; "installed" uses library().
## The "after" run also compares against pkg280_reference_before.rds and prints every difference.
## Expected differences after the change: the Stukel row of the batteries (joint score instead of the
## marginal sum) and one new row, DEF.sym. Everything else must be identical.
##
## Run after the 2.8.0 review fixes:  Rscript pkg280_reference.R after2 devtree
## The "after2" run compares against pkg280_reference_after.rds, element by element, and lists every
## number that differs by more than 1e-10 (relative to max(1, |value|)) and every other change.
##
## Run after the column-rule correction (2.7.0 filter kept, columns scaled):  Rscript pkg280_reference.R after3 devtree
## The same comparison as after2, against pkg280_reference_after.rds; the outputs go to pkg280_reference_after3.rds.
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

args   <- commandArgs(trailingOnly = TRUE)
mode   <- if (length(args) >= 1) args[1] else "before"
loader <- if (length(args) >= 2) args[2] else "devtree"
stopifnot(mode %in% c("before", "after", "after2", "after3"), loader %in% c("devtree", "installed"))

PKG <- Sys.getenv("EBRAHIM_GOF_SRC")
OUT <- edge_path("declarations")
options(width = 200, warn = 1)

if (loader == "devtree") {
  library(ebrahim.gof)  # archive: was load_all() of the author's development tree; install ebrahim.gof (>= 2.9.0) from CRAN instead
} else {
  library(ebrahim.gof)
}
ver  <- as.character(utils::packageVersion("ebrahim.gof"))
head <- tryCatch(system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE),
                 error = function(e) NA_character_)
cat("ebrahim.gof", ver, "loaded by", loader, "| dev tree HEAD", head, "\n")

tm <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  v  <- withCallingHandlers(expr, warning = function(w) {
    cat("  warning in", label, ":", conditionMessage(w), "\n"); invokeRestart("muffleWarning")
  })
  cat(sprintf("  %-40s %.1f s\n", label, proc.time()[["elapsed"]] - t0))
  v
}

## ---- fixtures -------------------------------------------------------------------------------------------
make_fit <- function(seed = 1, n = 500, link = "logit") {         # test-run-all-gof / test-def-gof style
  set.seed(seed)
  x <- runif(n, -3, 3)
  eta <- 0.6 * x
  p <- if (link == "cloglog") 1 - exp(-exp(eta)) else 1 / (1 + exp(-eta))
  glm(rbinom(n, 1, p) ~ x, family = binomial())
}
gof_demo <- ebrahim.gof::gof_demo
wrong <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial())
right <- glm(outcome ~ poly(age, 2) + bmi + sex + treatment, data = gof_demo, family = binomial())

set.seed(42)
x1 <- rnorm(400); x2 <- runif(400, -2, 2)
fit42 <- glm(rbinom(400, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2)) ~ x1 + x2, family = binomial())

set.seed(5)
x1 <- rnorm(400); d <- factor(sample(c("A", "B"), 400, replace = TRUE))
fitcat <- glm(rbinom(400, 1, plogis(0.3 + 0.6 * x1 + ifelse(d == "B", 0.5, 0))) ~ x1 + d,
              family = binomial())

fits <- list(mf1_500 = make_fit(1, 500), mf1_600 = make_fit(1, 600), mf2_600 = make_fit(2, 600),
             mf3_600 = make_fit(3, 600), cloglog7 = make_fit(7, 1500, "cloglog"),
             seed42 = fit42, cat5 = fitcat, demo_wrong = wrong, demo_right = right)

ref <- list(meta = list(version = ver, loader = loader, head = head, time = Sys.time(),
                        R = R.version.string))

## ---- run.all.gof fast battery ------------------------------------------------------------------------------
ref$battery <- lapply(names(fits), function(nm)
  tm(paste("battery", nm), as.data.frame(run.all.gof(fits[[nm]], include_slow = FALSE, install = "no"))))
names(ref$battery) <- names(fits)
ref$battery_pred_only <- tm("battery (y, ph) mf2_600",
  as.data.frame(suppressMessages(run.all.gof(as.numeric(fits$mf2_600$y), fitted(fits$mf2_600),
                                             include_slow = FALSE, install = "no"))))

## ---- def.gof unit form, edge.gof, def.ensemble.gof ------------------------------------------------------------
ref$def <- list()
for (nm in names(fits)) for (b in c("poly3", "poly2", "stukel")) for (m in c("satterthwaite", "imhof")) {
  key <- paste(nm, b, m, sep = "|")
  ref$def[[key]] <- tryCatch(def.gof(fits[[nm]], basis = b, method = m), error = function(e) conditionMessage(e))
}
for (G in c(5, 20)) ref$def[[paste("demo_wrong|poly3|G", G)]] <- def.gof(wrong, G = G)
ref$def[["mf2_600|poly2|yphX"]] <- def.gof(as.numeric(fits$mf2_600$y), fitted(fits$mf2_600),
                                           X = model.matrix(fits$mf2_600), basis = "poly2")
ref$def[["mf2_600|poly2|yph_naive"]] <- suppressWarnings(def.gof(as.numeric(fits$mf2_600$y),
                                                                 fitted(fits$mf2_600), basis = "poly2"))
ref$edge <- list(poly3 = edge.gof(fits$mf1_500), stukel = edge.gof(fits$mf1_500, basis = "stukel"),
                 ensemble = edge.gof(fits$mf1_500, basis = "ensemble"))

ref$ensemble <- list()
for (nm in names(fits)) {
  ref$ensemble[[paste(nm, "cct")]]    <- def.ensemble.gof(fits[[nm]])
  ref$ensemble[[paste(nm, "cct+EF")]] <- def.ensemble.gof(fits[[nm]], add_ef = TRUE)
  ref$ensemble[[paste(nm, "minp")]]   <- def.ensemble.gof(fits[[nm]], combine = "minp")
  ref$ensemble[[paste(nm, "fisher")]] <- def.ensemble.gof(fits[[nm]], combine = "fisher")
}
ref$ensemble[["mf1_600 extra"]] <- def.ensemble.gof(fits$mf1_600, extra_pvalues = c(Tsiatis = 0.2))
ref$ensemble[["mf1_600 edges"]] <- edges.gof(fits$mf1_600)
ref$ensemble[["mf1_600 basis=ensemble"]] <- def.gof(fits$mf1_600, basis = "ensemble")

## ---- gof.features, cdef.gof ---------------------------------------------------------------------------------
ref$features <- tm("gof.features gof_demo wrong", suppressWarnings(gof.features(wrong)))
ref$features_right <- tm("gof.features gof_demo right", suppressWarnings(gof.features(right)))
ref$cdef <- list(poly = suppressWarnings(cdef.gof(wrong, basis = "poly")),
                 spline = suppressWarnings(cdef.gof(wrong, basis = "spline")))

## ---- calm.gof, shrink.gof (roxygen examples) -------------------------------------------------------------
set.seed(1)
Xc <- matrix(rnorm(300 * 20), 300)
yc <- rbinom(300, 1, 1 / (1 + exp(-(Xc[, 1] - 0.5 * Xc[, 2]))))
ref$calm <- tm("calm.gof example", calm.gof(Xc, yc, lambda = 100))
set.seed(1)
Xs <- matrix(rnorm(400 * 5), 400, 5)
ys <- rbinom(400, 1, plogis(0.3 + Xs %*% c(0.8, -0.5, 0.3, 0, 0)))
ref$shrink <- tm("shrink.gof example", shrink.gof(Xs, ys, lambda = 40, basis = "decile", B = 99, seed = 1))

## ---- legoft, legoft.localize (roxygen examples) -------------------------------------------------------------
set.seed(1)
x1 <- runif(300, -3, 3); x2 <- rnorm(300)
y  <- rbinom(300, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.25 * x1^2))
lfit <- glm(y ~ x1 + x2, family = binomial())
ref$legoft <- tm("legoft example (B = 99)", unclass(legoft(lfit, B = 99, seed = 1)))
set.seed(1)
x1 <- runif(400, -3, 3); x2 <- rnorm(400)
y  <- rbinom(400, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.4 * x1 * x2))
lfit2 <- glm(y ~ x1 + x2, family = binomial())
ref$legoft_localize <- tm("legoft.localize example (B = 99)", unclass(legoft.localize(lfit2, B = 99, seed = 1)))

## ---- save, and compare when this is the "after" run ------------------------------------------------------------
file <- file.path(OUT, sprintf("pkg280_reference_%s.rds", mode))
saveRDS(ref, file)
cat("saved", file, "\n")

if (mode == "after") {
  bf <- file.path(OUT, "pkg280_reference_before.rds")
  if (!file.exists(bf)) stop("no before file to compare against")
  before <- readRDS(bf)
  cat("\n==== comparison against", bf, "(", before$meta$version, "->", ver, ") ====\n")
  same <- function(a, b) identical(a, b)
  for (sec in setdiff(names(before), "meta")) {
    a <- before[[sec]]; b <- ref[[sec]]
    if (sec == "battery") {
      for (nm in names(a)) {
        A <- a[[nm]]; B <- b[[nm]]
        extra <- setdiff(B$Test, A$Test); gone <- setdiff(A$Test, B$Test)
        common <- intersect(A$Test, B$Test)
        diffs <- common[!vapply(common, function(t) same(A[A$Test == t, , drop = FALSE],
                                                         B[B$Test == t, , drop = FALSE]) ||
                                  isTRUE(all.equal(`rownames<-`(A[A$Test == t, ], NULL),
                                                   `rownames<-`(B[B$Test == t, ], NULL), tolerance = 0)),
                                logical(1))]
        cat(sprintf("battery %-10s rows %d -> %d | new: %s | removed: %s | changed: %s\n", nm, nrow(A), nrow(B),
                    paste(extra, collapse = ","), paste(gone, collapse = ","), paste(diffs, collapse = ",")))
        for (t in diffs) {
          cat("   before:"); print(A[A$Test == t, ], row.names = FALSE)
          cat("   after: "); print(B[B$Test == t, ], row.names = FALSE)
        }
      }
    } else if (is.list(a) && !is.data.frame(a) && !is.null(names(a)) && sec %in% c("def", "edge", "ensemble", "cdef")) {
      bad <- names(a)[!vapply(names(a), function(k) same(a[[k]], b[[k]]), logical(1))]
      cat(sprintf("%-16s %d objects, identical: %s%s\n", sec, length(a), length(bad) == 0,
                  if (length(bad)) paste(" | differ:", paste(bad, collapse = ", ")) else ""))
    } else {
      ok <- same(a, b)
      cat(sprintf("%-16s identical: %s%s\n", sec, ok,
                  if (!ok) paste(" | all.equal:", paste(all.equal(a, b), collapse = "; ")) else ""))
    }
  }
}

if (mode %in% c("after2", "after3")) {
  af <- file.path(OUT, "pkg280_reference_after.rds")
  if (!file.exists(af)) stop("no after file to compare against")
  after <- readRDS(af)
  cat("\n==== comparison against", af, "(", after$meta$version, "HEAD", after$meta$head, "-> dev tree HEAD", head,
      ") ====\n")
  TOL <- 1e-10
  diffs <- character(0); maxd <- 0
  walk <- function(a, b, path) {
    if (is.list(a) && is.list(b)) {                                # lists and data frames
      if (!identical(names(a), names(b)) || length(a) != length(b)) {
        diffs <<- c(diffs, sprintf("%s: names or length differ (%s -> %s)", path,
                                   paste(names(a), collapse = ","), paste(names(b), collapse = ",")))
        return(invisible())
      }
      for (i in seq_along(a)) walk(a[[i]], b[[i]], paste0(path, "$", if (is.null(names(a))) i else names(a)[i]))
    } else if (is.numeric(a) && is.numeric(b)) {
      if (length(a) != length(b) || !identical(dim(a), dim(b))) {
        diffs <<- c(diffs, sprintf("%s: length or dim differ", path)); return(invisible())
      }
      if (!identical(is.na(a), is.na(b))) {
        diffs <<- c(diffs, sprintf("%s: NA pattern differs", path)); return(invisible())
      }
      ok <- !is.na(a)
      d  <- ifelse(a[ok] == b[ok], 0, abs(a[ok] - b[ok]))
      if (length(d)) {
        maxd <<- max(maxd, d)
        bad <- d > TOL * pmax(1, abs(a[ok]))
        if (any(bad)) diffs <<- c(diffs, sprintf("%s: %d value(s) differ, max abs diff %.3g (%s -> %s)", path,
                                               sum(bad), max(d[bad]), format(a[ok][bad][1], digits = 10),
                                               format(b[ok][bad][1], digits = 10)))
      }
    } else if (!identical(a, b)) {
      diffs <<- c(diffs, sprintf("%s: %s -> %s", path, paste(format(a), collapse = " "),
                                 paste(format(b), collapse = " ")))
    }
  }
  for (sec in setdiff(names(after), "meta")) {
    n0 <- length(diffs); m0 <- maxd; maxd <- 0
    walk(after[[sec]], ref[[sec]], sec)
    cat(sprintf("%-18s differences beyond %g: %d | max abs numeric diff %.3g\n", sec, TOL, length(diffs) - n0, maxd))
    maxd <- max(m0, maxd)
  }
  if (length(diffs)) { cat("\nDifferences:\n"); cat(paste(" ", diffs), sep = "\n") } else cat("\nNo differences.\n")
  cat("overall max abs numeric diff:", format(maxd, digits = 3), "\n")
}
