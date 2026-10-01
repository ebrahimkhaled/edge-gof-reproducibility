## run_M_blockR_realdata.R -- block R: corrupted predictors on the real cohort, semi-synthetic outcome.
## Contract: paper_EDGE/theory/PREDECLARATION_blockR_realdata_corruption.md (sha256 e80aaef5...), frozen before
## any scenario ran.
##
##   Rscript run_M_blockR_realdata.R [--workers 20] [--reps 500] [--out battery/R]
##
## The cohort, the exclusion, the split seed and the fitted model are run_I_bigdata.R's, reproduced here so the
## two sections describe the same data. The outcome is drawn from the frozen predictions, which makes the frozen
## model true for that replicate; k records then have num_lab_procedures multiplied by four and their prediction
## recomputed at the corrupted value, while the outcome stays as drawn. Every rejection is a false alarm.
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
suppressPackageStartupMessages({ library(parallel); library(data.table) })
## which covariate carries the error. Block R corrupted num_lab_procedures, whose fitted coefficient is
## -0.0004, so a fourfold error moved the linear predictor by 0.08 and nothing happened to any test.
## Block R2 (declaration d27b1d6b) repeats it on number_inpatient, the model's largest coefficient.
RR_VAR <- Sys.getenv("RR_VAR", "num_lab_procedures")
## Paper 3 (RR_COHORT=p3): the cohort of _cohort_p3.R (one admission per patient, no death or hospice discharge,
## split by patient) and the paper's default basis, a cubic in the group's mean risk plus the external intercept.
## Unset, the script reproduces the paper-2 blocks R and R2 exactly as deposited.
RR_COHORT <- Sys.getenv("RR_COHORT", "")

rr_prepare <- function() {
  if (identical(RR_COHORT, "p3")) {
    source(file.path(SIMDIR, "_cohort_p3.R"))
    CO <- cohort_p3(SIMDIR)
    fit <- stats::glm(CO$f, data = CO$Ddev, family = stats::binomial())
    Dval <- CO$Dval
    eta0 <- as.numeric(stats::predict(fit, newdata = Dval, type = "link"))
    b_lab <- unname(stats::coef(fit)[[RR_VAR]])
    chk <- sample.int(nrow(Dval), 100)
    Dc <- Dval[chk, , drop = FALSE]; Dc[[RR_VAR]] <- 4 * Dc[[RR_VAR]]
    ref <- as.numeric(stats::predict(fit, newdata = Dc, type = "response"))
    got <- stats::plogis(eta0[chk] + 3 * b_lab * Dval[[RR_VAR]][chk])
    if (max(abs(ref - got)) > 1e-10) stop("the linear-predictor shortcut disagrees with predict()")
    return(list(fit = fit, Dval = Dval, eta0 = eta0, b_lab = b_lab, p0 = stats::plogis(eta0)))
  }
  D <- read.csv(file.path(SIMDIR, "..", "data_large", "diabetic_data.csv"),
                stringsAsFactors = FALSE, na.strings = c("?", ""))
  D <- D[D$gender %in% c("Female", "Male"), ]
  D$y <- as.integer(D$readmitted == "<30")
  age_mid <- c("[0-10)" = 5, "[10-20)" = 15, "[20-30)" = 25, "[30-40)" = 35, "[40-50)" = 45, "[50-60)" = 55,
               "[60-70)" = 65, "[70-80)" = 75, "[80-90)" = 85, "[90-100)" = 95)
  D$age_n <- unname(age_mid[D$age])
  D$female <- as.integer(D$gender == "Female")
  D$insulin <- factor(D$insulin, levels = c("No", "Steady", "Up", "Down"))
  D$change <- as.integer(D$change == "Ch"); D$diabetesMed <- as.integer(D$diabetesMed == "Yes")
  f <- y ~ age_n + female + time_in_hospital + num_lab_procedures + num_procedures + num_medications +
           number_outpatient + number_emergency + number_inpatient + number_diagnoses + insulin +
           change + diabetesMed
  D <- D[stats::complete.cases(D[, all.vars(f)]), ]
  set.seed(20260917)                                   # run_I_bigdata.R's split, unchanged
  dev <- sample(nrow(D), floor(nrow(D) / 2))
  fit <- stats::glm(f, data = D[dev, ], family = stats::binomial())
  Dval <- D[-dev, ]
  eta0 <- as.numeric(stats::predict(fit, newdata = Dval, type = "link"))
  b_lab <- unname(stats::coef(fit)[[RR_VAR]])
  ## the shortcut used in rr_one() is checked here against predict() on a hundred records, once
  chk <- sample.int(nrow(Dval), 100)
  Dc <- Dval[chk, , drop = FALSE]; Dc[[RR_VAR]] <- 4 * Dc[[RR_VAR]]
  ref <- as.numeric(stats::predict(fit, newdata = Dc, type = "response"))
  got <- stats::plogis(eta0[chk] + 3 * b_lab * Dval[[RR_VAR]][chk])
  if (max(abs(ref - got)) > 1e-10)
    stop("the linear-predictor shortcut disagrees with predict() by ", max(abs(ref - got)))
  list(fit = fit, Dval = Dval, eta0 = eta0, b_lab = b_lab, p0 = stats::plogis(eta0))
}

## the external-mode grouped quantities: nothing is estimated on this half, so Omega = I
rr_grp <- function(p, G) {
  ## paper 3: the rank rule of def.gof() in ebrahim.gof 2.9.0; paper 2's blocks used quantile break points
  if (identical(RR_COHORT, "p3")) return(as.integer(pmin(ceiling(rank(p, ties.method = "first") / (length(p) / G)), G)))
  br <- stats::quantile(p, seq(0, 1, length.out = G + 1), type = 1)
  br[1] <- -Inf; br[length(br)] <- Inf
  as.integer(cut(p, unique(br), include.lowest = TRUE))
}
## group sums by rowsum(), not tapply(): with two thousand groups over fifty thousand records the
## difference is seconds a replicate, and the arithmetic is the same
rr_sums <- function(y, p, g) {
  M <- rowsum(cbind(y, p, p * (1 - p)), g, reorder = TRUE)
  list(o = M[, 1], e = M[, 2], v = M[, 3], K = nrow(M))
}
rr_edge <- function(y, p, G, g = rr_grp(p, G)) {
  s <- rr_sums(y, p, g)
  r <- (s$o - s$e) / sqrt(s$v)
  pb <- pmin(pmax(s$e / as.numeric(table(g)), 1e-8), 1 - 1e-8)
  eta <- stats::qlogis(pb)
  Z <- if (identical(RR_COHORT, "p3")) cbind(1, pb, pb^2, pb^3) else cbind(1, eta, eta^2, eta^3)
  Q <- qr(Z); Z <- Z[, Q$pivot[seq_len(Q$rank)], drop = FALSE]
  S <- sum(qr.fitted(qr(Z), r)^2)
  stats::pchisq(S, ncol(Z), lower.tail = FALSE)
}
rr_hl <- function(y, p, G, g = rr_grp(p, G)) {
  s <- rr_sums(y, p, g)
  S <- sum(((s$o - s$e) / sqrt(s$v))^2)
  stats::pchisq(S, s$K, lower.tail = FALSE)
}
## Stukel's joint score on the frozen linear predictor, and the cubic calibration LR, both record-level
rr_stukel <- function(y, p) {
  eta <- stats::qlogis(pmin(pmax(p, 1e-8), 1 - 1e-8))
  z1 <- 0.5 * eta^2 * (eta >= 0); z2 <- -0.5 * eta^2 * (eta < 0)
  f0 <- stats::glm(y ~ offset(eta), family = stats::binomial())
  f1 <- stats::glm(y ~ offset(eta) + z1 + z2, family = stats::binomial())
  stats::anova(f0, f1, test = "Rao")$`Pr(>Chi)`[2]
}
rr_cubic <- function(y, p) {
  eta <- stats::qlogis(pmin(pmax(p, 1e-8), 1 - 1e-8))
  f0 <- stats::glm(y ~ eta, family = stats::binomial())
  f1 <- stats::glm(y ~ eta + I(eta^2) + I(eta^3), family = stats::binomial())
  stats::anova(f0, f1, test = "LRT")$`Pr(>Chi)`[2]
}

RR_TESTS <- c("EDGE.Grule", "EDGE.G10", "HL.G10", "Stk.joint", "Cubic.LR")

rr_one <- function(rep, k, P) {
  set.seed(910000000L + 1000L * k + rep)
  n <- nrow(P$Dval)
  y <- stats::rbinom(n, 1L, P$p0)                      # the frozen model is true for this replicate
  p <- P$p0
  if (k > 0L) {
    idx <- sample.int(n, k)
    ## the model is linear in num_lab_procedures, so multiplying that covariate by four adds exactly
    ## 3 * beta * x to the record's linear predictor. This is predict(newdata) to machine precision,
    ## without rebuilding a model frame of fifty thousand rows for every replicate; the equality is
    ## checked once at start-up below.
    p[idx] <- stats::plogis(P$eta0[idx] + 3 * P$b_lab * P$Dval[[RR_VAR]][idx])
  }
  G <- max(10L, as.integer(ceiling(n / 25)))
  out <- c(EDGE.Grule = rr_edge(y, p, G), EDGE.G10 = rr_edge(y, p, 10L), HL.G10 = rr_hl(y, p, 10L),
           Stk.joint = tryCatch(rr_stukel(y, p), error = function(e) NA_real_),
           Cubic.LR  = tryCatch(rr_cubic(y, p), error = function(e) NA_real_))
  c(rep = rep, k = k, n = n, events = sum(y), out)
}

OPT <- list(); a <- commandArgs(TRUE); i <- 1L
while (i <= length(a)) { OPT[[sub("^--", "", a[i])]] <- a[i + 1L]; i <- i + 2L }
W <- if (is.null(OPT$workers)) 20L else as.integer(OPT$workers)
REPS <- if (is.null(OPT$reps)) 500L else as.integer(OPT$reps)
OUT <- normalizePath(if (!is.null(OPT$out)) OPT$out else edge_battery("R"),
                     winslash = "/", mustWork = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lg <- function(...) { s <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", sprintf(...))
                      cat(s, "\n"); cat(s, "\n", file = file.path(OUT, "_progress.log"), append = TRUE) }

P <- rr_prepare()
lg("block R: validation half n=%d, event rate %.4f, rule G=%d, %d workers, %d replicates, covariate %s (coefficient %.4f)",
   nrow(P$Dval), mean(P$p0), max(10L, as.integer(ceiling(nrow(P$Dval) / 25))), W, REPS, RR_VAR, P$b_lab)

cl <- makePSOCKcluster(W)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
on.exit(stopCluster(cl))
clusterExport(cl, c("rr_one", "rr_grp", "rr_sums", "rr_edge", "rr_hl", "rr_stukel", "rr_cubic", "RR_COHORT",
                    "RR_TESTS", "RR_VAR", "P"), envir = environment())
invisible(clusterEvalQ(cl, Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")))

## P holds the cohort and the fitted glm, tens of megabytes once serialised. Passing it as an argument
## to clusterApplyLB sends the whole object with EVERY replicate, which cost about a hundred seconds a
## replicate against the one and a half seconds the tests actually take. The workers already have their
## own copy from clusterExport, so the task closure reads that instead; setting its environment to the
## global one keeps the closure itself small, because R does not serialise the global environment.
rr_task <- function(rep, k) rr_one(rep, k, P)
environment(rr_task) <- globalenv()

for (k in c(0L, 10L, 50L, 100L)) {
  f <- file.path(OUT, sprintf("real_C1_k%03d_pvalues.csv.gz", k))
  if (file.exists(f)) { lg("k=%d: present, skipped", k); next }
  t0 <- Sys.time()
  M <- as.data.table(do.call(rbind, clusterApplyLB(cl, seq_len(REPS), rr_task, k = k)))
  fwrite(M, f, compress = "gzip")
  lg("k=%-3d  %s  %.1f min", k,
     paste(sprintf("%s %.3f", RR_TESTS, vapply(RR_TESTS, function(t) mean(M[[t]] <= 0.05, na.rm = TRUE), 0)),
           collapse = "  "),
     as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
lg("block R: run finished")
