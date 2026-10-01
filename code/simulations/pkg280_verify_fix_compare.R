## pkg280_verify_fix_compare.R -- compares the outputs of pkg280_verify_fix_engine.R for three copies of ebrahim.gof:
##   head = fix-stukel-joint-sym after the review fixes (906bd0c), pre = before them (0d63fba), main = 2.7.0 (83c7e11).
## head vs pre: every output. head vs main: the unit form (poly2 / poly3 / stukel), the (y, ph, X) unit path, the
## ensemble, and head's form = "marginal" Stukel row against main's Stukel row (numbers only).
## Numbers agree when |a - b| <= 1e-10 * max(1, |a|). Every other case is listed.
## Run: Rscript pkg280_verify_fix_compare.R > ../paper_EDGE/theory/pkg280_verify_fix_compare.log 2>&1
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

VF  <- file.path(tempdir(), "vfix")
OUT <- edge_path("declarations")
options(width = 250, warn = 1)
H <- readRDS(file.path(VF, "engine_head.rds"))
P <- readRDS(file.path(VF, "engine_pre.rds"))
M <- readRDS(file.path(VF, "engine_main.rds"))
for (x in list(H, P, M)) cat(sprintf("%-5s ebrahim.gof %s HEAD %s, run %s\n", x$meta$label, x$meta$version, x$meta$head,
                                     format(x$meta$time)))

TOL <- 1e-10
NUM <- c("Test_Statistic", "Statistic", "df", "p_value")
cmp1 <- function(a, b, numeric_only = FALSE) {
  if (is.null(a) || is.null(b)) return(list(status = "missing", d = NA_real_))
  ea <- is.character(a); eb <- is.character(b)
  if (ea && eb) return(list(status = if (identical(a, b)) "both stop, same message" else "both stop, message differs",
                            d = NA_real_))
  if (ea) return(list(status = "old stops, new runs", d = NA_real_))
  if (eb) return(list(status = "old runs, new stops", d = NA_real_))
  if (identical(a, b)) return(list(status = "identical", d = 0))
  d <- 0; bad <- FALSE
  for (f in intersect(names(a), NUM)) {
    x <- suppressWarnings(as.numeric(a[[f]])); y <- suppressWarnings(as.numeric(b[[f]]))
    if (length(x) != length(y) || !identical(is.na(x), is.na(y))) { bad <- TRUE; d <- Inf; next }
    if (all(is.na(x))) next
    dd <- abs(x - y); d <- max(d, dd, na.rm = TRUE)
    if (any(dd > TOL * pmax(1, abs(x)), na.rm = TRUE)) bad <- TRUE
  }
  oth <- setdiff(union(names(a), names(b)), NUM)
  same_oth <- numeric_only || all(vapply(oth, function(f) identical(a[[f]], b[[f]]), logical(1)))
  list(status = if (bad) "DIFFER" else if (same_oth) "equal to 1e-10" else "numbers equal, labels differ", d = d)
}
kclass <- function(k) {
  parts <- strsplit(k, "|", fixed = TRUE)[[1]]
  parts <- vapply(parts, function(p) if (grepl("^G[0-9]+$", p)) (if (p == "G10") "G10" else "Gauto") else p, "")
  paste(parts, collapse = "|")
}
fmt <- function(x) {
  if (is.character(x)) return(substr(x, 1, 90))
  g <- function(f) { v <- x[[f]]; if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(v[1])) }
  st <- if (!is.null(x$Test_Statistic)) g("Test_Statistic") else g("Statistic")
  sprintf("S %.10g | df %.6g | p %.10g", st, g("df"), g("p_value"))
}

GROUPS <- c("ref9", "exa7", "base", "c0m2", "s2", "extra")
ROWS <- list(); CASES <- list()
do_cmp <- function(label, OLD, NEW, keymap) {
  for (g in GROUPS) {
    nms <- names(NEW[[g]])
    st <- character(0); kc <- character(0); dd <- numeric(0)
    for (i in seq_along(NEW[[g]])) {
      nr <- NEW[[g]][[i]]; or <- OLD[[g]][[i]]
      for (k in names(nr)) {
        if (startsWith(k, "pbar|")) next
        ok <- keymap(k); if (is.na(ok)) next
        r <- cmp1(or[[ok]], nr[[k]], numeric_only = startsWith(k, "stukel|"))
        st <- c(st, r$status); kc <- c(kc, kclass(k)); dd <- c(dd, r$d)
        if (!r$status %in% c("identical", "equal to 1e-10", "both stop, same message")) {
          G  <- if (grepl("\\|G[0-9]+", k)) sub(".*\\|G([0-9]+).*", "\\1", k) else "10"
          pb <- nr[[paste0("pbar|G", G)]]
          CASES[[length(CASES) + 1]] <<- data.frame(
            cmp = label, group = g, fit = if (is.null(nms)) as.character(i) else nms[i], key = k, status = r$status,
            maxdiff = r$d, n_groups_ge_half = pb[["n_up"]], min_pbar_ge_half = pb[["min_up"]],
            max_pbar_ge_half = pb[["max_up"]], max_pbar_lt_half = pb[["max_lo"]],
            old = fmt(or[[ok]]), new = fmt(nr[[k]]), stringsAsFactors = FALSE)
        }
      }
    }
    ROWS[[length(ROWS) + 1]] <<- data.frame(cmp = label, group = g, key = kc, status = st, d = dd,
                                            stringsAsFactors = FALSE)
  }
}
do_cmp("head vs pre (0d63fba)", P, H, function(k) k)
do_cmp("head vs main (2.7.0)", M, H, function(k) {
  if (grepl("^(unit\\|(poly2|poly3|stukel)\\||yphX unit\\||ensemble\\|)", k)) k
  else if (k == "stukel|marginal") "stukel|default" else NA_character_
})
RR <- do.call(rbind, ROWS)

cat("\n==== status counts by comparison and output (all fit groups pooled) ====\n")
tab <- as.data.frame(table(cmp = RR$cmp, key = RR$key, status = RR$status), stringsAsFactors = FALSE)
tab <- tab[tab$Freq > 0, ]
mx <- aggregate(d ~ cmp + key + status, data = transform(RR, d = ifelse(is.finite(RR$d), RR$d, NA)), FUN = max,
                na.action = na.pass)
tab <- merge(tab, mx, all.x = TRUE)
tab <- tab[order(tab$cmp, tab$key, tab$status), ]
names(tab)[names(tab) == "d"] <- "max_abs_diff"
print(tab, row.names = FALSE)

cat("\n==== status counts by comparison and fit group ====\n")
print(as.data.frame.matrix(table(paste(RR$cmp, RR$group, sep = " | "), RR$status)))

cat("\n==== cases that are not identical / equal to 1e-10 / both-stop-same ====\n")
if (length(CASES)) {
  CC <- do.call(rbind, CASES)
  print(CC, row.names = FALSE)
  write.csv(CC, file.path(OUT, "pkg280_verify_fix_compare_cases.csv"), row.names = FALSE)
  cat("\nby comparison x status x key:\n")
  print(as.data.frame(table(cmp = CC$cmp, status = CC$status, key = vapply(CC$key, kclass, "")),
                      stringsAsFactors = FALSE) |> subset(Freq > 0), row.names = FALSE)
} else cat("none\n")

cat("\n==== the stukel window in the c0 = -2 fits (G = 20 and G = 10) ====\n")
for (g in c("c0m2", "extra")) for (G in c(10, 20)) {
  pb <- t(vapply(H[[g]], function(r) r[[paste0("pbar|G", G)]], numeric(4)))
  win <- pb[, "n_up"] >= 1 & pb[, "max_up"] < 0.5016
  cat(sprintf("  %-5s G %2d: fits %d | with a group mean >= 0.5: %d | all such groups below 0.5016: %d | below 0.500167: %d\n",
              g, G, nrow(pb), sum(pb[, "n_up"] >= 1), sum(win), sum(pb[, "n_up"] >= 1 & pb[, "max_up"] < 0.500167)))
}
cat("\ndone\n")
