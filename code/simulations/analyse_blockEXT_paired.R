## analyse_blockEXT_paired.R -- paired differences in size-adjusted power, EDGE at ten groups minus Cox's test and
## minus the GiViTI belt, per alternative cell of block EXT, with Monte Carlo standard errors on the same replicates.
## Added after the referee round of 2026-10-01; the critical values are those of analyse_blockEXT.R (cell 2).
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
suppressMessages(library(data.table))
IN <- edge_battery("EXT")
rd <- function(id) fread(list.files(IN, sprintf("^c%02d_.*gz$", id), full.names = TRUE))
crit <- rd(2)[, lapply(.SD, function(p) as.numeric(quantile(p, 0.05, type = 1))), .SDcols = c("EDGE_G10", "cox", "giviti")]
FAM <- c(`4` = "shift 0.2", `5` = "shift 0.4", `6` = "slope 0.8", `7` = "slope 0.6", `8` = "U-shape 0.4",
         `9` = "U-shape 0.8", `10` = "threshold 2", `11` = "threshold 4", `12` = "interaction 0.6", `13` = "interaction 1")
R <- rbindlist(lapply(4:13, function(id) {
  d <- rd(id)
  e <- d$EDGE_G10 <= crit$EDGE_G10; cx <- d$cox <= crit$cox; gv <- !is.na(d$giviti) & d$giviti <= crit$giviti
  data.table(id = id, cell = FAM[[as.character(id)]],
             d_cox = mean(e - cx), se_cox = sd(e - cx) / sqrt(nrow(d)),
             d_giviti = mean(e - gv), se_giviti = sd(e - gv) / sqrt(nrow(d)),
             same_p_cox_giviti = mean(abs(d$giviti - d$cox) < 1e-12, na.rm = TRUE))
}))
print(R, digits = 3)
cat(sprintf("\nmean difference vs Cox %.3f, vs belt %.3f; without the interaction cells %.3f and %.3f\n",
            mean(R$d_cox), mean(R$d_giviti), mean(R[id < 12, d_cox]), mean(R[id < 12, d_giviti])))
fwrite(R, file.path(IN, "_paired_EDGE_G10.csv"))
