## pkg280_verify_corr_rox.R -- are man/*.Rd and NAMESPACE what roxygen2 generates from the branch HEAD?
## roxygen2::roxygenise() runs on a throw-away clone; git then lists any file it changed.
## Run: Rscript pkg280_verify_corr_rox.R <throw-away clone dir> > ../paper_EDGE/theory/pkg280_verify_corr_rox.log 2>&1
args <- commandArgs(trailingOnly = TRUE)
PKG <- args[1]
cat("clone HEAD", system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE),
    "| roxygen2", as.character(packageVersion("roxygen2")), "| RoxygenNote",
    read.dcf(file.path(PKG, "DESCRIPTION"), fields = "RoxygenNote")[1], "|", format(Sys.time()), "\n")
suppressMessages(roxygen2::roxygenise(PKG))
st <- system2("git", c("-C", shQuote(PKG), "status", "--porcelain"), stdout = TRUE)
cat("files changed by roxygenise():", length(st), "\n")
if (length(st)) { cat(st, sep = "\n"); cat(system2("git", c("-C", shQuote(PKG), "diff", "--stat"), stdout = TRUE), sep = "\n") }
cat("done", format(Sys.time()), "\n")
