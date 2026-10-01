## pkg280_verify_corr_devtest.R -- devtools::test() on a clone of fix-stukel-joint-sym, with the totals.
## Run: Rscript pkg280_verify_corr_devtest.R <clone dir> > ../paper_EDGE/theory/pkg280_verify_corr_devtest.log 2>&1
args <- commandArgs(trailingOnly = TRUE)
PKG <- args[1]
cat("clone HEAD", system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE), "|", format(Sys.time()), "\n")
res <- devtools::test(PKG, reporter = testthat::SummaryReporter$new(show_praise = FALSE))
df <- as.data.frame(res)
cat(sprintf("\nblocks %d | expectations %d | failed %d | errors %d | skipped %d | warnings %d\n",
            nrow(df), sum(df$nb), sum(df$failed), sum(df$error), sum(df$skipped), sum(df$warning)))
print(aggregate(cbind(blocks = 1, nb, failed, skipped, warning, error = as.numeric(error)) ~ file, data = df, FUN = sum), row.names = FALSE)
cat("done", format(Sys.time()), "\n")
