## pkg280_verify_fix_devtest.R -- on a separate clean clone of the fixed branch: roxygen2::roxygenise() drift check
## (git status must stay empty), then devtools::test(). Nothing is installed.
## Run: Rscript pkg280_verify_fix_devtest.R > ../paper_EDGE/theory/pkg280_verify_fix_devtest.log 2>&1
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

VF <- file.path(tempdir(), "vfix")
P2 <- file.path(VF, "head2")
sha <- system2("git", c("-C", shQuote(P2), "rev-parse", "--short", "HEAD"), stdout = TRUE)
cat("clone", P2, "| HEAD", sha, "| roxygen2", as.character(packageVersion("roxygen2")), "| devtools",
    as.character(packageVersion("devtools")), "|", R.version.string, "\n")
suppressMessages(roxygen2::roxygenise(P2))
st <- system2("git", c("-C", shQuote(P2), "status", "--porcelain"), stdout = TRUE)
cat("git status --porcelain after roxygenise:", if (length(st)) paste0("\n  ", st) else "empty", "\n")
r <- devtools::test(P2, reporter = "summary", stop_on_failure = FALSE)
d <- as.data.frame(r)
cat("blocks", nrow(d), "; expectations", sum(d$nb), "; passed", sum(d$passed), "; failed", sum(d$failed),
    "; errors", sum(d$error), "; skipped", sum(d$skipped), "; warnings", sum(d$warning), "\n")
st2 <- system2("git", c("-C", shQuote(P2), "status", "--porcelain"), stdout = TRUE)
cat("git status --porcelain after tests:", if (length(st2)) paste0("\n  ", st2) else "empty", "\n")
