## pkg280_verify_corr_install.R -- is the installed ebrahim.gof the code of the branch HEAD?
## Every top-level assignment in R/*.R of a clone of fix-stukel-joint-sym is parsed without source references and
## compared, by deparse, with the object of the same name in the installed namespace.
## Verification of the package, not a paper result.
## Run: Rscript pkg280_verify_corr_install.R <clone dir> > ../paper_EDGE/theory/pkg280_verify_corr_install.log 2>&1

args <- commandArgs(trailingOnly = TRUE)
PKG  <- args[1]
options(width = 200)
suppressPackageStartupMessages(library(ebrahim.gof))
ns  <- asNamespace("ebrahim.gof")
sha <- system2("git", c("-C", shQuote(PKG), "rev-parse", "--short", "HEAD"), stdout = TRUE)
cat("installed ebrahim.gof", as.character(packageVersion("ebrahim.gof")), "| Packaged",
    packageDescription("ebrahim.gof")$Packaged, "| lib", dirname(system.file(package = "ebrahim.gof")), "\n")
cat("clone HEAD", sha, "|", R.version.string, "|", format(Sys.time()), "\n\n")

n_same <- 0L; diffs <- character(0); missing <- character(0)
for (f in list.files(file.path(PKG, "R"), pattern = "[.]R$", full.names = TRUE)) {
  for (ex in parse(f, keep.source = FALSE)) {
    if (!(is.call(ex) && identical(ex[[1]], as.name("<-")) && is.name(ex[[2]]))) next
    nm <- as.character(ex[[2]])
    e  <- new.env(parent = ns)
    val <- tryCatch(eval(ex[[3]], e), error = function(err) NULL)
    if (is.null(val)) next
    if (!exists(nm, envir = ns, inherits = FALSE)) { missing <- c(missing, nm); next }
    inst <- get(nm, envir = ns, inherits = FALSE)
    a <- paste(deparse(val, control = c("keepInteger", "niceNames")), collapse = "\n")
    b <- paste(deparse(inst, control = c("keepInteger", "niceNames")), collapse = "\n")
    if (identical(a, b)) n_same <- n_same + 1L else diffs <- c(diffs, paste0(basename(f), ": ", nm))
  }
}
cat("objects identical by deparse:", n_same, "\n")
cat("objects that differ:", length(diffs), "\n"); if (length(diffs)) cat(paste0("  ", diffs), sep = "\n")
cat("objects not in the installed namespace:", length(missing), "\n"); if (length(missing)) cat(paste0("  ", missing), sep = "\n")
cat("\nkey rules in the installed code:\n")
dg <- paste(deparse(get("def.gof", ns)), collapse = "\n")
gs <- paste(deparse(get("gof_stukel", ns)), collapse = "\n")
chk <- c("2.7.0 filter"        = grepl("colSums(abs(Z)) > 1e-08", dg, fixed = TRUE),
         "unit-length scaling" = grepl("Z/rep(sqrt(colSums(Z^2)), each = nrow(Z))", dg, fixed = TRUE),
         "score guard"         = grepl("diag(I) > 1e-10 * colSums(Zs^2)", dg, fixed = TRUE),
         "degenerate check"    = grepl("min(sum(y), n - sum(y)) == 0", dg, fixed = TRUE),
         "Stukel joint guard"  = grepl("diag(I) > 1e-10 * colSums(W * Z^2)", gs, fixed = TRUE),
         "Stukel raw risks"    = grepl("ph <- as.numeric(stats::fitted(ctx$model))", gs, fixed = TRUE),
         "no relative rule"    = !grepl("1e-06 * max(nz)", dg, fixed = TRUE))
print(chk)
cat("\ndone", format(Sys.time()), "\n")
