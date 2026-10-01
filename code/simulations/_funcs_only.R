# Benchmark EF vs ALL classical partition-based tests across null/link/quad scenarios.
# Tests: EF(chi2), HL, HL-equal-width, Pigeon-Heyse, Tsiatis, Xie, Pulkstenis-Robinson.
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
suppressMessages({library(parallel); library(ebrahim.gof); library(ResourceSelection)})
PROJ <- edge_path("code", "legacy_not_deposited")
sw <- function(e) suppressWarnings(suppressMessages(e))
load_tests <- function(){
  suppressMessages({library(ebrahim.gof);library(ResourceSelection);library(MASS);library(dplyr)})
  invisible(sw(source(file.path(PROJ,"pigeonheyse.R"))))
  invisible(sw(source(file.path(PROJ,"Hosmer (H) (equal width interval).R"))))
  invisible(sw(source(file.path(PROJ,"Tsiatis.R"))))
  invisible(sw(source(file.path(PROJ,"Xie.R"))))
  invisible(sw(source(file.path(PROJ,"PR_test_only.R"))))
}
load_tests()
G <- 10; alpha <- 0.05

gen <- function(scn,n){
  x<-runif(n,-3,3); d<-rbinom(n,1,.5); eta<-0.6*x+0.5*d
  p<-switch(scn, null=plogis(eta), cloglog=1-exp(-exp(eta)), loglog=exp(-exp(-eta)),
            probit=pnorm(eta), quad=plogis(eta+0.5*x^2-0.6), stop("bad"))
  data.frame(x=x,d=d,y=rbinom(n,1,p))
}
pv <- function(expr) tryCatch(expr, error=function(e) NA_real_)
one_rep <- function(scn,n){
  dat<-gen(scn,n); fit<-sw(glm(y~x+d,data=dat,family=binomial())); ph<-as.numeric(fitted(fit))
  z<-ef.gof(y=dat$y,predicted_probs=ph,G=G)$Test_Statistic; ef<-1-pchisq(z*sqrt(2*(G-2))+(G-2),G-2)
  f2<-fit; f2$predicted_probs<-ph
  c(EF      = ef,
    HL      = pv(hoslem.test(dat$y,ph,g=G)$p.value),
    HLeqw   = pv(hosmer_lemeshow_equal_intervals(ph,dat$y,num_groups=G)$p.value),
    PH      = pv(pigeon_heyse_test(data.frame(y=dat$y),fit,g=G)$p_value),
    Tsiatis = pv(score_gof_clustering(fit,num_groups=G,y=dat$y)$p_value),
    Xie     = pv(as.numeric(XieGoodnessOfFitTest(dat,f2))),
    PR      = pv(pr_test(dat,"y","d",ph)$p_value))
}

zz <- file(nullfile(), open="w"); sink(zz)
set.seed(123)
res <- lapply(c("null","cloglog","loglog","probit","quad"), function(scn)
  list(scn=scn, m=replicate(300, one_rep(scn, 1000))))
sink(); close(zz)
for (r in res) {
  pw <- round(rowMeans(r$m < 0.05, na.rm=TRUE), 3)
  na <- rowSums(is.na(r$m))
  cat(sprintf("%-8s power: %s | NA: %s\n", r$scn,
      paste(names(pw), pw, sep="=", collapse=" "), paste(na, collapse=",")))
}
