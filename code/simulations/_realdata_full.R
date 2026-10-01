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
suppressMessages({library(ResourceSelection);library(CompQuadForm);library(MASS);library(dplyr)})
set.seed(1); PROJ<-edge_path("code", "legacy_not_deposited"); sw<-function(x)suppressWarnings(suppressMessages(x))
invisible(sw(source(file.path(PROJ,"pigeonheyse.R")))); invisible(sw(source(file.path(PROJ,"Hosmer (H) (equal width interval).R"))))
invisible(sw(source(file.path(PROJ,"Tsiatis.R")))); invisible(sw(source(file.path(PROJ,"Xie.R")))); invisible(sw(source(file.path(PROJ,"PR_test_only.R"))))
ef_dir_cal<-function(y,ph,fit,basis,k){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])));U<-U/sqrt(Vg);Om<-diag(G)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  if(basis=="poly")Z<-as.matrix(poly(pbar,k)) else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9]
  p<-tryCatch(CompQuadForm::imhof(S,lam)$Qq,error=function(e)NA);min(max(p,0),1)}
ef_omni<-function(y,ph){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);g<-pmin(ceiling(rank(ph,ties.method="first")/(n/10)),10)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),8)}
hleqw<-function(y,ph){g<-cut(ph,seq(0,1,length.out=11),include.lowest=TRUE,labels=FALSE);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);nn<-tapply(y,g,length);k<-!is.na(o);o<-o[k];e<-e[k];nn<-nn[k];st<-sum((o-e)^2/(e+1e-10)+((nn-o)-(nn-e))^2/((nn-e)+1e-10));if(length(o)>2)1-pchisq(st,length(o)-2) else NA}
stuk<-function(fit){e<-predict(fit);d<-fit$data;d$za<-0.5*e^2*(e>=0);d$zb<- -0.5*e^2*(e<0);fa<-suppressWarnings(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial()));pchisq(deviance(fit)-deviance(fa),2,lower.tail=FALSE)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_);if(is.null(v)||length(v)!=1||!is.finite(v))NA_real_ else as.numeric(v)}
sa<-read.csv(file.path(PROJ,"SA_Data.csv"))
rrep<-function(fit,lab){ph<-as.numeric(fitted(fit));y<-sa$y;dat<-sa;f2<-fit;f2$predicted_probs<-ph;saved<-get(".Random.seed",envir=.GlobalEnv)
  out<-c(DEF.poly2=sc(ef_dir_cal(y,ph,fit,"poly",2)),DEF.poly3=sc(ef_dir_cal(y,ph,fit,"poly",3)),DEF.stk3=sc(ef_dir_cal(y,ph,fit,"stukel",3)),
    EF=sc(ef_omni(y,ph)),HL=sc(hoslem.test(y,ph,g=10)$p.value),HLeqw=sc(hleqw(y,ph)),
    PH=sc(pigeon_heyse_test(data.frame(y=y),fit,g=10)$p_value),Tsiatis=sc(score_gof_clustering(fit,num_groups=10,y=y)$p_value),
    Xie=sc(as.numeric(XieGoodnessOfFitTest(dat,f2))),PR=sc(pr_test(dat,"y",c("race2","race3","smoke"),ph)$p_value),Stukel=sc(stuk(fit)))
  assign(".Random.seed",saved,envir=.GlobalEnv);cat("\n==",lab,"==\n");print(round(out,3))}
rrep(glm(y~age+lwt+race2+race3+smoke,data=sa,family=binomial()),"Additive")
rrep(glm(y~age+lwt+race2+race3+smoke+age:lwt+smoke:lwt,data=sa,family=binomial()),"With AGE:LWT + SMOKE:LWT")
