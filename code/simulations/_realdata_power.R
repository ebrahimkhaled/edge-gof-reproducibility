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
suppressMessages({library(parallel)}); PROJ<-edge_path("code", "legacy_not_deposited"); alpha<-0.05
sw<-function(x)suppressWarnings(suppressMessages(x))
load_tests<-function(){suppressMessages({library(ResourceSelection);library(MASS);library(dplyr);library(aplore3)})
  invisible(sw(source(file.path(PROJ,"pigeonheyse.R"))));invisible(sw(source(file.path(PROJ,"Hosmer (H) (equal width interval).R"))))
  invisible(sw(source(file.path(PROJ,"Tsiatis.R"))));invisible(sw(source(file.path(PROJ,"Xie.R"))));invisible(sw(source(file.path(PROJ,"PR_test_only.R"))));invisible(NULL)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_);if(is.null(v)||length(v)!=1||!is.finite(v))NA_real_ else as.numeric(v)}
ef_dir_cal<-function(y,ph,fit,basis,k,G=10){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y)
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg);g0<-length(idx)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])))/sqrt(Vg);Om<-diag(g0)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  if(basis=="poly")Z<-as.matrix(poly(pbar,k)) else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  if(ncol(Z)<1)return(NA_real_);Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA_real_);cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}
ef_omni<-function(y,ph,G=10){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);g<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);gg<-length(o);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),gg-2)}
hleqw<-function(y,ph,B=10){g<-cut(ph,seq(0,1,length.out=B+1),include.lowest=TRUE,labels=FALSE);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);nn<-tapply(y,g,length);k<-!is.na(o);o<-o[k];e<-e[k];nn<-nn[k];st<-sum((o-e)^2/(e+1e-10)+((nn-o)-(nn-e))^2/((nn-e)+1e-10));if(length(o)>2)1-pchisq(st,length(o)-2) else NA_real_}
stuk<-function(fit){e<-predict(fit);d<-fit$data;d$za<-0.5*e^2*(e>=0);d$zb<- -0.5*e^2*(e<0);fa<-suppressWarnings(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial()));pchisq(deviance(fit)-deviance(fa),2,lower.tail=FALSE)}
nm<-c("DEF.poly2","DEF.poly3","DEF.stk3","EF","HL","HLeqw","PH","Tsiatis","Xie","PR","Stukel")
# ---- build the real GLOW design (fixed covariates) ----
gl<-aplore3::glow500; GD<-data.frame(y=as.integer(gl$fracture=="Yes"),age=gl$age,weight=gl$weight,priorfrac=as.integer(gl$priorfrac=="Yes"),momfrac=as.integer(gl$momfrac=="Yes"))
fit0<-glm(y~age+weight+priorfrac+momfrac,data=GD,family=binomial());eta0<-as.numeric(predict(fit0));pbar<-mean(GD$y)
z<-as.numeric(scale(GD$age));etac<-eta0-mean(eta0)
# cloglog truth with a steeper index (so the link departure is genuinely present, not logit-equivalent), intercept matched to prevalence
SLINK<-1.9
acll<-uniroot(function(a) mean(1-exp(-exp(a+SLINK*etac)))-pbar,c(-15,15))$root
truth<-function(kind,g_nl,g_int){
  if(kind=="null") plogis(eta0)
  else if(kind=="link") 1-exp(-exp(acll+SLINK*etac))
  else if(kind=="nonlin") plogis(eta0+g_nl*(z^2-1))
  else plogis(eta0+g_int*GD$priorfrac*z)}
one_rep<-function(kind,g_nl,g_int){p<-truth(kind,g_nl,g_int);y<-rbinom(length(p),1,p);d<-GD;d$y<-y
  fit<-suppressWarnings(glm(y~age+weight+priorfrac+momfrac,data=d,family=binomial()));ph<-as.numeric(fitted(fit))
  saved<-get(".Random.seed",envir=.GlobalEnv);f2<-fit;f2$predicted_probs<-ph
  v<-c(DEF.poly2=sc(ef_dir_cal(y,ph,fit,"poly",2)),DEF.poly3=sc(ef_dir_cal(y,ph,fit,"poly",3)),DEF.stk3=sc(ef_dir_cal(y,ph,fit,"stukel",3)),
    EF=sc(ef_omni(y,ph)),HL=sc(hoslem.test(y,ph,g=10)$p.value),HLeqw=sc(hleqw(y,ph)),
    PH=sc(pigeon_heyse_test(data.frame(y=y),fit,g=10)$p_value),Tsiatis=sc(score_gof_clustering(fit,num_groups=10,y=y)$p_value),
    Xie=sc(as.numeric(XieGoodnessOfFitTest(d,f2))),PR=sc(pr_test(d,"y","priorfrac",ph)$p_value),Stukel=sc(stuk(fit)))
  assign(".Random.seed",saved,envir=.GlobalEnv);v}
cct<-function(P){p<-pmin(pmax(P,1e-15),1-1e-15);0.5-atan(mean(tan((0.5-p)*pi),na.rm=TRUE))/pi}
REPS<-2000;g_nl<-0.85;g_int<-1.6
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L));clusterSetRNGStream(cl,7770)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl,c("PROJ","sw","load_tests","alpha","sc","ef_dir_cal","ef_omni","hleqw","stuk","nm","GD","eta0","etac","z","pbar","acll","SLINK","truth","one_rep","g_nl","g_int"))
invisible(clusterEvalQ(cl,load_tests()));on.exit(stopCluster(cl))
kinds<-c("null","link","nonlin","int");pl<-list();powtab<-list()
for(kd in kinds){clusterExport(cl,"kd",envir=environment());m<-parSapply(cl,1:REPS,function(i)one_rep(kd,g_nl,g_int))
  df<-as.data.frame(t(m[nm,]));names(df)<-nm
  # ensembles from per-rep p-values
  Vote<-apply(df[,c("DEF.poly2","DEF.poly3","DEF.stk3")],1,cct);Univ<-apply(df[,c("DEF.poly2","DEF.poly3","DEF.stk3","EF","Tsiatis")],1,cct)
  df$DEF.vote<-Vote;df$DEF.univ<-Univ;df<-cbind(kind=kd,rep=1:REPS,df);pl[[length(pl)+1]]<-df
  pw<-sapply(c(nm,"DEF.vote","DEF.univ"),function(c) mean(df[[c]]<0.05,na.rm=TRUE));powtab[[kd]]<-pw;cat(kd,"done\n")}
PW<-do.call(rbind,powtab)
write.csv(do.call(rbind,pl),"RealPower_pvalues.csv",row.names=FALSE)
write.csv(data.frame(kind=rownames(PW),round(PW,3)),"RealPower_power.csv",row.names=FALSE)
cat("\n== Semi-synthetic power on real GLOW covariates (n=500, ",REPS," reps) ==\n",sep="");print(round(PW,3))
cat("\nsaved RealPower_pvalues.csv, RealPower_power.csv\n")
