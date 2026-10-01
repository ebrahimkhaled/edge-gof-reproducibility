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
PROJ<-edge_path("code", "legacy_not_deposited"); alpha<-0.05
sw<-function(x)suppressWarnings(suppressMessages(x))
suppressMessages({library(ResourceSelection);library(MASS);library(dplyr);library(rpart);library(boot);library(catdata);library(aplore3)})
sw(source(file.path(PROJ,"pigeonheyse.R")));sw(source(file.path(PROJ,"Hosmer (H) (equal width interval).R")))
sw(source(file.path(PROJ,"Tsiatis.R")));sw(source(file.path(PROJ,"Xie.R")));sw(source(file.path(PROJ,"PR_test_only.R")))
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_);if(is.null(v)||length(v)!=1||!is.finite(v))NA_real_ else as.numeric(v)}
# link-general DEF (uses family's mu.eta so it is valid for cloglog fits too)
ef_dir_cal<-function(y,fit,basis,k,G=10){
  ph<-pmin(pmax(as.numeric(fitted(fit)),1e-6),1-1e-6);eta<-as.numeric(predict(fit,type="link"))
  dmu<-fit$family$mu.eta(eta);V<-ph*(1-ph);w<-dmu^2/V;n<-length(y)
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(V[I]))
  pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg);g0<-length(idx)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(dmu[I]*X[I,,drop=FALSE])))/sqrt(Vg)
  Om<-diag(g0)-U%*%solve(crossprod(X,w*X))%*%t(U)
  if(basis=="poly"){if(length(unique(round(pbar,8)))<k+1)return(NA_real_);Z<-as.matrix(poly(pbar,k))}else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  if(ncol(Z)<1)return(NA_real_);Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA_real_)
  cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}
ef_omni<-function(y,ph,G=10){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);g<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);gg<-length(o);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),gg-2)}
hleqw<-function(y,ph,B=10){g<-cut(ph,seq(0,1,length.out=B+1),include.lowest=TRUE,labels=FALSE);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);nn<-tapply(y,g,length);k<-!is.na(o);o<-o[k];e<-e[k];nn<-nn[k];st<-sum((o-e)^2/(e+1e-10)+((nn-o)-(nn-e))^2/((nn-e)+1e-10));if(length(o)>2)1-pchisq(st,length(o)-2) else NA_real_}
stuk<-function(fit){e<-predict(fit);d<-fit$data;d$za<-0.5*e^2*(e>=0);d$zb<- -0.5*e^2*(e<0);fa<-suppressWarnings(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial()));list(p=pchisq(deviance(fit)-deviance(fa),2,lower.tail=FALSE),bmax=max(abs(coef(fa)),na.rm=TRUE),conv=fa$converged)}
nm<-c("DEF.poly2","DEF.poly3","DEF.stk3","EF","HL","HLeqw","PH","Tsiatis","Xie","PR","Stukel")
run_all<-function(dat,form,fam=binomial(),catvar=NA,G=10){
  dat<-as.data.frame(dat);fit<-suppressWarnings(glm(form,data=dat,family=fam));ph<-as.numeric(fitted(fit));y<-dat$y
  islogit<-fit$family$link=="logit";f2<-fit;f2$predicted_probs<-ph;set.seed(7)
  st<-if(islogit) sc(stuk(fit)$p) else NA_real_
  # DEF and Stukel are diagnostics for a fitted *logit* model; for non-logit fits report only the link-agnostic tests
  d2<-if(islogit) sc(ef_dir_cal(y,fit,"poly",2,G)) else NA_real_
  d3<-if(islogit) sc(ef_dir_cal(y,fit,"poly",3,G)) else NA_real_
  ds<-if(islogit) sc(ef_dir_cal(y,fit,"stukel",3,G)) else NA_real_
  v<-c(DEF.poly2=d2,DEF.poly3=d3,DEF.stk3=ds,
    EF=sc(ef_omni(y,ph,G)),HL=sc(hoslem.test(y,ph,g=G)$p.value),HLeqw=sc(hleqw(y,ph,G)),
    PH=sc(pigeon_heyse_test(data.frame(y=y),fit,g=G)$p_value),Tsiatis=sc(score_gof_clustering(fit,num_groups=G,y=y)$p_value),
    Xie=sc(as.numeric(XieGoodnessOfFitTest(dat,f2))),PR=if(is.na(catvar))NA_real_ else sc(pr_test(dat,"y",catvar,ph)$p_value),Stukel=st)
  attr(v,"bmax")<-max(abs(coef(fit)),na.rm=TRUE);attr(v,"conv")<-fit$converged;attr(v,"event")<-mean(y);attr(v,"n")<-length(y);v}
ROWS<-list()
add<-function(label,v) ROWS[[length(ROWS)+1]]<<-data.frame(model=label,n=attr(v,"n"),event=round(attr(v,"event"),2),bmax=round(attr(v,"bmax"),1),conv=attr(v,"conv"),t(round(v,3)),check.names=FALSE)

## 1) BEETLE (Bliss 1935): logit is the wrong link; cloglog is correct ----
be<-data.frame(dose=c(1.6907,1.7242,1.7552,1.7842,1.8113,1.8369,1.8610,1.8839),n=c(59,60,62,56,63,59,62,60),k=c(6,13,18,28,52,53,61,60))
bdat<-do.call(rbind,Map(function(d,n,k) data.frame(dose=d,y=c(rep(1,k),rep(0,n-k))),be$dose,be$n,be$k))
add("Beetle: logit  y~dose",        run_all(bdat,y~dose,binomial("logit"),NA,8))
add("Beetle: cloglog y~dose [fix]",  run_all(bdat,y~dose,binomial("cloglog"),NA,8))

## 2) KYPHOSIS (rpart): non-monotone Age -> needs Age^2 ----
ky<-data.frame(y=as.integer(rpart::kyphosis$Kyphosis=="present"),Age=rpart::kyphosis$Age,Number=rpart::kyphosis$Number,Start=rpart::kyphosis$Start)
add("Kyphosis: additive",            run_all(ky,y~Age+Number+Start,binomial(),NA,8))
add("Kyphosis: + Age^2 [fix]",       run_all(transform(ky,Age2=scale(Age)^2),y~Age+Number+Start+Age2,binomial(),NA,8))

## 3) LOW BIRTH WEIGHT (HL): omitted interactions (mirror paper) ----
bw<-MASS::birthwt;lbw<-data.frame(y=bw$low,age=bw$age,lwt=bw$lwt,race=factor(bw$race),smoke=bw$smoke)
add("LBW: additive",                 run_all(lbw,y~age+lwt+race+smoke,binomial(),"smoke",10))
add("LBW: + age:lwt,smoke:lwt [fix]",run_all(transform(lbw,al=age*lwt,sl=smoke*lwt),y~age+lwt+race+smoke+al+sl,binomial(),"smoke",10))

## 4) NODAL involvement (Brown 1980; boot): sparse binary design ----
nd<-boot::nodal;nodf<-data.frame(y=nd$r,stage=nd$stage,grade=nd$grade,xray=nd$xray,acid=nd$acid,aged=nd$aged)
add("Nodal: 5 binary predictors",    run_all(nodf,y~stage+grade+xray+acid+aged,binomial(),"xray",6))

## 5) VASOCONSTRICTION (Finney 1947; catdata): separation-prone ----
data(vaso,package="catdata");va<-data.frame(y=vaso$vaso-min(vaso$vaso),vol=vaso$vol,rate=vaso$rate);va$y<-ifelse(va$y>0,1,0)
add("Vaso: y~vol+rate",              run_all(va,y~vol+rate,binomial(),NA,6))
add("Vaso: y~log(vol)+log(rate)[fix]",run_all(va,y~log(vol)+log(rate),binomial(),NA,6))

## 6) ICU (HL book; aplore3) ----
ic<-aplore3::icu;icd<-data.frame(y=as.integer(ic$sta=="Died"),age=ic$age,sys=ic$sys,typ=as.integer(ic$typ=="Emergency"),coma=as.integer(ic$loc!="Nothing"))
add("ICU: age+sys+typ+coma",         run_all(icd,y~age+sys+typ+coma,binomial(),"typ",10))

## 7) GLOW500 (HL book; aplore3) ----
gl<-aplore3::glow500;gld<-data.frame(y=as.integer(gl$fracture=="Yes"),age=gl$age,weight=gl$weight,priorfrac=as.integer(gl$priorfrac=="Yes"),momfrac=as.integer(gl$momfrac=="Yes"),raterisk=as.integer(gl$raterisk))
add("GLOW: age+weight+priorfrac+...", run_all(gld,y~age+weight+priorfrac+momfrac+raterisk,binomial(),"priorfrac",10))

OUT<-do.call(rbind,ROWS);rownames(OUT)<-NULL
write.csv(OUT,"RealData_direct.csv",row.names=FALSE)
print(OUT,row.names=FALSE)
cat("\nsaved RealData_direct.csv\n")
