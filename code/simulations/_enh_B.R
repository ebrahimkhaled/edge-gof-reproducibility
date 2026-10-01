# Enhancement B: DIRECTED (low-df) EF. Project group residuals onto a shape basis, test few df.
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
suppressMessages({library(parallel); library(ResourceSelection)})
PROJ<-edge_path("code", "legacy_not_deposited"); sw<-function(e) suppressWarnings(suppressMessages(e))
invisible(sw(source(file.path(PROJ,"Hosmer (H) (equal width interval).R")))); alpha<-0.05
inv_stukel<-function(ev,a1,a2){z<-numeric(length(ev));pos<-ev>=0
  if(any(pos))  z[pos]<- if(abs(a1)<1e-12) ev[pos]  else (-1+sqrt(pmax(0,1+2*a1*ev[pos])))/a1
  if(any(!pos)) z[!pos]<-if(abs(a2)<1e-12) ev[!pos] else (-1+sqrt(pmax(0,1+2*a2*ev[!pos])))/a2; plogis(z)}
gen<-function(scn,n){x<-runif(n,-3,3);d<-rbinom(n,1,.5);eta<-0.6*x+0.5*d
  p<-switch(scn, null=plogis(eta), cloglog=1-exp(-exp(eta)), loglog=exp(-exp(-eta)),
            stukel_heavy=inv_stukel(eta,-1,-1), stukel_asym=inv_stukel(eta,-1,1), cauchit=pcauchy(eta), stop("bad"))
  data.frame(x=x,d=d,y=rbinom(n,1,p))}
# group residuals at G deciles
grp_resid<-function(y,ph,G=10){ph<-pmin(pmax(ph,1e-6),1-1e-6); n<-length(y)
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)
  og<-tapply(y,grp,sum); eg<-tapply(ph,grp,sum); ng<-tapply(y,grp,length); pbar<-as.numeric(tapply(ph,grp,mean))
  Vg<-ng*pbar*(1-pbar); list(r=as.numeric((og-eg)/sqrt(Vg)), pbar=pbar)}
ef_omni<-function(y,ph,G=10){gr<-grp_resid(y,ph,G); Tef<-sum(gr$r^2)-sum((1-2*gr$pbar)*gr$r/sqrt(1)) # note: corr term in standardized units
  # use the standard EF: Pearson - correction, both in (o-e)/Vg units -> recompute properly:
  NA}
# proper EF omnibus (chi2_{G-2}) and directed projections share grp computation
ef_all<-function(y,ph,G=10){
  ph<-pmin(pmax(ph,1e-6),1-1e-6); n<-length(y); grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)
  og<-tapply(y,grp,sum); eg<-tapply(ph,grp,sum); ng<-tapply(y,grp,length); pbar<-as.numeric(tapply(ph,grp,mean))
  Vg<-ng*pbar*(1-pbar); oe<-as.numeric(og-eg)
  Tomni<-sum(oe^2/Vg)-sum((1-2*pbar)*oe/Vg); p_omni<-1-pchisq(Tomni,G-2)
  r<-oe/sqrt(Vg)                                   # standardized residuals
  proj<-function(Z){Z<-as.matrix(Z); Ztr<-crossprod(Z,r); S<-as.numeric(t(Ztr)%*%solve(crossprod(Z))%*%Ztr); 1-pchisq(S,ncol(Z))}
  eta<-qlogis(pbar)
  c(EF.omni=p_omni,
    EFd.poly1=proj(poly(pbar,1)), EFd.poly2=proj(poly(pbar,2)), EFd.poly3=proj(poly(pbar,3)),
    EFd.stk=proj(cbind(eta^2*(eta>=0), -eta^2*(eta<0))),
    EFd.stk3=proj(cbind(eta, eta^2*(eta>=0), -eta^2*(eta<0))))}
stukel_p<-function(fit){eta<-predict(fit); d<-fit$data; d$za<-0.5*eta^2*(eta>=0); d$zb<- -0.5*eta^2*(eta<0)
  fa<-suppressWarnings(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial())); pchisq(deviance(fit)-deviance(fa),df=2,lower.tail=FALSE)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_); if(is.null(v)||length(v)!=1||!is.finite(v)) NA_real_ else as.numeric(v)}
one_rep<-function(scn,n){dat<-gen(scn,n); fit<-suppressWarnings(glm(y~x+d,data=dat,family=binomial())); ph<-as.numeric(fitted(fit)); y<-dat$y
  ea<-tryCatch(ef_all(y,ph,10),error=function(e) rep(NA_real_,6)); names(ea)<-c("EF.omni","EFd.poly1","EFd.poly2","EFd.poly3","EFd.stk","EFd.stk3")
  c(sapply(ea,function(v) sc(v)), HL=sc(hoslem.test(y,ph,g=10)$p.value), HLeqw=sc(hosmer_lemeshow_equal_intervals(ph,y,num_groups=10)$p.value), Stukel=sc(stukel_p(fit)))}
scns<-c("null","cloglog","loglog","stukel_heavy","stukel_asym","cauchit"); ns<-c(500,1000); REPS<-2000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L)); clusterSetRNGStream(cl,20250911)
parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
clusterExport(cl,c("PROJ","sw","alpha","inv_stukel","gen","grp_resid","ef_all","stukel_p","sc","one_rep"))
invisible(clusterEvalQ(cl,{suppressMessages(library(ResourceSelection)); sw<-function(e) suppressWarnings(suppressMessages(e)); source(file.path(PROJ,"Hosmer (H) (equal width interval).R"))})); on.exit(stopCluster(cl))
cat(sprintf("%-13s %5s | %7s %8s %8s %8s %7s %8s | %5s %6s %6s\n","scenario","n","EF.omni","EFd.pol1","EFd.pol2","EFd.pol3","EFd.stk","EFd.stk3","HL","HLeqw","Stukel"))
out<-list()
for(scn in scns) for(n in ns){clusterExport(cl,c("scn","n"),envir=environment()); m<-parSapply(cl,1:REPS,function(i) one_rep(scn,n))
  pw<-rowMeans(m<alpha,na.rm=TRUE); out[[length(out)+1]]<-data.frame(scenario=scn,n=n,t(round(pw,4)))
  cat(sprintf("%-13s %5d | %7.3f %8.3f %8.3f %8.3f %7.3f %8.3f | %5.3f %6.3f %6.3f\n",scn,n,
      pw["EF.omni"],pw["EFd.poly1"],pw["EFd.poly2"],pw["EFd.poly3"],pw["EFd.stk"],pw["EFd.stk3"],pw["HL"],pw["HLeqw"],pw["Stukel"]))}
write.csv(do.call(rbind,out),"Enhancement_B_directed.csv",row.names=FALSE); cat("\nSaved Enhancement_B_directed.csv\n")
