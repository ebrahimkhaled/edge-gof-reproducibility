# WP-B1 — n x G sweep for the EDGES ensemble paper.
# Parameterized copy of _ensemble.R: the helper blocks (lg,PROJ,sw,load_tests,inv_stukel,
# sqb,linkp,linkscns,gen,Zp2,Zp3,Zsk,Zcb,ef_omni,sc,nm,scns) are copied VERBATIM from
# _ensemble.R (lines 5-39,45); ef_omni is parameterized for G and one_rep takes (scn,n,G).
# Do NOT source("_ensemble.R") -- that file auto-runs the whole 1000/10 job at its bottom.
# Stores per-rep component p-values (poly2,poly3,stk,combo,EF,Tsiatis) to wp_b1_pvalues.csv.
# Combiners (minP/Fisher/Stouffer/HMP/CCT) are formed afterwards in wp_b1_combine.py.
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
suppressMessages(library(parallel)); alpha<-0.05; lg<-function(p) log(p/(1-p))
PROJ<-edge_path("code", "legacy_not_deposited"); sw<-function(x)suppressWarnings(suppressMessages(x))
load_tests<-function(){suppressMessages({library(ResourceSelection);library(MASS);library(dplyr)})
  invisible(sw(source(file.path(PROJ,"Tsiatis.R"))));invisible(NULL)}
inv_stukel<-function(ev,a1,a2){z<-numeric(length(ev));pos<-ev>=0
  if(any(pos)) z[pos]<- if(abs(a1)<1e-12) ev[pos] else (-1+sqrt(pmax(0,1+2*a1*ev[pos])))/a1
  if(any(!pos))z[!pos]<-if(abs(a2)<1e-12) ev[!pos] else (-1+sqrt(pmax(0,1+2*a2*ev[!pos])))/a2; plogis(z)}
sqb<-function(J) solve(rbind(c(1,-1.5,2.25),c(1,3,9),c(1,-3,9)),c(lg(.05),lg(.95),lg(J)))
linkp<-function(scn,eta) switch(scn, cloglog=1-exp(-exp(eta)), loglog=exp(-exp(-eta)), cauchit=pcauchy(eta),
  stukel_heavy=inv_stukel(eta,-1,-1), stukel_asym=inv_stukel(eta,-1,1))
linkscns<-c("cloglog","loglog","cauchit","stukel_heavy","stukel_asym")
gen<-function(scn,n){
  if(scn=="null"){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d))),f=y~x+d)}
  else if(scn%in%linkscns){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,linkp(scn,0.6*x+0.5*d))),f=y~x+d)}
  else if(scn=="quad"){b<-sqb(0.02);x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*x^2))),f=y~x)}
  else if(scn=="cubic"){x<-runif(n,-2.5,2.5);list(d=data.frame(x=x,y=rbinom(n,1,plogis(0.3+0.7*x+0.12*(x^3-3*x)))),f=y~x)}
  else if(scn=="sawtooth"){x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(0.8*x+1.2*(2*(x/1.5-floor(x/1.5+0.5)))))),f=y~x)}
  else {x<-runif(n,-2.5,2.5);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.9*x-0.5*x*d))),f=y~x+d)}}  # crossover
# general DEF test: pass a function Zfun(pbar, etabar) returning the basis matrix
def_test<-function(y,ph,fit,Zfun,G=10){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y)
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])))/sqrt(Vg);Om<-diag(length(idx))-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  eb<-qlogis(pbar);Z<-tryCatch(Zfun(pbar,eb),error=function(e)NULL);if(is.null(Z))return(NA_real_)
  Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE];if(ncol(Z)<1)return(NA_real_)
  Zi<-tryCatch(solve(crossprod(Z)),error=function(e)NULL);if(is.null(Zi))return(NA_real_)
  Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA_real_)
  cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}
Zp2<-function(pb,eb) as.matrix(poly(pb,2)); Zp3<-function(pb,eb) as.matrix(poly(pb,3))
Zsk<-function(pb,eb) cbind(eb,eb^2*(eb>=0),-eb^2*(eb<0)); Zcb<-function(pb,eb) cbind(poly(pb,3),eb^2*(eb>=0),-eb^2*(eb<0))
# ef_omni PARAMETERIZED for G (BEFORE/AFTER per WP-B1.1 step 2): (n/10)->(n/G), ,10)->,G), ,8)->,G-2)
ef_omni<-function(y,ph,G=10){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);g<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),G-2)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_);if(is.null(v)||length(v)!=1||!is.finite(v))NA_real_ else as.numeric(v)}
nm<-c("poly2","poly3","stk","combo","EF","Tsiatis")
# one_rep PARAMETERIZED for (scn,n,G): pass G into the 4 def_test calls and ef_omni; .Random.seed
# save/restore around the 6 test calls kept VERBATIM from _ensemble.R lines 41,44.
one_rep<-function(scn,n,G){G0<-gen(scn,n);dat<-G0$d;fit<-suppressWarnings(glm(G0$f,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  saved<-get(".Random.seed",envir=.GlobalEnv)
  v<-c(poly2=sc(def_test(y,ph,fit,Zp2,G=G)),poly3=sc(def_test(y,ph,fit,Zp3,G=G)),stk=sc(def_test(y,ph,fit,Zsk,G=G)),
       combo=sc(def_test(y,ph,fit,Zcb,G=G)),EF=sc(ef_omni(y,ph,G=G)),Tsiatis=sc(score_gof_clustering(fit,num_groups=G,y=y)$p_value))
  assign(".Random.seed",saved,envir=.GlobalEnv);v}
scns<-c("null","cloglog","loglog","cauchit","stukel_heavy","stukel_asym","quad","cubic","sawtooth","crossover")

## --- cluster core cap (copied from _harness.R:9-10) ---
EK_NCORES <- max(1L, detectCores(logical = TRUE) - 1L)
{ .cap <- suppressWarnings(as.integer(Sys.getenv("EK_NCORES_CAP", ""))); if (!is.na(.cap) && .cap >= 1L) EK_NCORES <- min(EK_NCORES, .cap) }

## --- REPS_SCALE smoke helper (copied from _harness.R:13) ---
ek_reps <- function(B) max(2L, as.integer(round(B * as.numeric(Sys.getenv("REPS_SCALE", "1")))))

## --- incremental append (copied from _harness.R:25-28) ---
ek_append <- function(df, path) {
  ex <- file.exists(path)
  write.table(df, path, sep = ",", row.names = FALSE, col.names = !ex, append = ex, qmethod = "double")
}

OUT <- edge_path("code/simulations/wp_b1_pvalues.csv")
if (file.exists(OUT)) file.remove(OUT)   # fresh run: header written on first append

## --- cell grid (LOCKED, WP-B1.1 step 4) ---
cells <- rbind(
  data.frame(n=200,G=10), data.frame(n=500,G=10),
  data.frame(n=1000,G=10), data.frame(n=2000,G=10),
  data.frame(n=1000,G=8),  data.frame(n=1000,G=12))
GRID_SEED <- 20260711L

cat(sprintf("WP-B1 n x G sweep starting: %d cells, EK_NCORES=%d, REPS_SCALE=%s\n",
            nrow(cells), EK_NCORES, Sys.getenv("REPS_SCALE","1")))
t0 <- Sys.time()

for (cid in seq_len(nrow(cells))) {
  n <- cells$n[cid]; G <- cells$G[cid]
  cl <- makeCluster(EK_NCORES)
  parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
  clusterSetRNGStream(cl, GRID_SEED + cid*1000L)   # distinct stream per cell
  clusterExport(cl, c("PROJ","sw","load_tests","alpha","lg","inv_stukel","sqb","linkp","linkscns",
                      "gen","def_test","Zp2","Zp3","Zsk","Zcb","ef_omni","sc","one_rep","n","G","nm"),
                envir = environment())
  invisible(clusterEvalQ(cl, load_tests()))
  for (scn in scns) {
    REPS <- ek_reps(if (scn=="null") 5000L else 2000L)
    clusterExport(cl, "scn", envir = environment())
    m <- parSapply(cl, 1:REPS, function(i) one_rep(scn, n, G))
    df <- as.data.frame(t(m[nm,,drop=FALSE])); names(df) <- nm
    df <- cbind(n=n, G=G, scenario=scn, rep=1:REPS, df)
    ek_append(df, OUT)
    cat(sprintf("  cell(n=%d,G=%d) %-13s done  reps=%d  elapsed=%.1f min\n",
                n, G, scn, REPS, as.numeric(difftime(Sys.time(), t0, units="mins"))))
  }
  stopCluster(cl)
}
cat(sprintf("WP-B1 DONE. total wall-clock=%.2f min -> %s\n",
            as.numeric(difftime(Sys.time(), t0, units="mins")), OUT))
