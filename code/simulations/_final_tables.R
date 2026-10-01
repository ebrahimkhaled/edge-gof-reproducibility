# Publication-grade run: core 7 tests, 14 scenarios x n in {200,500,1000,2000}, 10000 reps.
# Stores per-replication p-values (DEF_pvalues_all.csv) + power summary (DEF_power_summary.csv).
suppressMessages(library(parallel)); alpha<-0.05; lg<-function(p) log(p/(1-p))
inv_stukel<-function(ev,a1,a2){z<-numeric(length(ev));pos<-ev>=0
  if(any(pos)) z[pos]<- if(abs(a1)<1e-12) ev[pos] else (-1+sqrt(pmax(0,1+2*a1*ev[pos])))/a1
  if(any(!pos))z[!pos]<-if(abs(a2)<1e-12) ev[!pos] else (-1+sqrt(pmax(0,1+2*a2*ev[!pos])))/a2; plogis(z)}
sqb<-function(J) solve(rbind(c(1,-1.5,2.25),c(1,3,9),c(1,-3,9)),c(lg(.05),lg(.95),lg(J)))
sib<-function(I){e0<-lg(.1);e1<-lg(.2);e2<-lg(.2+I);b0<-(e0+e1)/2;b1<-(e1-e0)/6;b3<-(e2-(b0+3*b1))/6;c(b0,b1,3*b3,b3)}
scb<-function(K) solve(rbind(c(1,-3,-2,6),c(1,-3,0,0),c(1,3,0,0),c(1,3,2,6)),c(lg(.1),lg(.1),lg(.2),lg(.2+K)))
linkp<-function(scn,eta) switch(scn, cloglog=1-exp(-exp(eta)), loglog=exp(-exp(-eta)), probit=pnorm(eta), cauchit=pcauchy(eta),
  robit_t4=pt(eta,df=4), scobit2=plogis(eta)^2, scobit_half=plogis(eta)^0.5,
  stukel_heavy=inv_stukel(eta,-1,-1), stukel_light=inv_stukel(eta,1,1), stukel_asym=inv_stukel(eta,-1,1))
linkscns<-c("cloglog","loglog","probit","cauchit","robit_t4","scobit2","scobit_half","stukel_heavy","stukel_light","stukel_asym")
gen<-function(scn,n){
  if(scn=="null"){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d))),f=y~x+d)}
  else if(scn%in%linkscns){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,linkp(scn,0.6*x+0.5*d))),f=y~x+d)}
  else if(scn=="omit_quad"){b<-sqb(0.02);x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*x^2))),f=y~x)}
  else if(scn=="omit_int_bin"){b<-sib(0.5);x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*d+b[4]*x*d))),f=y~x+d)}
  else {b<-scb(0.5);x<-runif(n,-3,3);z<-rnorm(n);list(d=data.frame(x=x,z=z,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*z+b[4]*x*z))),f=y~x+z)}}
ef_dir_cal<-function(y,ph,fit,basis,k){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])));U<-U/sqrt(Vg);Om<-diag(G)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  if(basis=="poly")Z<-as.matrix(poly(pbar,k)) else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  if(ncol(Z)<1)return(NA_real_);Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA_real_)
  cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}
ef_omni<-function(y,ph){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);g<-pmin(ceiling(rank(ph,ties.method="first")/(n/10)),10)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),8)}
hleqw<-function(y,ph){g<-cut(ph,seq(0,1,length.out=11),include.lowest=TRUE,labels=FALSE);o<-tapply(y,g,sum);e<-tapply(ph,g,sum);nn<-tapply(y,g,length);k<-!is.na(o);o<-o[k];e<-e[k];nn<-nn[k];st<-sum((o-e)^2/(e+1e-10)+((nn-o)-(nn-e))^2/((nn-e)+1e-10));if(length(o)>2)1-pchisq(st,length(o)-2) else NA_real_}
stuk<-function(fit){e<-predict(fit);d<-fit$data;d$za<-0.5*e^2*(e>=0);d$zb<- -0.5*e^2*(e<0);fa<-suppressWarnings(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial()));pchisq(deviance(fit)-deviance(fa),2,lower.tail=FALSE)}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)NA_real_);if(is.null(v)||length(v)!=1||!is.finite(v))NA_real_ else as.numeric(v)}
one_rep<-function(scn,n){G0<-gen(scn,n);dat<-G0$d;fit<-suppressWarnings(glm(G0$f,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  c(DEF.poly2=sc(ef_dir_cal(y,ph,fit,"poly",2)),DEF.poly3=sc(ef_dir_cal(y,ph,fit,"poly",3)),DEF.stk3=sc(ef_dir_cal(y,ph,fit,"stukel",3)),
    EF.omni=sc(ef_omni(y,ph)),HL=sc(hoslem<-ResourceSelection::hoslem.test(y,ph,g=10)$p.value),HLeqw=sc(hleqw(y,ph)),Stukel=sc(stuk(fit)))}
suppressMessages(library(ResourceSelection))
scns<-c("null",linkscns,"omit_quad","omit_int_bin","omit_int_cont"); ns<-c(200,500,1000,2000); REPS<-10000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L));clusterSetRNGStream(cl,20250911)
clusterEvalQ(cl,suppressMessages(library(ResourceSelection)))
clusterExport(cl,c("alpha","lg","inv_stukel","sqb","sib","scb","linkp","linkscns","gen","ef_dir_cal","ef_omni","hleqw","stuk","sc","one_rep"))
on.exit(stopCluster(cl)); tests<-c("DEF.poly2","DEF.poly3","DEF.stk3","EF.omni","HL","HLeqw","Stukel")
pall<-list(); summ<-list()
for(scn in scns) for(n in ns){clusterExport(cl,c("scn","n"),envir=environment())
  m<-parSapply(cl,1:REPS,function(i) one_rep(scn,n))               # 7 x REPS
  df<-as.data.frame(t(m)); names(df)<-tests; df<-cbind(scenario=scn,n=n,rep=1:REPS,df)
  pall[[length(pall)+1]]<-df
  pw<-rowMeans(m<alpha,na.rm=TRUE); summ[[length(summ)+1]]<-data.frame(scenario=scn,n=n,t(round(pw,4)))
  cat(sprintf("%-14s n=%-5d done\n",scn,n))}
write.csv(do.call(rbind,pall),"DEF_pvalues_all.csv",row.names=FALSE)
S<-do.call(rbind,summ); write.csv(S,"DEF_power_summary.csv",row.names=FALSE)
cat("\nSaved DEF_pvalues_all.csv (per-rep p-values) and DEF_power_summary.csv\n"); print(S,row.names=FALSE)
