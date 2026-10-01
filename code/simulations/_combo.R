suppressMessages(library(parallel)); alpha<-0.05; lg<-function(p) log(p/(1-p))
inv_stukel<-function(ev,a1,a2){z<-numeric(length(ev));pos<-ev>=0
  if(any(pos)) z[pos]<- if(abs(a1)<1e-12) ev[pos] else (-1+sqrt(pmax(0,1+2*a1*ev[pos])))/a1
  if(any(!pos))z[!pos]<-if(abs(a2)<1e-12) ev[!pos] else (-1+sqrt(pmax(0,1+2*a2*ev[!pos])))/a2; plogis(z)}
sqb<-function(J) solve(rbind(c(1,-1.5,2.25),c(1,3,9),c(1,-3,9)),c(lg(.05),lg(.95),lg(J)))
linkp<-function(scn,eta) switch(scn, cloglog=1-exp(-exp(eta)), loglog=exp(-exp(-eta)), cauchit=pcauchy(eta),
  scobit2=plogis(eta)^2, stukel_heavy=inv_stukel(eta,-1,-1), stukel_asym=inv_stukel(eta,-1,1))
linkscns<-c("cloglog","loglog","cauchit","scobit2","stukel_heavy","stukel_asym")
gen<-function(scn,n){
  if(scn=="null"){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,plogis(0.6*x+0.5*d))),f=y~x+d)}
  else if(scn%in%linkscns){x<-runif(n,-3,3);d<-rbinom(n,1,.5);list(d=data.frame(x=x,d=d,y=rbinom(n,1,linkp(scn,0.6*x+0.5*d))),f=y~x+d)}
  else if(scn=="quad"){b<-sqb(0.02);x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*x^2))),f=y~x)}
  else {x<-runif(n,-2.5,2.5);list(d=data.frame(x=x,y=rbinom(n,1,plogis(0.3+0.7*x+0.12*(x^3-3*x)))),f=y~x)}}
ef_dir_cal<-function(y,ph,fit,basis){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])))/sqrt(Vg);Om<-diag(G)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  e<-qlogis(pbar)
  if(basis=="poly2")Z<-as.matrix(poly(pbar,2))
  else if(basis=="poly3")Z<-as.matrix(poly(pbar,3))
  else if(basis=="stukel"){Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  else {Zp<-as.matrix(poly(pbar,3));Zs<-cbind(e,e^2*(e>=0),-e^2*(e<0));Zs<-Zs[,colSums(abs(Zs))>1e-8,drop=FALSE];Zc<-cbind(Zp,Zs);q0<-qr(Zc,tol=1e-7);Z<-qr.Q(q0)[,seq_len(q0$rank),drop=FALSE]}
  Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(c(p=NA_real_,k=ncol(Z)))
  cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);c(p=1-pchisq(S/cc,nu),k=ncol(Z))}
sc<-function(e){v<-tryCatch(suppressWarnings(e),error=function(x)c(p=NA_real_,k=NA));if(is.null(v))c(p=NA_real_,k=NA) else v}
one_rep<-function(scn,n){G0<-gen(scn,n);dat<-G0$d;fit<-suppressWarnings(glm(G0$f,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  cb<-sc(ef_dir_cal(y,ph,fit,"combo"))
  c(DEF.poly2=sc(ef_dir_cal(y,ph,fit,"poly2"))["p"],DEF.poly3=sc(ef_dir_cal(y,ph,fit,"poly3"))["p"],DEF.stk3=sc(ef_dir_cal(y,ph,fit,"stukel"))["p"],DEF.combo=cb["p"],combo.k=cb["k"])}
scns<-c("null","cloglog","loglog","cauchit","scobit2","stukel_heavy","stukel_asym","quad","cubic");n<-1000;REPS<-4000
cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L));clusterSetRNGStream(cl,30303)
clusterExport(cl,c("alpha","lg","inv_stukel","sqb","linkp","linkscns","gen","ef_dir_cal","sc","one_rep","n"));on.exit(stopCluster(cl))
nm<-c("DEF.poly2.p","DEF.poly3.p","DEF.stk3.p","DEF.combo.p");pall<-list();summ<-list()
for(scn in scns){clusterExport(cl,"scn",envir=environment());m<-parSapply(cl,1:REPS,function(i)one_rep(scn,n))
  df<-as.data.frame(t(m[c(nm,"combo.k.k"),]));names(df)<-c("DEF.poly2","DEF.poly3","DEF.stk3","DEF.combo","combo.k");df<-cbind(scenario=scn,rep=1:REPS,df);pall[[length(pall)+1]]<-df
  pw<-rowMeans(m[nm,,drop=FALSE]<alpha,na.rm=TRUE);summ[[length(summ)+1]]<-data.frame(scenario=scn,combo.k=round(mean(m["combo.k.k",]),2),DEF.poly2=round(pw[1],3),DEF.poly3=round(pw[2],3),DEF.stk3=round(pw[3],3),DEF.combo=round(pw[4],3));cat(scn,"done\n")}
write.csv(do.call(rbind,pall),"Combo_pvalues.csv",row.names=FALSE)
S<-do.call(rbind,summ);write.csv(S,"Combo_power.csv",row.names=FALSE);print(S,row.names=FALSE)
