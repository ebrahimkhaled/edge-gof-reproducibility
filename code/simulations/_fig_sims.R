suppressMessages({library(parallel); library(ResourceSelection)}); alpha<-0.05; lg<-function(p) log(p/(1-p))
inv_stukel<-function(ev,a1,a2){z<-numeric(length(ev));pos<-ev>=0
  if(any(pos)) z[pos]<- if(abs(a1)<1e-12) ev[pos] else (-1+sqrt(pmax(0,1+2*a1*ev[pos])))/a1
  if(any(!pos))z[!pos]<-if(abs(a2)<1e-12) ev[!pos] else (-1+sqrt(pmax(0,1+2*a2*ev[!pos])))/a2; plogis(z)}
sqb<-function(J) solve(rbind(c(1,-1.5,2.25),c(1,3,9),c(1,-3,9)),c(lg(.05),lg(.95),lg(J)))
ef_dir_cal<-function(y,ph,fit,basis,k){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);G<-10
  grp<-pmin(ceiling(rank(ph,ties.method="first")/(n/G)),G);idx<-split(seq_len(n),grp);Wii<-ph*(1-ph)
  og<-sapply(idx,function(I)sum(y[I]));eg<-sapply(idx,function(I)sum(ph[I]));Vg<-sapply(idx,function(I)sum(Wii[I]));pbar<-sapply(idx,function(I)mean(ph[I]));r<-(og-eg)/sqrt(Vg)
  X<-model.matrix(fit);U<-t(sapply(idx,function(I)colSums(Wii[I]*X[I,,drop=FALSE])));U<-U/sqrt(Vg);Om<-diag(G)-U%*%solve(crossprod(X,Wii*X))%*%t(U)
  if(basis=="poly")Z<-as.matrix(poly(pbar,k)) else{e<-qlogis(pbar);Z<-cbind(e,e^2*(e>=0),-e^2*(e<0));Z<-Z[,colSums(abs(Z))>1e-8,drop=FALSE]}
  if(ncol(Z)<1)return(NA);Zi<-solve(crossprod(Z));Ztr<-crossprod(Z,r);S<-as.numeric(t(Ztr)%*%Zi%*%Ztr)
  lam<-Re(eigen(Zi%*%(t(Z)%*%Om%*%Z),only.values=TRUE)$values);lam<-lam[lam>1e-9];if(!length(lam))return(NA)
  cc<-sum(lam^2)/sum(lam);nu<-sum(lam)^2/sum(lam^2);1-pchisq(S/cc,nu)}
ef_omni<-function(y,ph){ph<-pmin(pmax(ph,1e-6),1-1e-6);n<-length(y);g<-pmin(ceiling(rank(ph,ties.method="first")/(n/10)),10)
  o<-tapply(y,g,sum);e<-tapply(ph,g,sum);ng<-tapply(y,g,length);pb<-as.numeric(tapply(ph,g,mean));V<-ng*pb*(1-pb);oe<-as.numeric(o-e);1-pchisq(sum(oe^2/V)-sum((1-2*pb)*oe/V),8)}
stuk<-function(fit){e<-predict(fit);d<-fit$data;d$za<-0.5*e^2*(e>=0);d$zb<- -0.5*e^2*(e<0);fa<-suppressWarnings(glm(update(formula(fit),.~.+za+zb),data=d,family=binomial()));pchisq(deviance(fit)-deviance(fa),2,lower.tail=FALSE)}
sc<-function(e){v<-tryCatch(e,error=function(x)NA);if(is.null(v)||length(v)!=1||!is.finite(v))NA else as.numeric(v)}
pv5<-function(dat,frm){fit<-suppressWarnings(glm(frm,data=dat,family=binomial()));ph<-as.numeric(fitted(fit));y<-dat$y
  c(DEF.poly3=sc(ef_dir_cal(y,ph,fit,"poly",3)),DEF.stk3=sc(ef_dir_cal(y,ph,fit,"stukel",3)),EF.omni=sc(ef_omni(y,ph)),HL=sc(hoslem.test(y,ph,g=10)$p.value),Stukel=sc(stuk(fit)))}
genQ<-function(J,n){b<-sqb(J);x<-runif(n,-3,3);list(d=data.frame(x=x,y=rbinom(n,1,plogis(b[1]+b[2]*x+b[3]*x^2))),f=y~x)}
genA<-function(a,n){x<-runif(n,-3,3);p<-if(a==0)plogis(0.8*x) else inv_stukel(0.8*x,a,a);list(d=data.frame(x=x,y=rbinom(n,1,p)),f=y~x)}
genW<-function(w,n){x<-runif(n,-3,3);p<-plogis(0.8*x+if(w==0)0 else 1.5*sin(w*x));list(d=data.frame(x=x,y=rbinom(n,1,p)),f=y~x)}
n<-1000;REPS<-3000;cl<-makeCluster(max(1L,detectCores(logical=TRUE)-1L));clusterSetRNGStream(cl,20250911)
clusterEvalQ(cl,suppressMessages(library(ResourceSelection)));clusterExport(cl,c("alpha","lg","inv_stukel","sqb","ef_dir_cal","ef_omni","stuk","sc","pv5","genQ","genA","genW","n"))
on.exit(stopCluster(cl)); tests<-c("DEF.poly3","DEF.stk3","EF.omni","HL","Stukel")
run<-function(genname,grid){do.call(rbind,lapply(grid,function(v){clusterExport(cl,c("genname","v"),envir=environment())
  m<-parSapply(cl,1:REPS,function(i){g<-get(genname)(v,n);pv5(g$d,g$f)});data.frame(sev=v,t(rowMeans(m<alpha,na.rm=TRUE)))}))}
Q<-run("genQ",c(0.01,0.02,0.03,0.05,0.10,0.20)); Q$family<-"Omitted quadratic (J)"; names(Q)[2:6]<-tests
A<-run("genA",c(0,0.3,0.6,1.0,1.5,2.0)); A$family<-"Stukel link (alpha)"; names(A)[2:6]<-tests
W<-run("genW",c(0,0.5,1,1.5,2,3,4,6,8)); names(W)[2:6]<-tests
write.csv(rbind(Q,A),"Fig_severity.csv",row.names=FALSE); write.csv(W,"Fig_omega.csv",row.names=FALSE)
cat("Saved Fig_severity.csv and Fig_omega.csv\n"); print(rbind(Q,A),row.names=FALSE); print(W,row.names=FALSE)
