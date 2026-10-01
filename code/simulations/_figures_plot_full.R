suppressMessages({library(ggplot2); library(tidyr); library(dplyr)})
ord<-c("DEF.poly2","DEF.poly3","DEF.stk3","EF.omni","HL","HLeqw","PH","Tsiatis","Xie","Stukel")
lab<-c(DEF.poly2="DEF-poly2",DEF.poly3="DEF-poly3",DEF.stk3="DEF-stk",EF.omni="EF",HL="HL",HLeqw="HL-w",PH="Pigeon-Heyse",Tsiatis="Tsiatis",Xie="Xie",Stukel="Stukel")
cols<-c("DEF-poly2"="#D55E00","DEF-poly3"="#B22222","DEF-stk"="#CC79A7","EF"="#0072B2","HL"="#56B4E9","HL-w"="#999999","Pigeon-Heyse"="#444444","Tsiatis"="#009E73","Xie"="#66C2A5","Stukel"="#000000")
fam<-function(t) ifelse(grepl("DEF",t),"DEF",ifelse(t=="Stukel","Stukel","Other partition test"))
lwd<-c("DEF"=1.2,"Stukel"=1.0,"Other partition test"=0.55)
prep<-function(df) {d<-pivot_longer(df,all_of(ord),names_to="t",values_to="Power");d$Test<-factor(lab[d$t],lab[ord]);d$Fam<-fam(d$Test);d}

sev<-prep(read.csv("Fig_severity_full.csv"))
f2<-ggplot(sev,aes(sev,Power,color=Test,linewidth=Fam,group=Test))+geom_line()+geom_point(size=1.1,show.legend=FALSE)+
  facet_wrap(~family,scales="free_x")+scale_color_manual(values=cols)+scale_linewidth_manual(values=lwd,guide="none")+
  scale_y_continuous(limits=c(0,1))+labs(x="Misspecification severity",y="Power",color=NULL,
    title="Power vs. severity: DEF (thick) tracks Stukel and dominates all other partition tests")+
  guides(color=guide_legend(nrow=2,override.aes=list(linewidth=1.1)))+
  theme_minimal(base_size=10.5)+theme(legend.position="bottom",plot.title=element_text(face="bold",size=9.5),strip.text=element_text(face="bold"))
ggsave("Fig_severity.png",f2,width=8,height=4.4,dpi=200)

om<-prep(read.csv("Fig_omega_full.csv"))
f3<-ggplot(om,aes(sev,Power,color=Test,linewidth=Fam,group=Test))+
  annotate("rect",xmin=2.5,xmax=8.5,ymin=0,ymax=1,alpha=.05,fill="red")+geom_line()+geom_point(size=1.1,show.legend=FALSE)+
  scale_color_manual(values=cols)+scale_linewidth_manual(values=lwd,guide="none")+scale_y_continuous(limits=c(0,1))+
  annotate("text",x=5.6,y=0.42,label="high-frequency regime:\nomnibus & covariate-space\ntests stay powerful;\ndirected tests (DEF, Stukel) fail",size=2.8,color="grey25")+
  labs(x=expression(paste("Oscillation frequency  ",omega,"     ( logit = 0.8x + 1.5 sin(",omega,"x) )")),y="Power",color=NULL,
    title="Boundary: under rough, high-frequency misfit the directed tests (DEF, Stukel) lose to the omnibus tests")+
  guides(color=guide_legend(nrow=2,override.aes=list(linewidth=1.1)))+
  theme_minimal(base_size=10.5)+theme(legend.position="bottom",plot.title=element_text(face="bold",size=9),legend.text=element_text(size=8.5))
ggsave("Fig_omega.png",f3,width=8,height=4.6,dpi=200)
cat("Rebuilt Fig_severity.png and Fig_omega.png with full partition set\n")
