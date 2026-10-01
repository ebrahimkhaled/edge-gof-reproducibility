## Writes headline_recount.csv — the per-cell enumeration behind fn:cells (19/22, 17/22).
## Rule: detectable = some compared test reaches size-adjusted power 0.15;
## non-saturated = not every compared test >= 0.97; tie margin 0.014.
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
SIM <- edge_path("code/simulations"); setwd(SIM)
pb <- read.csv("sim_power_broad.csv"); pb <- pb[pb$n==1000 & pb$alpha==0.05,]
PART<-c("EF","HL","HL-equalwidth","Pigeon-Heyse"); ALL6<-c(PART,"Tsiatis","Xie"); TIE<-0.014
g <- function(fam,par,t){v<-pb$power_size_adj[pb$family==fam&pb$param==par&pb$test==t]; if(length(v)) v[1] else NA}
cells <- unique(pb[pb$family!="null",c("family","param")])
rows <- lapply(seq_len(nrow(cells)), function(i){
  fam<-cells$family[i]; par<-cells$param[i]; p3<-g(fam,par,"DEF.poly3")
  rp<-sapply(PART,g,fam=fam,par=par); r6<-sapply(ALL6,g,fam=fam,par=par)
  cmp <- c(p3, rp)
  data.frame(source="sim_power_broad(size-adj)", scenario=paste(fam,par,sep="/"),
    edge_poly3=round(p3,3), best_partition=round(max(rp,na.rm=TRUE),3), best_partition_name=names(which.max(rp)),
    best_all6=round(max(r6,na.rm=TRUE),3), best_all6_name=names(which.max(r6)),
    detectable = max(c(p3,r6),na.rm=TRUE) >= 0.15, saturated = min(cmp,na.rm=TRUE) >= 0.97,
    beats_partition = p3 > max(rp,na.rm=TRUE)-TIE, beats_all6 = p3 > max(r6,na.rm=TRUE)-TIE)})
broad <- do.call(rbind, rows)
## tab_extra rows (raw power as printed in Table tab:extra; PR excluded from partition-rival set)
ex <- read.csv(text='scenario,edge_poly3,EF,HL,HLw,PH,Tsi,Xie
omit log x,0.957,0.842,0.847,0.860,0.790,0.817,0.803
omit x^2 (chi2_4 design),0.842,0.604,0.609,0.585,0.533,0.562,0.519
omit d1*d2,0.982,0.821,0.956,0.986,NA,NA,NA
omit x1*x2 (rho=.5),0.999,0.984,0.984,0.996,0.975,0.984,0.987
omit x^2 & x^3,0.954,0.805,0.792,0.768,0.721,0.769,0.759
omit x*d & x*z,0.878,0.846,0.842,0.906,0.805,1.000,1.000
omit d & x*d,0.182,0.112,0.116,0.111,0.083,0.111,0.096
omit z1 z2,0.052,0.057,0.057,0.045,0.034,0.048,0.038', stringsAsFactors=FALSE)
exr <- do.call(rbind, lapply(seq_len(nrow(ex)), function(i){ r<-ex[i,]
  rp<-unlist(r[c("EF","HL","HLw","PH")]); r6<-unlist(r[c("EF","HL","HLw","PH","Tsi","Xie")])
  cmp <- c(r$edge_poly3, rp)
  data.frame(source="tab_extra(raw)", scenario=r$scenario, edge_poly3=r$edge_poly3,
    best_partition=max(rp,na.rm=TRUE), best_partition_name=names(which.max(rp)),
    best_all6=max(r6,na.rm=TRUE), best_all6_name=names(which.max(r6)),
    detectable = max(c(r$edge_poly3,r6),na.rm=TRUE) >= 0.15, saturated = min(cmp,na.rm=TRUE) >= 0.97,
    beats_partition = r$edge_poly3 > max(rp,na.rm=TRUE)-TIE, beats_all6 = r$edge_poly3 > max(r6,na.rm=TRUE)-TIE)}))
out <- rbind(broad, exr)
out$in_headline_22 <- out$detectable & !out$saturated
write.csv(out, "headline_recount.csv", row.names=FALSE)
inh <- out[out$in_headline_22,]
cat(sprintf("cells in headline: %d | beats partition: %d | beats all6: %d\n", nrow(inh), sum(inh$beats_partition), sum(inh$beats_all6)))
print(inh[!inh$beats_partition, c("scenario","edge_poly3","best_partition","best_partition_name")], row.names=FALSE)
