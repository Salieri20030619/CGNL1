args<-commandArgs(trailingOnly=TRUE);p<-args[1]
d<-read.csv(gzfile(file.path(p,'FEATURES.csv.gz')),check.names=FALSE)
stopifnot(nrow(d)>=1000,!any(grepl('Q0VF96',d$feature,fixed=TRUE)))
a<-log2(as.matrix(d[,3:18]));b<-log2(as.matrix(d[,19:23]))
rp<-expand.grid(i=1:16,j=1:16);rp<-rp[rp$i!=rp$j,]
pp<-t(combn(1:5,2));out<-list()
for(s in 0:1){
 m<-d$split==s;stopifnot(sum(m)>=500)
 rd<-vapply(seq_len(nrow(rp)),function(k)a[m,rp$i[k]]-a[m,rp$j[k]],numeric(sum(m)))
 pd<-vapply(seq_len(nrow(pp)),function(k)b[m,pp[k,1]]-b[m,pp[k,2]],numeric(sum(m)))
 cc<-cor(rd,pd,method='spearman')
 for(k in seq_len(nrow(pp))){
  o<-order(cc[,k],decreasing=TRUE);ii<-o[1]
  out[[length(out)+1]]<-data.frame(split=s,sample_i=pp[k,1],sample_j=pp[k,2],channel_i=rp$i[ii],channel_j=rp$j[ii],rho=cc[ii,k],runnerup=cc[o[2],k],margin=cc[ii,k]-cc[o[2],k])
 }
}
write.csv(do.call(rbind,out),file.path(p,'CHANNEL_INDEPENDENT_R.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(p,'CHANNEL_INDEPENDENT_R_SESSION.txt'))

map<-c(10,12,13,14,15)
sc<-do.call(rbind,lapply(2:5,function(j){v<-(b[,j]-b[,1])-(a[,map[j]]-a[,map[1]]);q<-median(v[d$split==0]);h<-v[d$split==1]-q;data.frame(sample=names(d)[18+j],offset=q,heldout_median=median(h),heldout_abs_p95=as.numeric(quantile(abs(h),.95)))}))
write.csv(sc,file.path(p,"SCALING_R.csv"),row.names=FALSE)
