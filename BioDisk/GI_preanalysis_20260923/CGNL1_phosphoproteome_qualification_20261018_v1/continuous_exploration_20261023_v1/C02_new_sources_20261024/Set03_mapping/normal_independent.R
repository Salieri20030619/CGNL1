args<-commandArgs(trailingOnly=TRUE);p<-args[1]
d<-read.csv(gzfile(file.path(p,'NORMAL_MAPPING_FEATURES.csv.gz')),check.names=FALSE)
a<-as.matrix(d[,3:18]);b<-as.matrix(d[,19:21]);normal<-as.matrix(d[,22:30]);anchors<-c(10,12,13);out<-list()
for(n in 1:9)for(k in 1:3)for(s in 0:1){
 m<-d$split==s & is.finite(normal[,n]);channels<-setdiff(1:16,anchors[k])
 rd<-vapply(channels,function(c)a[m,c]-a[m,anchors[k]],numeric(sum(m)))
 cc<-as.numeric(cor(rd,normal[m,n]-b[m,k],method='spearman'));o<-order(cc,decreasing=TRUE)
 out[[length(out)+1]]<-data.frame(normal=paste0('Normal',n),anchor_channel=anchors[k],split=s,best_channel=channels[o[1]],rho=cc[o[1]],margin=cc[o[1]]-cc[o[2]])
}
write.csv(do.call(rbind,out),file.path(p,'NORMAL_MAPPING_INDEPENDENT_R.csv'),row.names=FALSE)
