args<-commandArgs(trailingOnly=TRUE);p<-args[1]
d<-read.csv(gzfile(file.path(p,'CALIBRATION_INPUT.csv.gz')),check.names=FALSE)
t<-read.csv(file.path(p,'TARGET_INPUT.csv'),check.names=FALSE)
stopifnot(!any(grepl('Q0VF96',d$feature,fixed=TRUE)))
keys<-unique(d[,c('set','layer','sample')]);factors<-list()
for(i in seq_len(nrow(keys))){
 k<-keys[i,];x<-d[d$set==k$set & d$layer==k$layer & d$sample==k$sample,]
 train<-x[x$split==0,];held<-x[x$split==1,];stopifnot(nrow(train)>=500,nrow(held)>=500)
 b<-median(train$export_relative_log-train$raw_relative_log)
 factors[[i]]<-data.frame(set=k$set,layer=k$layer,sample=k$sample,A=0,B=b,C=-median(train$raw_relative_log),heldout_median=median(held$export_relative_log-held$raw_relative_log-b))
}
f<-do.call(rbind,factors);write.csv(f,file.path(p,'FACTORS_R.csv'),row.names=FALSE)
res<-list()
for(s in unique(t$set))for(m in c('A','B','C')){
 x<-t[t$set==s,];v<-log2(x$S486)-log2(x$CGNL1_protein)
 for(i in seq_len(nrow(x))){
  pp<-f[f$set==s & f$layer=='phospho' & f$sample==x$sample[i],]
  gp<-f[f$set==s & f$layer=='protein' & f$sample==x$sample[i],]
  stopifnot(nrow(pp)==1,nrow(gp)==1)
  if(m=='B' && x$group[i] %in% c('tumor','normal'))stopifnot(abs(pp$heldout_median)<=.01,abs(gp$heldout_median)<=.01)
  v[i]<-v[i]+pp[[m]]-gp[[m]]
 }
 for(i in which(x$group=='tumor'))for(j in which(x$group=='normal')){
  res[[length(res)+1]]<-data.frame(set=s,model=m,tumor=x$sample[i],normal=x$sample[j],delta=v[i]-v[j])
 }
}
write.csv(do.call(rbind,res),file.path(p,'COMPARISONS_R.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(p,'R_SESSION.txt'))
