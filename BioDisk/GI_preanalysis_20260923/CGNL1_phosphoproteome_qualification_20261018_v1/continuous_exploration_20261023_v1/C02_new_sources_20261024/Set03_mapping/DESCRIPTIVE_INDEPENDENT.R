args<-commandArgs(trailingOnly=TRUE);p<-args[1]
d<-read.csv(file.path(p,'DESCRIPTIVE_INPUT.csv'));r<-log2(d$S486)-log2(d$S283)
stopifnot(max(abs(r-d$log2ratio))<1e-12)
out<-do.call(rbind,lapply(1:3,function(i)data.frame(tumor=d$sample[i],normal=d$sample[4:5],delta=r[i]-r[4:5])))
write.csv(out,file.path(p,'DESCRIPTIVE_R.csv'),row.names=FALSE)
