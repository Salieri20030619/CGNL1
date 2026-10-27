args<-commandArgs(trailingOnly=TRUE);p<-args[1]
d<-read.csv(file.path(p,'SET04_DESCRIPTIVE_INPUT.csv'))
r<-log2(d$S486)-log2(d$S283);dd<-r[1]-r[-1]
write.csv(data.frame(normal=d$sample[-1],delta=dd),file.path(p,'SET04_DESCRIPTIVE_R.csv'),row.names=FALSE)
stopifnot(max(abs(r-d$log2ratio))<1e-12)
