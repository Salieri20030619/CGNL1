args<-commandArgs(trailingOnly=TRUE);p<-args[1]
d<-read.csv(gzfile(file.path(p,"NORMAL_MAPPING_FEATURES.csv.gz")),check.names=FALSE)
samples<-c("pCCA13","iCCA32","dCCA2","Normal7","Normal8","Normal9");proc<-c("processed2","processed3","processed4","normal7","normal8","normal9");raw<-c(2,3,10,13,14,15)
g<-read.csv(file.path(p,"Set04_protein_mapping/SCALING_R.csv"))
out<-do.call(rbind,lapply(1:6,function(i){v<-d[[proc[i]]]-d$processed1-(d[[paste0("raw",raw[i])]]-d$raw1);fit<-median(v[d$split==0],na.rm=TRUE);test<-median(v[d$split==1]-fit,na.rm=TRUE);gp<-g$offset[g$sample==samples[i]];data.frame(sample=samples[i],phospho=fit,heldout=test,protein=gp,difference=fit-gp)}))
write.csv(out,file.path(p,"CROSS_LAYER_SCALING_R.csv"),row.names=FALSE)
