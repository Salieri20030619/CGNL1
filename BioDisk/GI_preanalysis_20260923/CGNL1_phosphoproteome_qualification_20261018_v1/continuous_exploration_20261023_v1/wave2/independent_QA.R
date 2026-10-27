o<-commandArgs(trailingOnly=TRUE)[1];d<-read.csv(file.path(o,"patient_contrasts.csv"),stringsAsFactors=FALSE)
out<-lapply(unique(d$test),function(k){x<-d$value[d$test==k];n<-length(x);p<-if(n>=10) wilcox.test(x,exact=TRUE,correct=FALSE)$p.value else NA_real_;data.frame(test=k,n=n,mean=mean(x),median=median(x),P=p)})
write.csv(do.call(rbind,out),file.path(o,"R_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
