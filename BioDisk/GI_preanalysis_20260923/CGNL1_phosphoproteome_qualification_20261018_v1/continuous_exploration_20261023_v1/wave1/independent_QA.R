args<-commandArgs(trailingOnly=TRUE);o<-args[1];d<-read.csv(file.path(o,"patient_contrasts.csv"),stringsAsFactors=FALSE)
exactp<-function(x) {
 x<-x[x!=0];if(!length(x))return(NA_real_)
 if(!anyDuplicated(abs(x)))return(wilcox.test(x,mu=0,alternative="two.sided",exact=TRUE,correct=FALSE)$p.value)
 ranks<-as.integer(round(2*rank(abs(x))));obs<-sum(ranks[x>0]);cut<-min(obs,sum(ranks)-obs)
 p<-c(1,rep(0,cut))
 for(r in ranks){old<-p;p<-old/2;if(r<=cut)p[(r+1):(cut+1)]<-p[(r+1):(cut+1)]+old[1:(cut-r+1)]/2}
 min(1,2*sum(p))
}
out<-lapply(unique(d$test),function(k){z<-d[d$test==k,];data.frame(test=k,n=nrow(z),mean=mean(z$value),median=median(z$value),P=exactp(z$value),alternative_P=exactp(z$alternative))})
write.csv(do.call(rbind,out),file.path(o,"independent_R_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
writeLines(capture.output(sessionInfo()),file.path(o,"R_session.txt"),useBytes=TRUE)
