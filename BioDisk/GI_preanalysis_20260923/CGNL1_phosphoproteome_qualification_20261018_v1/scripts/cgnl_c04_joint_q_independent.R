args<-commandArgs(trailingOnly=TRUE);w<-args[1]
f<-file.path(w,"PSM_QVALUES.tsv.gz")
if(!file.exists(f))stop("Joint input is not ready; no q values computed")
d<-read.delim(gzfile(f),check.names=FALSE,quote="",stringsAsFactors=FALSE)
decoy<-tolower(as.character(d$is_decoy))=="true"
phospho<-tolower(as.character(d$is_phospho))=="true"
calcq<-function(score,dd){
 if(!length(score))return(numeric())
 o<-order(score);ss<-score[o];ends<-c(which(diff(ss)!=0),length(ss))
 nd<-cumsum(as.integer(dd[o]))[ends];nt<-ends-nd
 q<-rev(cummin(rev(pmin(1,(nd+1)/pmax(nt,1)))))
 v<-rep(q,diff(c(0,ends)));ans<-numeric(length(v));ans[o]<-v;ans
}
qa<-rep(NA_real_,nrow(d));qp<-qa
for(fam in unique(d$family)){
 ii<-which(d$family==fam);qa[ii]<-calcq(d$e_value[ii],decoy[ii])
 jj<-ii[phospho[ii]];qp[jj]<-calcq(d$e_value[jj],decoy[jj])
}
err1<-max(abs(qa-d$q_global));err2<-if(any(phospho))max(abs(qp[phospho]-d$q_phospho[phospho])) else 0
stopifnot(is.finite(err1),is.finite(err2),err1<1e-12,err2<1e-12)
write.csv(data.frame(rows=nrow(d),families=length(unique(d$family)),global_max_error=err1,phospho_max_error=err2),file.path(w,"INDEPENDENT_R_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
writeLines(capture.output(sessionInfo()),file.path(w,"INDEPENDENT_R_SESSION.txt"),useBytes=TRUE)
