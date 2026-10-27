library(jsonlite)
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
w<-args[1];x<-read.delim(gzfile(file.path(w,'PSM_QVALUES.tsv.gz')),check.names=FALSE,stringsAsFactors=FALSE)
o<-order(x$e_value);d<-as.logical(x$is_decoy[o]);score<-x$e_value[o]
ends<-which(!duplicated(score,fromLast=TRUE));dc<-cumsum(d)[ends];tc<-cumsum(!d)[ends]
fdr<-pmin(1,(dc+1)/pmax(tc,1));q<-rev(cummin(rev(fdr)));calculated<-rep(q,diff(c(0,ends)))
error<-max(abs(calculated-x$q_pilot[o]));stopifnot(error<1e-12)
pass<-!as.logical(x$is_decoy)&!as.logical(x$is_contaminant)&x$q_pilot<=.01
r<-fromJSON(file.path(w,'PSM_QA_RESULTS.json'));stopifnot(nrow(x)==r$PSMs,sum(pass)==r$noncontaminant_target_PSM_q01,length(unique(x$peptide[pass]))==r$unique_passing_peptides)
qa<-list(pass=TRUE,rows=nrow(x),max_q_error=error,passing_noncontaminant_PSMs=sum(pass),scope='single fraction pilot, not full cohort or protein FDR')
write_json(qa,file.path(w,'INDEPENDENT_R_PSM_QA.json'),auto_unbox=TRUE,pretty=TRUE,digits=17);cat(toJSON(qa,auto_unbox=TRUE,digits=17),'\n')
