o<-commandArgs(trailingOnly=TRUE)[1];d<-read.csv(file.path(o,"inputs.csv"),stringsAsFactors=FALSE)
out<-lapply(unique(d$proxy),function(k){z<-d[d$proxy==k,];X<-cbind(1,z$purity*10);y<-z$value;b<-solve(crossprod(X),crossprod(X,y));e<-as.vector(y-X%*%b);bread<-solve(crossprod(X));h<-rowSums((X%*%bread)*X);V<-bread%*%crossprod(X,X*(e/(1-h))^2)%*%bread;se<-sqrt(V[2,2]);data.frame(proxy=k,n=nrow(z),beta=b[2],SE=se,P=2*pt(-abs(b[2]/se),df=nrow(z)-2),rho=cor(z$purity,z$value,method="spearman"))})
write.csv(do.call(rbind,out),file.path(o,"R_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
