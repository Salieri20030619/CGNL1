args <- commandArgs(trailingOnly=TRUE)
o <- args[1]
d <- read.csv(file.path(o,"MODEL_INPUT.csv"),stringsAsFactors=FALSE)
out <- list()
for(sid in c("S0","S1","S2","S3")) {
 z <- if(sid=="S1") d[d$repeated==0,] else d
 f <- y ~ site + protein + age10 + male + size_cm
 if(sid %in% c("S2","S3")) {
  z$batch <- factor(if(sid=="S2") z$batch_original else z$batch_22,levels=sort(unique(d$batch_original)))
  f <- update(f,. ~ . + batch)
 }
 fit <- glm(f,data=z,family=binomial(),control=glm.control(epsilon=1e-12,maxit=200))
 b <- coef(summary(fit))["site",]
 out[[sid]] <- data.frame(scenario=sid,n=nrow(z),beta=b[1],se=b[2],P=b[4],converged=fit$converged,rank=fit$rank)
}
write.csv(do.call(rbind,out),file.path(o,"R_INDEPENDENT_QA.csv"),row.names=FALSE,fileEncoding="UTF-8")
writeLines(capture.output(sessionInfo()),file.path(o,"R_SESSION_INFO.txt"),useBytes=TRUE)
