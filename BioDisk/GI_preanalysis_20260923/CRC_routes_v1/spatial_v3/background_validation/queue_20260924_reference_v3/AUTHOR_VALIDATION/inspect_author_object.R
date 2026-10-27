args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(Matrix))
x <- readRDS(file.path(args[1], "merged_NDE_22261_fisher.rds"))
sink(file.path(args[1],"OBJECT_STRUCTURE.txt"))
cat("class:", class(x), "\n")
cat("S4:", isS4(x), "\n")
if(isS4(x)) {
 cat("attributes:",names(attributes(x)),"\n")
 for(n in setdiff(names(attributes(x)), "class")){
  y <- attr(x,n);cat("\nATTRIBUTE",n,"class",class(y),"dim",dim(y),"length",length(y),"\n")
  if(is.list(y)){cat("names:",head(names(y),30),"\n");str(y,max.level=1,give.attr=FALSE)}
  if(!is.null(rownames(y)))cat("rowhead:",head(rownames(y),8),"\n")
  if(!is.null(colnames(y)))cat("colhead:",head(colnames(y),8),"\n")
 }
} else str(x,max.level=2,give.attr=FALSE)
sink()
cat("Object structure recorded\n")
