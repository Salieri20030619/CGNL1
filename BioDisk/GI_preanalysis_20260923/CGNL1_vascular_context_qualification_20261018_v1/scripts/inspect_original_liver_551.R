suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly=TRUE)
x <- readRDS(args[1])
info <- list(class=class(x), slots=slotNames(x))
if (inherits(x,"Seurat")) {
 meta <- x@meta.data
 info$dimensions <- dim(x)
 info$metadata_columns <- colnames(meta)
 info$metadata_categories <- lapply(meta, function(v) {v<-as.character(v); if(length(unique(v))<=80) as.list(table(v,useNA="ifany")) else list(unique_count=length(unique(v)))})
 info$assays <- names(x@assays)
 write.table(cbind(cell=rownames(meta),meta), file=args[3], sep="\t", row.names=FALSE, quote=TRUE, fileEncoding="UTF-8")
}
writeLines(toJSON(info,pretty=TRUE,auto_unbox=TRUE,null="null"),args[2],useBytes=TRUE)
