args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(Matrix))
D <- args[1];x <- readRDS(file.path(D,"merged_NDE_22261_fisher.rds"))
cts <- attr(x,"counts");ge <- attr(x,"gene_names");ba <- attr(x,"batch_ID");nc <- attr(x,"num_cells");co <- attr(x,"coord")
stopifnot(nrow(cts)==length(ba),ncol(cts)==length(ge),identical(colnames(cts),ge),all(cts@x>=0),all(abs(cts@x-round(cts@x))<1e-8))
writeMM(t(cts),file.path(D,"author_counts.mtx"))
write.table(data.frame(gene=ge),file.path(D,"author_genes.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
tab <- data.frame(batch=ba,barcode=attr(x,"cell_names"),imagerow=co[,1],imagecol=co[,2],library_provided_genes=Matrix::rowSums(cts),genes_detected=Matrix::rowSums(cts>0),as.data.frame(nc),check.names=FALSE)
write.table(tab,file.path(D,"author_spots.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
targets <- intersect(c("CGNL1","CREB3L1","PCSK5","SLC9B2"),ge)
rows <- list()
for(direction in c("pos","neg")){
 obj <- attr(x,paste0("niche_DE_pval_",direction));p <- obj$gene_level;p[!is.finite(p)] <- 1;q <- p.adjust(p,"BH")
 for(g in targets){i<-match(g,ge);rows[[length(rows)+1]]<-data.frame(gene=g,direction=direction,stored_gene_level_P=obj$gene_level[i],recomputed_all_gene_BH=q[i],family_size=length(p))}
}
write.table(do.call(rbind,rows),file.path(D,"author_gene_level_novelty_audit.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
sink(file.path(D,"EXPORT_NOTES.txt"));print(table(ba));cat("scale:",attr(x,"scale"),"spot_distance:",attr(x,"spot_distance"),"sigma:",attr(x,"sigma"),"\n");print(summary(rowSums(nc)))
for(g in targets){
 cat("\nGENE",g,"\n")
 for(s in names(attr(x,"niche_DE"))){cat("SIGMA",s,"\n");str(attr(x,"niche_DE")[[s]][[match(g,ge)]],max.level=2)}
}
sink();writeLines('{"exit_code":0,"integer_counts_verified":true}',file.path(D,"EXPORT_COMPLETE.json"))
cat("Author raw matrix and metadata exported\n")
