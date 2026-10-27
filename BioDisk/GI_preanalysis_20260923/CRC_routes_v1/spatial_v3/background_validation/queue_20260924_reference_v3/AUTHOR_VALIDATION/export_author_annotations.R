a <- commandArgs(trailingOnly=TRUE);D<-a[1];suppressPackageStartupMessages(library(Matrix));x<-readRDS(file.path(D,"merged_NDE_22261_fisher.rds"));ge<-attr(x,"gene_names");ct<-attr(x,"cell_types")
write.table(data.frame(batch=attr(x,"batch_ID"),barcode=rownames(attr(x,"num_cells"))),file.path(D,"author_barcodes.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
rows<-list()
for(g in c("CGNL1","CREB3L1")){
 gi<-match(g,ge)
 for(s in names(attr(x,"niche_DE"))){
  y<-attr(x,"niche_DE")[[s]][[gi]];b<-y$betas;t<-y$T_stat
  if(!is.matrix(b))next
  for(i in seq_along(ct))for(j in seq_along(ct)){
   rows[[length(rows)+1]]<-data.frame(gene=g,sigma=s,index=ct[i],niche=ct[j],beta=b[i,j],T_stat=t[i,j],stored_P_pos=attr(x,"niche_DE_pval_pos")$interaction_level[i,j,gi],stored_P_neg=attr(x,"niche_DE_pval_neg")$interaction_level[i,j,gi])
  }
 }
}
write.table(do.call(rbind,rows),file.path(D,"author_target_interactions.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
