suppressPackageStartupMessages(library(edgeR))
r <- '/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923'
out <- file.path(r,'04_results','CRC');dir.create(out,showWarnings=FALSE,recursive=TRUE)
x <- read.delim(gzfile(file.path(r,'02_raw','GSE145356_counts.txt.gz')),check.names=FALSE)
dup <- duplicated(x$Ensembl.Gene.ID); write.table(x[dup,1:3],file.path(out,'duplicate_annotation_rows.tsv'),sep='\t',row.names=FALSE,quote=FALSE); x <- x[!dup,]
 cnt <- as.matrix(x[,-(1:3)]);storage.mode(cnt)<-'numeric';rownames(cnt)<-x$Ensembl.Gene.ID
stopifnot(!anyNA(cnt),all(cnt>=0),all(cnt==round(cnt)),!anyDuplicated(rownames(cnt)))
g <- rep(NA_character_,ncol(cnt));names(g)<-colnames(cnt)
g[grepl('BRCD_T1-',names(g))]<-'Saline';g[grepl('BRCD_T2-',names(g))]<-'DMSO_Barcode';g[grepl('BRCD_T3-',names(g))]<-'5FU_LV';g[grepl('BRCD_T5-',names(g))]<-'CPT_Regrowth';g[grepl('BRCD_T7-',names(g))]<-'FOLFIRI_Regrowth';g[grepl('_TR_',names(g))]<-'CPT_Resistant';g[grepl('^CPT_',names(g))]<-'CPT_DTP';g[grepl('^CTRL',names(g))]<-'DMSO_CTRL';stopifnot(!anyNA(g))
d <- DGEList(cnt);keep<-rowSums(cpm(d)>=1)>=3;d<-d[keep,,keep.lib.sizes=FALSE];d<-calcNormFactors(d,method='TMM');lc<-cpm(d,log=TRUE,prior.count=1)
meta<-data.frame(sample=colnames(cnt),group=g,donor='POP66',patient_n_for_main_contrasts=1,raw_library_size=colSums(cnt),detected_genes=colSums(cnt>0),TMM_factor=d$samples$norm.factors)
write.table(meta,file.path(out,'sample_qc.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(gene_id=rownames(lc),symbol=x$Gene.symbol[match(rownames(lc),x$Ensembl.Gene.ID)],lc),gzfile(file.path(r,'03_processed','CRC_TMM_logCPM.tsv.gz')),sep='\t',row.names=FALSE,quote=FALSE)
mu<-sapply(unique(g),function(z)rowMeans(lc[,g==z,drop=FALSE]))
eff<-data.frame(gene_id=rownames(lc),symbol=x$Gene.symbol[match(rownames(lc),x$Ensembl.Gene.ID)],DTP_vs_CTRL=mu[,'CPT_DTP']-mu[,'DMSO_CTRL'],Regrowth_vs_DTP=mu[,'CPT_Regrowth']-mu[,'CPT_DTP'],Regrowth_vs_CTRL=mu[,'CPT_Regrowth']-mu[,'DMSO_CTRL'],Regrowth_vs_DMSO_Barcode=mu[,'CPT_Regrowth']-mu[,'DMSO_Barcode'],Regrowth_vs_Saline=mu[,'CPT_Regrowth']-mu[,'Saline'])
eff$state_reversal<-sign(eff$DTP_vs_CTRL)!=sign(eff$Regrowth_vs_DTP)&abs(eff$DTP_vs_CTRL)>=1
write.table(eff,file.path(out,'descriptive_gene_effects.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
v<-order(apply(lc,1,var),decreasing=TRUE);p<-prcomp(t(lc[v[1:min(2000,length(v))],]),scale.=FALSE)
write.table(data.frame(sample=rownames(p$x),group=g,p$x[,1:3]),file.path(out,'pca.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
cc<-cor(lc,method='spearman');write.table(cc,file.path(out,'sample_spearman.tsv'),sep='\t',quote=FALSE)
summary<-c(raw_genes=nrow(cnt),filtered_genes=nrow(lc),libraries=ncol(cnt),independent_patient_models_main_contrasts=1,DTP_abs_effect_ge1=sum(abs(eff$DTP_vs_CTRL)>=1),regrowth_CTRL_abs_effect_ge1=sum(abs(eff$Regrowth_vs_CTRL)>=1),median_abs_DTP=median(abs(eff$DTP_vs_CTRL)),median_abs_regrowth_CTRL=median(abs(eff$Regrowth_vs_CTRL)),DTP_regrowth_effect_correlation=cor(eff$DTP_vs_CTRL,eff$Regrowth_vs_DTP,method='spearman'))
write.table(data.frame(metric=names(summary),value=summary),file.path(out,'summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines('Descriptive effects only. No inferential P/FDR: principal contrasts derive from one patient model. Mouse samples are not independent human patients. No enhancer causality established.',file.path(out,'STATISTICAL_SCOPE.txt'))
capture.output(sessionInfo(),file=file.path(r,'logs','CRC_R_sessionInfo.txt'))
writeLines('exit_code=0; count_integrity=passed; outputs_written',file.path(out,'COMPLETE.txt'));print(summary)

