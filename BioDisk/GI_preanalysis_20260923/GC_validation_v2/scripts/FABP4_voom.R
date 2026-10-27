suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(limma))
B<-'/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923';O<-file.path(B,'GC_validation_v2','04_results','FABP4_audit')
writeLines(c('Follow-up prompted by OLS result: replicate original limma-voom primary model then add adipocyte proxy; additionally adipocyte plus myeloid.', 'Use original retained gene universe and original gene IDs; BH across all retained genes within each model.', 'Report all models; no selection by significance. These are post-hoc sensitivities, not external validation.'),file.path(O,'VOOM_PLAN.txt'))
m<-read.delim(file.path(O,'sample_values.tsv'),check.names=FALSE);m$patient<-factor(m$patient);m$site<-factor(m$site,levels=c('GCPT','GCPM'));m$preservation<-factor(m$preservation)
x<-read.delim(gzfile(file.path(B,'02_raw','GSE314812_counts.tsv.gz')),check.names=FALSE);cnt<-as.matrix(x[,m$sample]);rownames(cnt)<-x$gene_id
orig<-read.delim(file.path(B,'04_results','GC_bulk','paired_preservation_adjusted_DE.tsv'),check.names=FALSE)
d<-calcNormFactors(DGEList(cnt[match(orig$gene_id,rownames(cnt)),,drop=FALSE]))
forms<-list(primary=~patient+preservation+site,adipocyte=~patient+preservation+adipocyte+site,adipocyte_myeloid=~patient+preservation+adipocyte+myeloid+site)
out<-list()
for(k in names(forms)){de<-model.matrix(forms[[k]],m);stopifnot(qr(de)$rank==ncol(de));v<-voom(d,de,plot=FALSE);f<-eBayes(lmFit(v,de));tt<-topTable(f,coef='siteGCPM',number=Inf,sort.by='none',confint=TRUE);tt$gene_id<-rownames(tt);tt$symbol<-orig$symbol[match(tt$gene_id,orig$gene_id)];write.table(tt,gzfile(file.path(O,paste0('voom_',k,'_all_genes.tsv.gz'))),sep='\t',row.names=FALSE,quote=FALSE);ii<-which(tt$symbol=='FABP4');out[[k]]<-data.frame(model=k,tt[ii,],n_genes=nrow(tt));if(k=='primary'){z<-match(orig$gene_id,tt$gene_id);stopifnot(max(abs(tt$logFC[z]-orig$logFC))<1e-8,max(abs(tt$adj.P.Val[z]-orig$adj.P.Val))<1e-8)}}
y<-do.call(rbind,out);write.table(y,file.path(O,'FABP4_voom_model_comparison.tsv'),sep='\t',row.names=FALSE,quote=FALSE);print(y);writeLines('exit_code=0; primary_reproduced_logFC_and_FDR_tolerance_1e-8; all_models_full_rank',file.path(O,'VOOM_COMPLETE.txt'))
