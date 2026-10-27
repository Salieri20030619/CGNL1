r <- '/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923'
cat(R.version.string, '\n')
for (p in c('edgeR','limma','ggplot2','data.table')) cat(p,requireNamespace(p,quietly=TRUE),'\n')
