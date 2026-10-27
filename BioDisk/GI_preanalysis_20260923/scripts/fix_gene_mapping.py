import pathlib,pandas as pd
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
x=pd.read_csv(r/'02_raw'/'GSE145356_counts.txt.gz',sep='\t',usecols=['Ensembl.Gene.ID','Gene.symbol']);x=x.drop_duplicates('Ensembl.Gene.ID');x.columns=['gene_id','symbol'];x.to_csv(r/'01_metadata'/'human_gene_symbols_from_GSE145356.tsv',sep='\t',index=False)
p=r/'scripts'/'gc_bulk_preanalysis.R';s=p.read_text();s=s.replace('human_gene_symbols_from_multiome.tsv','human_gene_symbols_from_GSE145356.tsv');p.write_text(s)
