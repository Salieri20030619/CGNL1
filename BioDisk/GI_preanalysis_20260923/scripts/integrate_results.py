import pathlib,pandas as pd
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');o=r/'04_results'/'GC_bulk'
print(pd.read_csv(o/'composition_sensitivity_panels.tsv',sep='\t').to_string(index=False))
a=pd.read_csv(o/'paired_preservation_adjusted_DE.tsv',sep='\t');b=pd.read_csv(o/'composition_sensitivity_DE.tsv',sep='\t');c=pd.read_csv(o/'exclude_GC05_sensitivity_DE.tsv',sep='\t');f=pd.read_csv(o/'fresh_frozen_pair_effects.tsv',sep='\t');z=a.merge(b[['gene_id','logFC','adj.P.Val']],on='gene_id',suffixes=('','_composition')).merge(c[['gene_id','logFC','adj.P.Val']],on='gene_id',suffixes=('','_exclude05')).merge(f[['gene_id','mean_paired_log2CPM_difference']],on='gene_id');z.to_csv(o/'integrated_sensitivity_gene_table.tsv',sep='\t',index=False)
g=['CD36','FABP4','FABP5','SLC27A1','SLC27A2','SLC27A4','LDLR','LPL','AGPAT2','HACD1','ME1','ITGA5','ITGB1','PTK2','FASN','SCD','CPT1A']
print('CANDIDATE_PANEL',z[z.symbol.isin(g)][['symbol','logFC','P.Value','adj.P.Val','positive_pairs','negative_pairs','logFC_composition','adj.P.Val_composition','mean_paired_log2CPM_difference']].to_string(index=False))
p=r/'01_metadata'/'GC_ascites_SourceData_MOESM5.xlsx';x=pd.ExcelFile(p);print('SOURCE_DATA_SHEETS',x.sheet_names)
for sh in x.sheet_names:
 d=pd.read_excel(p,sheet_name=sh,header=None);print(sh,d.shape)
