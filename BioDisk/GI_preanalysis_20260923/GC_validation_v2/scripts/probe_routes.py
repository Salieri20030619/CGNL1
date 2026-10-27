import pathlib,pandas as pd,json,gzip
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');V=B/'GC_validation_v2';O=V/'candidate_routes';O.mkdir(exist_ok=True)
m=pd.read_csv(V/'03_processed/GSE234129_annotation_qc.tsv.gz',sep='\t');print(m.columns.tolist());print(m.groupby(['patient','sample','site']).size().to_string());print(pd.crosstab(m.celltype,m.site).to_string());print(m.groupby(['sample','broad_type']).size().unstack(fill_value=0).to_string());print(m.head(2).to_string(index=False))
q=pd.read_csv(V/'03_processed/GSM7133741_P1_gating_refined.tsv.gz',sep='\t');print(q.columns.tolist());print(q.head(2).to_string(index=False))
