import pathlib,gzip,pandas as pd
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
p=r/'02_raw'/'GSE145356_counts.txt.gz'
if not p.exists():p=r/'01_metadata'/'GSE145356_counts.txt.gz'
df=pd.read_csv(p,sep='\t');print('CRC',df.shape,df.columns.tolist(),df.head(3).to_string(index=False))
m=pd.read_csv(r/'02_raw'/'GSE234129_meta.tsv.gz',sep='\t');b=pd.read_csv(r/'02_raw'/'GSE234129_barcodes.tsv.gz',sep='\t',header=None);print('GC_META',m.shape,'BARCODES',b.shape,'UNANNOTATED',len(set(b[0])-set(m.cell_barcodes)))
p=r/'01_metadata'/'GSE145356.soft.txt';t=p.read_text();print('SOFT_TAIL',t[-1000:])
