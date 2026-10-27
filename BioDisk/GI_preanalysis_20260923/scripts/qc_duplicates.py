import pandas as pd,pathlib,gzip
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
x=pd.read_csv(r/'02_raw'/'GSE145356_counts.txt.gz',sep='\t');d=x[x.duplicated('Ensembl.Gene.ID',keep=False)];print('DUPLICATES',d.shape,d.iloc[:12,:7].to_string(index=False));print('DUP_COUNT_DISAGREE',sum(z.iloc[:,3:].drop_duplicates().shape[0]>1 for _,z in d.groupby('Ensembl.Gene.ID')))
with gzip.open(r/'02_raw'/'GSM6042428_1-1.txt.gz','rt') as f:
 lines=f.readlines();idx=next(i for i,l in enumerate(lines) if l.startswith('FEATURES\t'));cols=lines[idx].strip().split('\t');print('RAW_PROBES',[[z.strip().split('\t')[cols.index(k)] for k in ['ProbeName','SystematicName']] for z in lines[idx+10:idx+20]])
