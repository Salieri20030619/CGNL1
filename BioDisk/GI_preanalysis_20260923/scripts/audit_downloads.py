import pathlib,json,gzip,pandas as pd,xml.etree.ElementTree as ET
r=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
for p in (r/'logs').glob('*.download.json'):
 j=json.loads(p.read_text());print(p.name,j['status'],j.get('bytes'),j.get('exit_code'))
for p in (r/'logs').glob('*.fetch.log'):
 print(p.name,p.read_text()[-250:].replace('\r','\n'))
p=r/'01_metadata'/'GSE145356.soft.txt'
lines=p.read_text().splitlines()
for i,l in enumerate(lines):
 if l.startswith('!sample_table_begin'):print('SAMPLE_TABLE_HEAD',lines[i:i+4]);break
m=pd.read_csv(r/'02_raw'/'GSE234129_meta.tsv.gz',sep='\t');print('GC_COUNTS',pd.crosstab(m['sample'],m['celltype']).to_string());print('GC_TYPES',m.celltype.value_counts().to_string())
