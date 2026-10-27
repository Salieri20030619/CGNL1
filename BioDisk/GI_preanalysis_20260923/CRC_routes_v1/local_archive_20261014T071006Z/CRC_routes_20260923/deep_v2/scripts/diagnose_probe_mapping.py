import pathlib,pandas as pd,re
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2');a=pd.read_csv(O/'metadata/GPL6370_annotation.tsv.gz',sep='\t',dtype=str);print(a.SPOT_ID.head(35).tolist());print('PREFIX',a.SPOT_ID.str.extract(r'^([^_]+)')[0].value_counts().head(15).to_dict());print('REFSEQ',a.SPOT_ID.str.contains(r'[NX][MPR]_\d+',na=False).sum());print(a[a.SPOT_ID.str.contains('CGNL1|KIAA|115573|032866',na=False)].head(15).to_string(index=False));
import gzip
s=gzip.open(O/'metadata/GSE14297.soft.gz','rt',encoding='utf-8').read();print('\n'.join(x for x in s.splitlines() if x.startswith(('!Platform_title','!Platform_description','!Platform_manufacturer'))))
