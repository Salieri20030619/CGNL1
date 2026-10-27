from pathlib import Path
import gzip,json
O=Path.cwd()/'mechanism_investigation_v2';p=O/'raw/GSE228955_counts_ftp.gz'
with gzip.open(p,'rt',encoding='utf-8') as f:
 header=f.readline();row=f.readline();n=1+sum(1 for _ in f)
print('HEADER',header);print('FIRST_ROW',row);print('ROWS',n)
(O/'metadata/KO_MATRIX_SCHEMA.json').write_text(json.dumps({'header':header,'first_row':row,'rows':n},indent=2),encoding='utf-8')
