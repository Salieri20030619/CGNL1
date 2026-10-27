from pathlib import Path
import json,gzip,io,re,subprocess,datetime
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
O=Path.cwd();B=O.parent
def read_table(p,nrows=None,header='infer'):
 b=p.read_bytes();isgz=b[:2]==b'\x1f\x8b';raw=gzip.decompress(b) if isgz else b;txt=raw.decode('utf-8',errors='strict');print('FORMAT',p.name,'gzip' if isgz else 'plain_text_despite_suffix',txt[:180]);return pd.read_csv(io.StringIO(txt),sep='\t',header=header,nrows=nrows),isgz
f,gz=read_table(O/'metadata/GSE239676_features.tsv.gz',header=None);hit=f[f.astype(str).apply(lambda c:c.str.contains('CGNL1|ENSG00000128849')).any(axis=1)];print('FEATURES',len(f),'CGNL1',hit.to_dict('records'))
fq={'feature_rows':len(f),'target_rows':hit.to_dict('records'),'actual_gzip':gz,'CGNL1_measured':len(hit)>0};(O/'metadata/GSE239676_TARGET_QUALIFICATION.json').write_text(json.dumps(fq,indent=2),encoding='utf-8')
try:
 m,gz=read_table(O/'metadata/GSE239676_meta.tsv.gz');print('239META',m.shape,list(m),m.head(2).to_string(index=False));print('239_CATEGORIES',{c:m[c].value_counts().to_dict() for c in m if m[c].nunique()<80});m.head(5).to_csv(O/'metadata/GSE239676_preview.tsv',sep='\t',index=False)
except Exception as e:print('239META_PENDING',str(e))
for p in (B/'01_metadata').iterdir():
 if '314812' in p.name or 'GC' in p.name:print('BULK_META',p.name)
for rel in ['02_raw/GSE314812_counts.tsv.gz','03_processed/GC314812_TMM_logCPM.tsv.gz','GC_validation_v2/candidate_discovery/anatomy_proxy_sample_scores.tsv','GC_validation_v2/04_results/paired_scRNA/pseudobulk_design.tsv']:
 p=B/rel
 if p.exists():
  d=pd.read_csv(p,sep='\t',nrows=2);print('SCHEMA',rel,d.to_string(index=False))
sources={'mendeley_files_v1':'https://data.mendeley.com/public-api/datasets/9t86hd78sj/files?version=1','mendeley_files_root':'https://data.mendeley.com/public-api/datasets/9t86hd78sj/files?folder_id=root&version=1','mendeley_v2_files':'https://data.mendeley.com/api/datasets-v2/datasets/9t86hd78sj/files?folder_id=root&version=1'}
def fetch(kv):
 name,url=kv;p=O/'metadata'/name;r=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','10','--max-time','40','-sS','-w','%{http_code}','-o',str(p),url],capture_output=True,timeout=45);rec={'url':url,'exit':r.returncode,'http':r.stdout.decode('utf-8'),'file':str(p)};print('MENDELEY',json.dumps(rec),p.read_text(encoding='utf-8')[:6000] if p.exists() else '');return rec
with ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(fetch,sources.items()))
(O/'metadata/MENDELEY_ROUTE_RECEIPTS.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
