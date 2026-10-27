from pathlib import Path
import json,re,gzip
import pandas as pd
from bs4 import BeautifulSoup
O=Path.cwd();B=O.parent
page=(O/'metadata/mendeley_page').read_text(encoding='utf-8');soup=BeautifulSoup(page,'html.parser')
for s in soup.find_all('script'):
 text=s.get_text()
 if '9t86hd78sj' in text:
  (O/'metadata/mendeley_embedded_state.txt').write_text(text,encoding='utf-8');print('EMBEDDED_HEAD',text[:1200])
  for pat in ['download','files','filename','content_details','9t86hd78sj','zip','csv','xlsx']:
   print('PATTERN',pat,[text[max(0,m.start()-100):m.start()+450] for m in list(re.finditer(re.escape(pat),text))[:5]])
f=pd.read_csv(O/'metadata/GSE239676_features.tsv.gz',sep='\t',header=None);print('FEATURES',f.shape, f.head().to_string(index=False));print('CGNL1',f[f.astype(str).apply(lambda c:c.str.contains('CGNL1|ENSG00000128849')).any(axis=1)].to_string(index=False))
try:
 m=pd.read_csv(O/'metadata/GSE239676_meta.tsv.gz',sep='\t');print('239META',m.shape,list(m),m.head(3).to_string(index=False));print('239_CATEGORIES',{c:m[c].value_counts().to_dict() for c in m if m[c].nunique()<80});m.groupby([c for c in m if any(x in c.lower() for x in ['patient','sample','tissue','celltype','type','site'])]).size().to_csv(O/'metadata/GSE239676_coverage.tsv',sep='\t')
except Exception as e:print('239_META_NOT_READY',str(e))
for p in (B/'01_metadata').iterdir():
 if '314812' in p.name or 'GC' in p.name:print('BULK_META',p.name)
for rel in ['02_raw/GSE314812_counts.tsv.gz','03_processed/GC314812_TMM_logCPM.tsv.gz','GC_validation_v2/candidate_discovery/anatomy_proxy_sample_scores.tsv','GC_validation_v2/04_results/paired_scRNA/pseudobulk_design.tsv']:
 p=B/rel
 if p.exists():
  d=pd.read_csv(p,sep='\t',nrows=3);print('SCHEMA',rel,d.to_string(index=False))
