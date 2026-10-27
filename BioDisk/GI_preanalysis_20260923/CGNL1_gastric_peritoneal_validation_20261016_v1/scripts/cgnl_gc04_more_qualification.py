from pathlib import Path
import json,re,subprocess,datetime,gzip,io
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from bs4 import BeautifulSoup
O=Path.cwd();B=O.parent
d=pd.read_csv(B/'01_metadata/GSE314812_paired_design.tsv',sep='\t');print('PAIRED_DESIGN',d.to_string(index=False))
m=pd.read_csv(B/'GC_validation_v2/candidate_discovery/anatomy_proxy_sample_scores.tsv',sep='\t');print('PRESERVATION',pd.crosstab(m.site,m.preservation).to_string());print('SAMPLE_METADATA',list(m))
meta=pd.read_csv(B/'02_raw/GSE234129_meta.tsv.gz',sep='\t');print('234_META',meta.shape,list(meta));print('234_CATS',{c:meta[c].value_counts().to_dict() for c in meta if meta[c].nunique()<40})
src={'mendeley_folders':'https://data.mendeley.com/public-api/datasets/9t86hd78sj/folders?version=1','mendeley_empty_root':'https://data.mendeley.com/public-api/datasets/9t86hd78sj/files?folder_id=&version=1','GSE314812_geo':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314812&targ=self&form=text&view=full','GSE234129_geo':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234129&targ=self&form=text&view=full'}
page=BeautifulSoup((O/'metadata/mendeley_page').read_text(encoding='utf-8'),'html.parser');scripts=[x.get('src') for x in page.find_all('script') if x.get('src')];print('MENDELEY_JS',scripts)
def fetch(kv):
 name,url=kv;p=O/'metadata'/name;r=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','10','--max-time','35','-sS','-w','%{http_code}','-o',str(p),url],capture_output=True,timeout=40);rec={'name':name,'url':url,'exit':r.returncode,'http':r.stdout.decode('utf-8')};print(json.dumps(rec),p.read_text(encoding='utf-8')[:8000] if p.exists() else '');return rec
with ThreadPoolExecutor(max_workers=4) as ex:rr=list(ex.map(fetch,src.items()))
(O/'metadata/DOWNLOAD_RECEIPTS_04.json').write_text(json.dumps(rr,indent=2),encoding='utf-8')
