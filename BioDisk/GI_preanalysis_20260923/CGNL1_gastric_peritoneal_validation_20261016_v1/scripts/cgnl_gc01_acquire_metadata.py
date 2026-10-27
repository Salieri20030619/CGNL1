from pathlib import Path
import json,subprocess,hashlib,gzip,re,datetime
from concurrent.futures import ThreadPoolExecutor
from bs4 import BeautifulSoup
O=Path.cwd();src={
'mendeley_page':'https://data.mendeley.com/datasets/9t86hd78sj/1',
'mendeley_api':'https://data.mendeley.com/public-api/datasets/9t86hd78sj/versions/1',
'mendeley_files':'https://data.mendeley.com/public-api/datasets/9t86hd78sj/versions/1/files',
'GSE239676_geo':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239676&targ=self&form=text&view=full',
'GSE239676_features.tsv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239676/suppl/GSE239676_features.tsv.gz',
'GSE239676_meta.tsv.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE239nnn/GSE239676/suppl/GSE239676_meta.tsv.gz',
'GSE228598_geo':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228598&targ=self&form=text&view=full',
'GSE239676_author_readme':'https://raw.githubusercontent.com/Coolgenome/GastricCancer/main/README.md',
}
def fetch(kv):
 name,url=kv;p=O/'metadata'/name;proc=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','15','--max-time','90','--retry','1','-sS','-w','%{http_code}','-o',str(p),url],capture_output=True,timeout=195);r={'name':name,'url':url,'exit':proc.returncode,'http':proc.stdout.decode('utf-8'),'stderr':proc.stderr.decode('utf-8'),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 if p.exists():b=p.read_bytes();r.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
 print(json.dumps(r),flush=True);return r
with ThreadPoolExecutor(max_workers=4) as ex:rec=list(ex.map(fetch,src.items()))
(O/'metadata/DOWNLOAD_RECEIPTS_01.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
for name in ['mendeley_api','mendeley_files','GSE239676_geo','GSE228598_geo']:
 p=O/'metadata'/name;txt=p.read_text(encoding='utf-8');print('FILE',name,txt[:7000])
page=(O/'metadata/mendeley_page').read_text(encoding='utf-8');soup=BeautifulSoup(page,'html.parser');data=[x.get_text() for x in soup.find_all('script') if '9t86hd78sj' in x.get_text()];(O/'metadata/mendeley_embedded.txt').write_text('\n'.join(data),encoding='utf-8');print('MENDELEY_EMBEDDED',str(data)[:5000])
import pandas as pd
for name in ['GSE239676_features.tsv.gz','GSE239676_meta.tsv.gz']:
 try:
  d=pd.read_csv(O/'metadata'/name,sep='\t',nrows=5,header=None if 'features' in name else 'infer');print('PREVIEW',name,d.to_string(index=False))
  if 'features' in name:
   d=pd.read_csv(O/'metadata'/name,sep='\t',header=None);print('CGNL1',d[d.astype(str).apply(lambda x:x.str.contains('CGNL1|ENSG00000128849',regex=True)).any(axis=1)].to_string(index=False));print('FEATURE_ROWS',len(d))
 except Exception as e:print('PARSE_FAILURE',name,str(e))
old=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2');print('OLD_SCRIPTS',[p.name for p in (old/'scripts').iterdir()]);print('OLD_METADATA',[p.name for p in (old/'01_metadata').iterdir()]);print('OLD_RAW_MANIFEST',json.loads((old/'logs/GSE228598_complete.json').read_text()))
