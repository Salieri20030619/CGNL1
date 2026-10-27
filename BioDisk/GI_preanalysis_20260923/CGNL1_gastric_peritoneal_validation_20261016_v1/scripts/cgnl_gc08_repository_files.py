from pathlib import Path
import json,subprocess,datetime,hashlib
from concurrent.futures import ThreadPoolExecutor
O=Path.cwd()
sources={'31patient_folders':'https://data.mendeley.com/public-api/datasets/9t86hd78sj/folders/1','Gut2025_folders':'https://data.mendeley.com/public-api/datasets/jwkc5t6r55/folders/1','Gut2025_root':'https://data.mendeley.com/public-api/datasets/jwkc5t6r55/files?folder_id=root&version=1&$start=0&$limit=1000','31patient_snapshot':'https://data.mendeley.com/public-api/datasets/9t86hd78sj/snapshot/1'}
def fetch(item):
 name,url=item;p=O/'metadata'/name;r=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','10','--max-time','45','-H','Accept: application/vnd.mendeley-public-dataset.1+json','-sS','-w','%{http_code}','-o',str(p),url],capture_output=True,timeout=50);rec={'name':name,'url':url,'exit':r.returncode,'http':r.stdout.decode('utf-8')}
 if p.exists():rec['sha256']=hashlib.sha256(p.read_bytes()).hexdigest();print(name,p.read_text(encoding='utf-8')[:13000])
 return rec
with ThreadPoolExecutor(max_workers=4) as ex:rec=list(ex.map(fetch,sources.items()))
(O/'metadata/REPOSITORY_FILES_RECEIPTS.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
