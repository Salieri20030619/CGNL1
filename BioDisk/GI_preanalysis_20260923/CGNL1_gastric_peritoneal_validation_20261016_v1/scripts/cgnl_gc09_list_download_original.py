from pathlib import Path
import json,subprocess,hashlib,datetime
from concurrent.futures import ThreadPoolExecutor
O=Path.cwd();folders=json.loads((O/'metadata/31patient_folders').read_text());records=[]
def listing(folder):
 u='https://data.mendeley.com/public-api/datasets/9t86hd78sj/files?folder_id='+folder['id']+'&version=1&$start=0&$limit=1000';p=O/'metadata'/('original31_files_'+folder['id']+'.json');r=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','10','--max-time','45','-H','Accept: application/vnd.mendeley-public-dataset.1+json','-sS','-w','%{http_code}','-o',str(p),u],capture_output=True,timeout=50);a=json.loads(p.read_text()) if r.returncode==0 else [];print('FOLDER',folder['name'],'FILES',len(a));print(json.dumps(a[:2],ensure_ascii=False));return {'folder':folder,'url':u,'exit':r.returncode,'http':r.stdout.decode(),'files':a}
with ThreadPoolExecutor(max_workers=3) as ex:lists=list(ex.map(listing,folders))
(O/'metadata/ORIGINAL31_FILE_LIST.json').write_text(json.dumps(lists,ensure_ascii=False,indent=2),encoding='utf-8')
allfiles=[dict(f,folder_name=r['folder']['name'],folder_id=r['folder']['id']) for r in lists for f in r['files']];print('ALL_FILES',json.dumps([{'name':x.get('filename'), 'size':x.get('size'), 'folder':x['folder_name'],'keys':list(x)} for x in allfiles],ensure_ascii=False))
