from pathlib import Path
import json,re,subprocess,datetime,hashlib
O=Path.cwd();txt=(O/'metadata/mendeley_js_7').read_text(encoding='utf-8');snips=[]
for pat in ['fetchFolders','folders?','/folders','folder-tree','folderTree','folder_id=null','root', 'publicApiBaseUrl']:
 snips += [txt[max(0,m.start()-180):m.start()+350] for m in list(re.finditer(re.escape(pat),txt))[:5]]
(O/'metadata/MENDELEY_JS_ROUTES.json').write_text(json.dumps(snips,indent=2),encoding='utf-8');print('ROUTES',json.dumps(snips))
url='https://data.mendeley.com/public-api/datasets/9t86hd78sj/files?folder_id=root&version=1&$start=0&$limit=1000';p=O/'metadata/mendeley_official_headers_files';r=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','10','--max-time','40','-H','Accept: application/vnd.mendeley-public-dataset.1+json','-sS','-w','%{http_code}','-o',str(p),url],capture_output=True,timeout=45);print('OFFICIAL_FILES',r.returncode,r.stdout.decode(),p.read_text(encoding='utf-8'))
