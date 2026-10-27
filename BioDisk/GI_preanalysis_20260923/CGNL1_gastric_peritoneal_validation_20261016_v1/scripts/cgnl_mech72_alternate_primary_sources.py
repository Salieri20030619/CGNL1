from pathlib import Path
import requests,json,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2';urls={'Yonemura_Elsevier.xml':'https://api.elsevier.com/content/article/PII:S221112472301625X?httpAccept=text/xml','DRA014736_DDBJ_listing.html':'https://ddbj.nig.ac.jp/public/ddbj_database/dra/fastq/DRA014/DRA014736/'}
def get(kv):
 n,u=kv;r={'file':n,'url':u}
 try:
  a=requests.get(u,timeout=(15,60));(O/'metadata'/n).write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content));print(n,a.status_code,a.content.decode('utf-8')[:1800],flush=True)
 except Exception as e:r['error']=repr(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:r=list(ex.map(get,urls.items()))
(O/'metadata/ALTERNATE_PRIMARY_SOURCES_72.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
