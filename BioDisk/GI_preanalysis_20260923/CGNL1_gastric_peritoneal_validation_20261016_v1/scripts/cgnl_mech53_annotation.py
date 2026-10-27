from pathlib import Path
import requests,json,concurrent.futures,hashlib
O=Path.cwd()/'mechanism_investigation_v2'
tasks={'GPL570_full.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570&targ=self&form=text&view=full','HPA_CGNL1.xml':'https://www.proteinatlas.org/ENSG00000128849.xml'}
def get(kv):
 n,u=kv;r={'file':n,'url':u}
 try:
  a=requests.get(u,timeout=(15,90));(O/'metadata'/n).write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content),sha256=hashlib.sha256(a.content).hexdigest());print(n,a.status_code,len(a.content),flush=True)
 except Exception as e:r['error']=repr(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:r=list(ex.map(get,tasks.items()))
(O/'metadata/ACQUISITION_53.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
