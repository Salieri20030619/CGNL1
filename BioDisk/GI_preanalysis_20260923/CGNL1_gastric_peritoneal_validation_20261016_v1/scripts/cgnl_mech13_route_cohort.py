from pathlib import Path
import requests,json,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2'
urls={'GSE237876.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237876&targ=self&form=text&view=full','GSE237876_samples.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237876&targ=gsm&form=text&view=full','GSE237876_files.html':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE237nnn/GSE237876/suppl/','GSE165211.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165211&targ=self&form=text&view=full'}
def get(k,u):
 r={'file':k,'url':u}
 try:
  a=requests.get(u,timeout=(15,50));(O/'metadata'/k).write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content));print(k,a.status_code,len(a.content),flush=True)
 except Exception as e:r['error']=str(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(lambda kv:get(*kv),urls.items()))
(O/'metadata/ACQUISITION_13.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
