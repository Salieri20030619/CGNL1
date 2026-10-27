from pathlib import Path
import requests,json,hashlib,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2'
urls={
'GSE228955_counts_retry.gz':'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE228955&format=file&file=GSE228955_Normalized_counts_cpm.txt.gz',
'GSE228955_WT_retry.xlsx':'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE228955&format=file&file=GSE228955_Normalized_counts_cpm_project22.xlsx',
'GSE251950.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251950&targ=self&form=text&view=full',
}
def get(k,u):
 r={'file':k,'url':u}
 try:
  a=requests.get(u,timeout=(15,90));(O/'raw'/k).write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content),magic_hex=a.content[:8].hex(),sha256=hashlib.sha256(a.content).hexdigest());print(json.dumps(r),flush=True)
 except Exception as e:r['error']=str(e);print(k,str(e),flush=True)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(lambda kv:get(*kv),urls.items()))
(O/'metadata/ACQUISITION_09.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
