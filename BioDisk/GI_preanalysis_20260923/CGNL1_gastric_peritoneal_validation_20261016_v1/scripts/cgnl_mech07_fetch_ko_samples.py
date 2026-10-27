from pathlib import Path
import requests,json,concurrent.futures,hashlib
O=Path.cwd()/'mechanism_investigation_v2'
urls={
 'GSE228955_samples.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228955&targ=gsm&form=text&view=full',
 'GSE228955_Normalized_counts_cpm_project22.xlsx':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228955/suppl/GSE228955_Normalized_counts_cpm_project22.xlsx',
 'GSE228955_Normalized_counts_cpm.txt.gz':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228955/suppl/GSE228955_Normalized_counts_cpm.txt.gz',
 'GSE150833.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150833&targ=self&form=text&view=full',
 'GSE76162.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76162&targ=self&form=text&view=full',
}
def fetch(k,u):
 r={'file':k,'url':u}
 try:
  a=requests.get(u,timeout=(15,60));(O/'raw'/k).write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content),sha256=hashlib.sha256(a.content).hexdigest());print(k,a.status_code,len(a.content),flush=True)
 except Exception as e:r['error']=str(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:r=list(ex.map(lambda kv:fetch(*kv),urls.items()))
(O/'metadata/ACQUISITION_07.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
