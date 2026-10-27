from pathlib import Path
import requests,json,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2'
tasks={'DRA014736_runs.tsv':('https://www.ebi.ac.uk/ena/portal/api/search',{'result':'read_run','query':'submission_accession="DRA014736"','fields':'run_accession,study_accession,sample_accession,experiment_accession,sample_title,library_strategy,library_layout,fastq_ftp,fastq_bytes,fastq_md5','format':'tsv','limit':'0'}),'DRA014736.xml':('https://www.ebi.ac.uk/ena/browser/api/xml/DRA014736',None)}
def fetch(kv):
 n,(u,p)=kv;r={'file':n,'url':u,'params':p}
 try:
  a=requests.get(u,params=p,timeout=(15,60));(O/'metadata'/n).write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content));print(n,a.status_code,a.content.decode('utf-8')[:5500],flush=True)
 except Exception as e:r['error']=str(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:r=list(ex.map(fetch,tasks.items()))
(O/'metadata/ACQUISITION_36.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
