from pathlib import Path
import requests,json,re,hashlib,concurrent.futures
ROOT=Path(__file__).resolve().parent
(ROOT/'sources').mkdir(exist_ok=True)
BASE='https://singlecell.broadinstitute.org'
urls={'portal':BASE+'/single_cell/study/SCP2079/combined-pd-1-braf-and-mek-inhibition-in-braf-v600e-colorectal-cancer-a-phase-2-trial','study_api':BASE+'/single_cell/api/v1/studies/SCP2079','files_api':BASE+'/single_cell/api/v1/studies/SCP2079/study_files','site_api':BASE+'/single_cell/api/v1/site/studies/SCP2079'}
def fetch(k,u):
 try:
  r=requests.get(u,timeout=35)
  (ROOT/'sources'/f'{k}.txt').write_bytes(r.content)
  d=dict(name=k,url=u,final_url=r.url,status=r.status_code,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest())
  print(json.dumps(d),flush=True)
  if k=='portal':
   for line in r.text.splitlines():
    if re.search('epi_|bucket|studyFiles|study_files|download|firecloud',line,re.I):print(line[:1500],flush=True)
  else:print(r.text[:1800],flush=True)
  return d
 except Exception as e:return dict(name=k,error=str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: rows=list(ex.map(lambda kv:fetch(*kv),urls.items()))
(ROOT/'SOURCE_FETCH.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')