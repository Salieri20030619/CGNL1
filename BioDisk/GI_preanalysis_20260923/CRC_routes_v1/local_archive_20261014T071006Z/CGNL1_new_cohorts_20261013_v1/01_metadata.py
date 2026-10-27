from pathlib import Path
import requests,json,concurrent.futures,datetime
L=Path(__file__).resolve().parent;M=L/'metadata';M.mkdir(exist_ok=True)
jobs={
 'GSE200997_index':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200997/suppl/',
 'GSE132465_index':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132465/suppl/',
 'atlas_tree':'https://api.github.com/repos/icbi-lab/crc-atlas/git/trees/main?recursive=1',
 'atlas_zenodo':'https://zenodo.org/api/records/17878122',
 'EMTAB12022':'https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-12022',
 'GSE200997_soft':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200997/soft/GSE200997_family.soft.gz'}
def fetch(item):
 k,url=item
 try:
  p=requests.get(url,timeout=(12,45));p.raise_for_status()
  dest=M/(k+('.gz' if k.endswith('_soft') else '.txt'));dest.write_bytes(p.content)
  return dict(key=k,url=url,status=p.status_code,bytes=len(p.content),path=str(dest),preview=p.text[:2500] if k.endswith('_index') else '')
 except Exception as e:return dict(key=k,url=url,error=repr(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(fetch,jobs.items()))
(M/'FETCH_AUDIT.json').write_text(json.dumps(results,indent=2))
for result in results:print(json.dumps(result),flush=True)
