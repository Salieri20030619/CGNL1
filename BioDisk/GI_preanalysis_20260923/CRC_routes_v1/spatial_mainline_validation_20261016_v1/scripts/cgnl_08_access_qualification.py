from pathlib import Path
import json,re,subprocess,hashlib,datetime
from concurrent.futures import ThreadPoolExecutor
O=Path.cwd();D=O/'downloads';L=O/'literature'
text=(L/'HGP_project.response').read_text(encoding='utf-8');m=re.search(r'jsonState\s*=\s*(\{.*?\})\s*;?\s*projectUtils.loadProject',text,re.S);j=json.loads(m.group(1));(O/'metadata/HGP_PROJECT_DECODED.json').write_text(json.dumps(j,indent=2),encoding='utf-8')
sources={'zenodo_current':'https://zenodo.org/api/records/19043057/versions/latest','zenodo_initial':'https://zenodo.org/api/records/17559845','HGP_worker':'https://crclm2.serve.scilifelab.se/static/js/utils/h5Utils_worker.js','HGP_h5Utils':'https://crclm2.serve.scilifelab.se/static/js/utils/h5Utils.js','HGP_projectUtils':'https://crclm2.serve.scilifelab.se/static/js/utils/projectUtils.js','HGP_code_tree':'https://api.github.com/repos/axanderssonumu/iss_crclm/git/trees/main?recursive=1','CRC_source_tree':'https://api.github.com/repos/yliuup/CRC_micromets_ST/git/trees/main?recursive=1','HGP_EJSO_2026':'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pubmed.cgi/BioC_xml/42296573/unicode'}
def fetch(item):
 name,url=item;p=D/(name+'.response');r={'name':name,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 proc=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','10','--max-time','45','--retry','0','-sS','-w','%{http_code}','-o',str(p),url],capture_output=True,timeout=50)
 r.update(exit_code=proc.returncode,http=proc.stdout.decode('utf-8'),error=proc.stderr.decode('utf-8'))
 if p.exists():
  b=p.read_bytes();r.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
  if proc.returncode==0:
   if name.startswith('zenodo'):
    try:a=json.loads(b);r.update(id=a.get('id'),access=a.get('metadata',{}).get('access_right',a.get('access')),files=a.get('files',[]))
    except Exception as e:r['parse_error']=str(e)
   if name.endswith('tree'):
    try:a=json.loads(b);r['paths']=[x['path'] for x in a.get('tree',[])];r['truncated']=a.get('truncated')
    except Exception as e:r['parse_error']=str(e)
   if name.startswith('HGP_') and name in ['HGP_worker','HGP_h5Utils','HGP_projectUtils']:
    txt=b.decode('utf-8');r['route_lines']=[line for line in txt.splitlines() if any(k in line for k in ['fetch(','.h5','download','h5ad','/h5','createFile','relative_root'])][:35]
 print('FETCH',name,proc.returncode,r['http'],flush=True);return r
with ThreadPoolExecutor(max_workers=4) as pool:receipts=list(pool.map(fetch,sources.items()))
(O/'metadata/ACCESS_QUALIFICATION.json').write_text(json.dumps(receipts,ensure_ascii=False,indent=2),encoding='utf-8')
for r in receipts:print(json.dumps(r,ensure_ascii=False))
