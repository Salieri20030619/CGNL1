from pathlib import Path
import json,subprocess,hashlib,datetime,re,os,signal
from concurrent.futures import ThreadPoolExecutor
O=Path.cwd();D=O/'downloads';records=[]
# Retire only the stalled task-owned requests collector; bounded curl probes supersede it.
pid=2120424;proc=Path('/proc')/str(pid)
if proc.exists():
 args=(proc/'cmdline').read_bytes().decode('utf-8').replace('\0',' ')
 if str(O/'scripts/cgnl_06_literature_repair.py') in args:
  os.kill(pid,signal.SIGTERM);records.append({'action':'terminate_stalled_own_fetch_collector','pid':pid,'reason':'unbounded DNS/TLS behavior of requests; response files preserved, bounded curl records available','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
urls={'HGP_cells_range':'https://crclm2.serve.scilifelab.se/iss_panel_2/prepared_cells_panel_2.h5ad','HGP_cells_download_range':'https://crclm2.serve.scilifelab.se/iss_panel_2/prepared_cells_panel_2.h5ad?dl=1','HGP_regions_range':'https://crclm2.serve.scilifelab.se/iss_panel_2/regions_panel_2.json','HGP_archive_range':'https://crclm2.serve.scilifelab.se/iss_panel_2.tmap?dl=1','HGP_repo_readme':'https://raw.githubusercontent.com/axanderssonumu/iss_crclm/main/README.md','HGP_repo_analysis':'https://raw.githubusercontent.com/axanderssonumu/iss_crclm/main/analyses.ipynb','GSE206552_current':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206552&targ=self&form=text&view=full','GSE254364_current':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254364&targ=self&form=text&view=full'}
def fetch(item):
 name,url=item;p=D/(name+'.response');headers=D/(name+'.headers');cmd=['curl','--noproxy','*','-L','--connect-timeout','10','--max-time','40','--max-filesize','3000000','-sS','-D',str(headers),'-w','%{http_code}','-o',str(p)]
 if 'range' in name:cmd+=['-r','0-8191']
 r=subprocess.run(cmd+[url],capture_output=True,timeout=45);rec={'name':name,'url':url,'exit_code':r.returncode,'http':r.stdout.decode('utf-8'),'error':r.stderr.decode('utf-8'),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 if p.exists():
  b=p.read_bytes();rec.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),HDF5_magic=b.startswith(b'\x89HDF\r\n\x1a\n'),HTML_response=b.lstrip().lower().startswith((b'<!doctype html',b'<html')))
  if name.startswith('HGP_repo'):
   txt=b.decode('utf-8');rec['CGNL1_mentions']=len(re.findall('CGNL1',txt,re.I));rec['data_paths']=re.findall(r'[^\s\"\']*\.h5ad',txt)[:20]
  if name.startswith('GSE'):rec['source_description']='\n'.join(x for x in b.decode('utf-8').splitlines() if any(k in x for k in ['Series_title','Series_summary','Series_overall_design','Series_pubmed_id','Series_sample_id']))
 return rec
with ThreadPoolExecutor(max_workers=4) as pool:records.extend(pool.map(fetch,urls.items()))
(O/'metadata/PUBLIC_FILE_PROBES.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
for r in records:print(json.dumps(r,ensure_ascii=False))
