from pathlib import Path
import json,re,os,subprocess,hashlib,time,concurrent.futures
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/v3_pilot_20261006');D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006/downloads')
env=os.environ.copy()
for k in list(env):
 if 'proxy' in k.lower():env.pop(k,None)
env.update(NO_PROXY='*',no_proxy='*')
s=(R/'metadata/iss_panel_2.tmap').read_text();m=re.search(r'jsonState\s*=\s*',s)
if m:
 d,_=json.JSONDecoder().raw_decode(s[m.end():]);(R/'metadata/ISS_PROJECT.json').write_text(json.dumps(d))
 print('ISS_PROJECT',{'mpp':d.get('mpp'),'layers':[x.get('name') for x in d.get('layers',[])]},flush=True)
z=json.loads((R/'metadata/zenodo_hd.json').read_text());print('ZENODO_LINKS',z.get('links'),flush=True)
latest=z.get('links',{}).get('latest')
tasks=[
('iss_cells.h5ad','https://crclm2.serve.scilifelab.se/iss_panel_2/prepared_cells_panel_2.h5ad'),
('iss_regions.json','https://crclm2.serve.scilifelab.se/iss_panel_2/regions_panel_2.json'),
('quiche_graph.py','https://raw.githubusercontent.com/jranek/quiche/main/src/quiche/tools/graph.py'),
('quiche_model.py','https://raw.githubusercontent.com/jranek/quiche/main/src/quiche/tools/quiche.py'),
('nichecompass_model.py','https://raw.githubusercontent.com/Lotfollahi-lab/nichecompass/main/src/nichecompass/models/nichecompass.py'),
('nichecompass_gp.py','https://raw.githubusercontent.com/Lotfollahi-lab/nichecompass/main/src/nichecompass/utils/gene_programs.py'),
('CGNL1_ncbi.html','https://www.ncbi.nlm.nih.gov/gene/84952')
]
if latest:tasks.append(('zenodo_latest.json',latest))
def run(it):
 n,u=it;p=D/n;tmp=p.with_suffix(p.suffix+'.part');start=time.time()
 cmd=['curl','-q','--proxy','','--noproxy','*','-L','--fail','--connect-timeout','12','--max-time','300','--max-filesize','2147483648','--retry','0','-o',str(tmp),u]
 with (R/'logs'/(n+'.curl.log')).open('w') as f:
  try:rc=subprocess.run(cmd,env=env,stdout=f,stderr=f,timeout=315).returncode
  except subprocess.TimeoutExpired:rc=124
 rec={'file':n,'url':u,'network_mode':'direct_only','curl_exit_code':rc,'elapsed_s':round(time.time()-start,2)}
 if rc==0:
  tmp.replace(p);rec.update(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),status='downloaded')
 else:rec['status']='failed'
 (R/'metadata'/(n+'.receipt.json')).write_text(json.dumps(rec,indent=2));print(json.dumps(rec),flush=True)
 return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:out=list(ex.map(run,tasks))
(R/'metadata/DIRECT_DATA_DOWNLOADS.json').write_text(json.dumps(out,indent=2))
print('DONE',flush=True)
