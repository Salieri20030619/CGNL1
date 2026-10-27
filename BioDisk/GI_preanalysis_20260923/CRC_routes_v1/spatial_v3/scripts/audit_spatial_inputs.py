import urllib.request,pathlib,json,gzip,concurrent.futures,sys,re
import pandas as pd
sys.stdout.reconfigure(encoding='utf-8')
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3');B=O.parent
op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897','http':'http://127.0.0.1:7897'}))
jobs=[(g+'.soft.gz',f'https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}/soft/{g}_family.soft.gz') for g in ['GSE217414','GSE335552']]
jobs += [('BJC2026.html','https://www.nature.com/articles/s41416-026-03567-y'),('iss_panel_2.tmap','https://crclm2.serve.scilifelab.se/iss_panel_2.tmap')]
def run(j):
 name,url=j;p=O/'metadata'/name
 try:
  if not p.exists():
   with op.open(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=40) as f:b=f.read()
   p.write_bytes(b)
  return {'name':name,'url':url,'bytes':p.stat().st_size,'ok':True}
 except Exception as e:return {'name':name,'url':url,'ok':False,'error':str(e)}
res=list(concurrent.futures.ThreadPoolExecutor(max_workers=4).map(run,jobs));(O/'metadata/FETCH.json').write_text(json.dumps(res,indent=2))
for r in res:print(r)
for g in ['GSE217414','GSE335552']:
 p=O/'metadata'/(g+'.soft.gz')
 if not p.exists():continue
 s=gzip.decompress(p.read_bytes()).decode();keep=[x for x in s.splitlines() if x.startswith(('^SAMPLE','!Sample_title','!Sample_characteristics','!Sample_supplementary_file','!Series_supplementary_file','!Series_overall_design'))];print(g+'\n'+'\n'.join(keep))
if (O/'metadata/BJC2026.html').exists():
 s=(O/'metadata/BJC2026.html').read_text(encoding='utf-8');links=sorted(set(re.findall(r'https?[^\s<>"\']+\.(?:xlsx|pdf|docx)',s)));print('BJC SUPPLEMENTS',json.dumps(links))
m=pd.read_csv(B/'metadata/spatial_patient_map.tsv',sep='\t');v=m[m.visium & m.sample_type.str.contains('Liver')];print('EXISTING MAIN COHORT',len(v),'liver sections',v.patient.nunique(),'patients');v.to_csv(O/'results/GSE294385_liver_sections.tsv',sep='\t',index=False)
j=json.loads((B/'metadata/spatial_table2.json').read_text(encoding='utf-8'));print('CLINICAL TABLE KEYS',list(j))
