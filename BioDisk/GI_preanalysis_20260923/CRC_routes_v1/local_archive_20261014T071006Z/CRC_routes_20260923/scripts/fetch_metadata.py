import pathlib,json,subprocess,concurrent.futures,time
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923')
p={'date':'2026-09-23','route1':'Colon primary to liver metastatic malignant epithelial state and regulation','route2':'Primary colon invasive-front tumor-stroma interaction','scope':'First-stage data audit and empirical feasibility. No automatic three-gene commitment. Identify independent patients, colon vs rectum, sample matches, treatment, true malignant and pathology annotations before inferential claims.','priorities':['Download processed matrices and author annotations, use raw reads only if essential','Colon only primary analyses; unknown location and rectal kept separately','Patient is independent unit; n<10 exploratory effect/concordance without inferential tests','Genomewide BH for any genomewide tested family; no selected-gene FDR rescue','Within-patient/within-cell-type or pathology-matched contrasts; do not infer cell intrinsic mechanism from bulk site differences','Spatial spots/bins are not independent patients; ligand receptor predictions do not prove signaling','Avoid recycling previously published micro-metastasis stemness/dormancy signatures as novelty'],'paused':['GC candidate genes','CRC drug persister/regrowth'],'deliverables':['Data availability and donor eligibility audit','Reproducible small-scale analysis for each executable route','Evidence-based priority with concrete remaining gaps']}
(O/'PLAN.json').write_text(json.dumps(p,indent=2))
jobs=[]
for g in ['GSE294385','GSE231559','GSE225857','GSE285991','GSE178318']:
 b='https://ftp.ncbi.nlm.nih.gov/geo/series/'+g[:-3]+'nnn/'+g
 jobs.extend([(g+'.soft.gz',b+'/soft/'+g+'_family.soft.gz'),(g+'.files.html',b+'/suppl/')])
jobs += [('zenodo19043057.json','https://zenodo.org/api/records/19043057'),('spatial2026.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13390483/fullTextXML'),('GSE231559.html','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231559')]
def fetch(v):
 name,u=v;f=O/'metadata'/name;t=time.time();r=subprocess.run(['curl.exe','-x','http://127.0.0.1:7897','-L','--fail','--connect-timeout','20','--max-time','120','--retry','1','-o',str(f),u],capture_output=True);j={'file':name,'url':u,'exit_code':r.returncode,'bytes':f.stat().st_size if f.exists() else 0,'seconds':time.time()-t};print(json.dumps(j),flush=True);return j
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:r=list(ex.map(fetch,jobs))
(O/'metadata_downloads.json').write_text(json.dumps(r,indent=2));print('METADATA COMPLETE',flush=True)
