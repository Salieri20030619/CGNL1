import pathlib,json,subprocess,gzip,re,concurrent.futures,time,hashlib
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');p=O/'scripts/spatial_eligibility.py';s=p.read_text(encoding='utf-8-sig').replace("read_text()","read_text(encoding='utf-8')");p.write_text(s,encoding='utf-8')
subprocess.run([r'C:\Users\31937\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe',str(p)],check=True)
plan={'scRNA':'GSE225857 five matched colon patients, author Tu labels; assess QC and author subtype coverage before differential interpretation. All treated before surgery; no n<10 P/FDR.','clinical_map':{'s0107':{'site':'Left colon','regimen':'FOLFOX','cycles':3},'s0115':{'site':'Left colon','regimen':'FOLFOX','cycles':4},'s0813':{'site':'Right colon','regimen':'FOLFOXIRI','cycles':5},'s0920':{'site':'Left colon','regimen':'FOLFOX','cycles':5},'s1231':{'site':'Left colon','regimen':'FOLFOX','cycles':3}},'spatial':'GSE294385 eight primary sections from five patients. Author Primary Colon is pathology label; individual colon/rectum clinical mapping not available, provisional CRC discovery only.','spatial_contrast':'Author front tumor vs other tumor regions in SAME section; each region>=20 spots; average section differences within patient. Sensitivity comparator restricted to muscularis/subserosa. Separate front-stroma versus nonfront-stroma comparison.','scRNA_min_cells':20,'scRNA_sensitivity':'per-cell CP10K and author subtype comparisons; omit Tu10_COL3A1 stromal-like cells as sensitivity, no claiming confirmed malignant without independent CNV','programs':{'proliferation':['MKI67','TOP2A','CENPF','TYMS','PCNA'],'hypoxia':['CA9','VEGFA','SLC2A1','LDHA','PDK1','EGLN3'],'invasion_ECM':['LAMC2','ITGA5','TGFBI','SERPINE1','FN1','VIM'],'intestinal_differentiation':['MUC2','TFF3','CLCA1','SPINK4','KRT20','FABP1'],'stromal_proxy':['COL1A1','COL1A2','DCN','LUM'],'epithelial_proxy':['EPCAM','KRT8','KRT18','KRT19']},'screen':'Descriptive genomewide effect tables; require expressed>=1CPM in >=3 relevant samples and effect consistency to shortlist; n<10 no inferential ranking. No selected-gene validation claim.'}
(O/'ANALYSIS_PLAN.json').write_text(json.dumps(plan,indent=2))
s=gzip.open(O/'metadata/GSE294385.soft.gz','rt',encoding='utf-8').read();wanted={'M-ST-13','M-ST-15','M-ST-29','M-ST-30','M-ST-31','M-ST-32','M-ST-33','M-ST-34'};jobs=[]
for c in s.split('^SAMPLE = ')[1:]:
 title=re.search(r'!Sample_title = (.+)',c).group(1).strip()
 if title in wanted:
  u=re.search(r'!Sample_supplementary_file_\d+ = (ftp.+)',c).group(1).strip().replace('ftp://','https://');jobs.append((title,u))
assert len(jobs)==8
(O/'spatial_download_plan.json').write_text(json.dumps(jobs,indent=2))
def fetch(v):
 name,u=v;p=O/'raw'/(name+'.tar.gz');t=time.time()
 with (O/'raw'/(name+'.download.log')).open('w') as f:r=subprocess.run(['curl.exe','-x','http://127.0.0.1:7897','-L','--fail','--max-time','1800','--retry','2','-C','-','-o',str(p)+'.part',u],stdout=f,stderr=f)
 j={'sample':name,'url':u,'exit_code':r.returncode,'seconds':time.time()-t}
 if r.returncode==0:
  pathlib.Path(str(p)+'.part').replace(p);j['bytes']=p.stat().st_size;j['sha256']=hashlib.file_digest(p.open('rb'),'sha256').hexdigest();j['transfer_exit_code']=subprocess.run(['scp',str(p),'sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/raw/']).returncode
 (O/'raw'/(name+'.complete.json')).write_text(json.dumps(j,indent=2));print(json.dumps(j),flush=True);return j
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rr=list(ex.map(fetch,jobs))
(O/'raw/SPATIAL_DOWNLOAD_COMPLETE.json').write_text(json.dumps(rr,indent=2));print('SPATIAL DOWNLOAD COMPLETE',flush=True)
