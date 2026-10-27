from pathlib import Path
import json,gzip,re,hashlib,datetime,sys,shutil
sys.stdout.reconfigure(encoding='utf-8')
S=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3');B=S/'background_validation';rid='run_'+datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');R=B/rid;R.mkdir(parents=True);(R/'scripts').mkdir();remote='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/spatial_v3/background_validation/'+rid
samples=[];files=[]
for g in ['GSE231559','GSE285991']:
 raw=S.parent/'metadata'/(g+'.soft.gz');shutil.copy2(raw,R/raw.name);s=gzip.decompress(raw.read_bytes()).decode()
 for block in s.split('^SAMPLE = ')[1:]:
  gsm=block.splitlines()[0].strip();title=re.search(r'!Sample_title = (.+)',block).group(1).strip();species=re.search(r'!Sample_organism_ch1 = (.+)',block).group(1).strip()
  if species!='Homo sapiens':continue
  if g=='GSE285991' and not re.fullmatch(r'P\d+_LM',title):continue
  chars=re.findall(r'!Sample_characteristics_ch1 = (.+)',block);urls=re.findall(r'!Sample_supplementary_file(?:_\d+)? = (ftp[^\n]+)',block)
  pid=(re.match(r'[CL]\d+',title).group() if g=='GSE231559' else title.split('_')[0]);sample={'dataset':g,'sample':gsm,'title':title,'patient_label_provisional':g+':'+pid,'species':species,'characteristics':chars,'pairing':'Not assumed between C and L labels','colon_rectum':'not individually confirmed','files':[]}
  for u in urls:
   u=u.strip().replace('ftp://','https://');name=u.rsplit('/',1)[-1];rel='inputs/'+g+'/'+name;files.append({'dataset':g,'sample':gsm,'name':name,'url':u,'relative_path':rel});sample['files'].append(rel)
  samples.append(sample)
for ident in [24803571,24803616]:
 j=json.loads((S/'metadata'/('figshare_'+str(ident)+'.json')).read_text());
 for f in j['files']:files.append({'dataset':'NicheDE','sample':'patient4' if ident==24803571 else 'patient5','name':f['name'],'url':f['download_url'],'relative_path':'inputs/NicheDE/'+f['name'],'expected_bytes':f['size']})
config={'run_id':rid,'remote_run':remote,'remote_spatial_root':remote.rsplit('/background_validation/',1)[0],'genes':['CGNL1','CREB3L1','PCSK5','SLC9B2'],'V01':{'baseline_only':True,'neighbour_distance_um':100,'distance_tolerance_um':1e-4,'require_six_annotated_neighbours':True,'primary_min_tumour_neighbours':5,'sensitivity_min_tumour_neighbours':6,'near_um':[100,300],'near_lower_open':True,'far_min_um':600,'min_spots_per_group':20,'min_patients_inference':10,'single_primary':'CGNL1:5_of_6:pseudobulk','score_covariates':['hepatic_score','stromal_score','epithelial_score','endothelial_score','log1pUMI'],'families':{'secondary_genes':3,'strict6_sensitivity':4,'score_sensitivity':16}},'data_stage':{'max_download_threads':3,'total_byte_budget':12*1024**3,'deadline_hours':8,'retry_count':3,'per_file_timeout_seconds':1800,'proxy':'http://127.0.0.1:17901'},'scope':'V01 full internal robustness; V02/V04 input and matrix-QC preparation. V03/V05 wait for reference/cell-identity eligibility review; not falsely marked running.'}
(R/'CONFIG.json').write_text(json.dumps(config,indent=2));(R/'SAMPLES.json').write_text(json.dumps(samples,ensure_ascii=False,indent=2),encoding='utf-8');(R/'DOWNLOADS.json').write_text(json.dumps(files,indent=2));(B/'LATEST_RUN.json').write_text(json.dumps({'run_id':rid,'local_run':str(R),'remote_run':remote},indent=2));shutil.copy2(B/'后台验证计划.md',R/'APPROVED_PLAN.md');print(json.dumps({'run_id':rid,'samples':len(samples),'input_files':len(files),'remote':remote},ensure_ascii=False))
