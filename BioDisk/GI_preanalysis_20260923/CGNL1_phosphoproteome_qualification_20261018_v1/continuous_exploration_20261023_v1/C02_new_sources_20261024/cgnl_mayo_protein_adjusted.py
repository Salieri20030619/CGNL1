from pathlib import Path
import json,csv,gzip,hashlib,math,statistics,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1/C02_new_sources_20261024');O=P/'protein_adjusted_models';O.mkdir(exist_ok=True)
spec=P/'PROTEIN_ADJUSTED_SPEC_BEFORE_RESULTS.json';assert hashlib.sha256(spec.read_bytes()).hexdigest()=='6f07ec2574b14328837d52861fc88b459d32c7d80a1f71c5e65b274e264f745d'
def read(f):return json.loads(f.read_text())
configs={3:[('iCCA27',10,'processed1','tumor'),('iCCA29',12,'processed2','tumor'),('iCCA30',13,'processed3','tumor'),('Normal5',14,'normal5','normal'),('Normal6',15,'normal6','normal')],4:[('dCCA1',1,'processed1','anchor'),('pCCA13',2,'processed2','other'),('iCCA32',3,'processed3','tumor'),('dCCA2',10,'processed4','other'),('Normal7',13,'normal7','normal'),('Normal8',14,'normal8','normal'),('Normal9',15,'normal9','normal')]}
calibration=[];targets=[];factors=[];results=[]
for batch,cfg in configs.items():
 pm=P/f'Set{batch:02}_protein_mapping';ph=P/'Set03_mapping' if batch==3 else P
 assert read(pm/'INDEPENDENT_QA.json')['pass'] and read(pm/'RESULTS.json')['accepted']
 assert read(ph/'CHANNEL_INDEPENDENT_QA.json')['pass'] and read(ph/'NORMAL_MAPPING_INDEPENDENT_QA.json')['pass']
 gate=read(P/f'Set{batch:02}_EVIDENCE_GATE.json');assert gate['pass']
 gp=P/f'Set{batch:02}_GP/proteinGroups.txt';rec=read(gp.with_name(gp.name+'.receipt.json'));assert hashlib.sha256(gp.read_bytes()).hexdigest()==rec['sha256']
 protein=read(P/f'Set{batch:02}_GP_CGNL1_ROWS.json')['rows'];assert len(protein)==1;protein=protein[0]
 assert protein['Protein IDs']=='Q0VF96' and not any(protein[k] for k in ['Reverse','Potential contaminant','Only identified by site'])
 with gzip.open(pm/'FEATURES.csv.gz','rt') as h:gprows=list(csv.DictReader(h))
 with gzip.open(ph/'NORMAL_MAPPING_FEATURES.csv.gz','rt') as h:phrows=list(csv.DictReader(h))
 for layer,rows in [('protein',gprows),('phospho',phrows)]:
  for sample,ch,col,group in cfg:
   values=[]
   for x in rows:
    assert 'Q0VF96' not in x['feature']
    if layer=='protein':raw=math.log2(float(x[f'raw{ch}']))-math.log2(float(x[f'raw{cfg[0][1]}']));export=math.log2(float(x[sample]))-math.log2(float(x[cfg[0][0]]))
    else:
     if not x[col]:continue
     raw=float(x[f'raw{ch}'])-float(x[f'raw{cfg[0][1]}']);export=float(x[col])-float(x[cfg[0][2]])
    if not(math.isfinite(raw) and math.isfinite(export)):continue
    item={'set':batch,'layer':layer,'sample':sample,'feature':x['feature'],'split':int(x['split']),'raw_relative_log':raw,'export_relative_log':export};values.append(item);calibration.append(item)
   train=[x for x in values if x['split']==0];hold=[x for x in values if x['split']==1];assert min(len(train),len(hold))>=500
   B=statistics.median(x['export_relative_log']-x['raw_relative_log'] for x in train);C=-statistics.median(x['raw_relative_log'] for x in train);validation=statistics.median(x['export_relative_log']-x['raw_relative_log']-B for x in hold)
   factors.append({'set':batch,'layer':layer,'sample':sample,'A':0.,'B':B,'C':C,'heldout_median':validation,'B_pass':abs(validation)<=.01})
 for sample,ch,col,group in cfg:
  a=gate['sites']['486']['site_row_mono_reporters'][ch-1];b=float(protein[f'Reporter intensity corrected {ch}']);assert a>0 and b>0 and math.isfinite(a) and math.isfinite(b)
  targets.append({'set':batch,'sample':sample,'group':group,'channel':ch,'S486':a,'CGNL1_protein':b,'raw_log_ratio':math.log2(a/b)})
 for model in ['A','B','C']:
  relevant=[x for x in factors if x['set']==batch and x['sample'] in [s for s,c,col,g in cfg if g in ['tumor','normal']]]
  if model=='B' and not all(x['B_pass'] for x in relevant):results.append({'set':batch,'model':model,'evaluable':False,'reason':'heldout normalization gate failed'});continue
  normalized={}
  for row in [x for x in targets if x['set']==batch]:
   fp=next(x[model] for x in factors if x['set']==batch and x['layer']=='phospho' and x['sample']==row['sample']);fg=next(x[model] for x in factors if x['set']==batch and x['layer']=='protein' and x['sample']==row['sample']);normalized[row['sample']]=row['raw_log_ratio']+fp-fg
  for tumor in [x for x in targets if x['set']==batch and x['group']=='tumor']:
   dd={n['sample']:normalized[tumor['sample']]-normalized[n['sample']] for n in targets if n['set']==batch and n['group']=='normal'}
   results.append({'set':batch,'model':model,'evaluable':True,'tumor':tumor['sample'],'differences':dd,'median':statistics.median(dd.values()),'range':[min(dd.values()),max(dd.values())]})
with gzip.open(O/'CALIBRATION_INPUT.csv.gz','wt',encoding='utf-8',newline='') as h:
 w=csv.DictWriter(h,fieldnames=list(calibration[0]));w.writeheader();w.writerows(calibration)
for name,rows in [('TARGET_INPUT.csv',targets),('FACTORS.csv',factors)]:
 with (O/name).open('w',encoding='utf-8',newline='') as h:w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
out={'spec_sha256':hashlib.sha256(spec.read_bytes()).hexdigest(),'results':results,'R_QA_pending':True,'new_biological_tests':0,'P':None,'FDR':None,'whole_project_complete':False};(O/'RESULTS.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(results))
