
from pathlib import Path
import json,subprocess,zipfile,csv,datetime
p=Path.cwd();old=p.parent/'CGNL1_gastric_peritoneal_validation_20261016_v1'
launch=json.loads((p/'metadata/DOWNLOAD_LAUNCH_215.json').read_text(encoding='utf-8'))
if not Path(launch['receipt']).exists():subprocess.run(['/usr/bin/python3',str(old/'scripts/cgnl_wait_pidfd_system.py'),str(launch['wrapper_pid'])],check=True)
assert Path(launch['receipt']).exists() and json.loads(Path(launch['receipt']).read_text())['exit_code']==0
assert json.loads((p/'metadata/DOWNLOAD_COMPLETE_215.json').read_text())['all_pass']
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'response_effects_not_inspected':True}
with zipfile.ZipFile(p/'raw/rnaseq_20191101.zip') as z:
 result['expression_members']=[{'name':i.filename,'bytes':i.file_size} for i in z.infolist()]
 result['expression_previews']={}
 for name in z.namelist():
  if name.endswith('.csv'):
   import io
   with io.TextIOWrapper(z.open(name),encoding='utf-8') as h:
    rd=csv.reader(h);header=next(rd);rows=list(__import__('itertools').islice(rd,2))
   result['expression_previews'][name]={'header_first15':header[:15],'columns':len(header),'row_first15':[r[:15] for r in rows]}
with (p/'raw/model_list_20230608.csv').open(encoding='utf-8',newline='') as h:
 reader=csv.DictReader(h);models=list(reader);result['model_columns']=reader.fieldnames;result['model_count']=len(models);result['model_first2']=models[:2]
aliases={'paclitaxel':'paclitaxel','docetaxel':'docetaxel','vinblastine':'vinblastine','vincristine':'vincristine','cisplatin':'cisplatin','oxaliplatin':'oxaliplatin','5-fluorouracil':'fluorouracil','fluorouracil':'fluorouracil','gemcitabine':'gemcitabine'}
result['drug_coverage']={}
for name in ['GDSC1_fitted_dose_response_24Jul22.csv','GDSC2_fitted_dose_response_24Jul22.csv']:
 with (p/'raw'/name).open(encoding='utf-8',newline='') as h:
  rd=csv.DictReader(h);headers=rd.fieldnames;groups={}
  for r in rd:
   drug=aliases.get(r['DRUG_NAME'].lower())
   if drug:
    key=(drug,r['DRUG_ID']);groups.setdefault(key,set()).add(r['COSMIC_ID'])
 result['drug_coverage'][name]={'columns':headers,'drug_ID_counts':[{'drug':k[0],'drug_id':k[1],'unique_models':len(v)} for k,v in sorted(groups.items())]}
(p/'metadata/INPUT_STRUCTURE_QUALIFICATION_217.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2),flush=True)
