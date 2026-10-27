from pathlib import Path
import json,hashlib,datetime,zipfile
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1')
qa=json.loads((R/'NUMERICAL_QA.json').read_text());assert qa['all_pass']
statuses={}
for name in ['03_extract_existing.py','04_extract_primary.py','05_analyze.py','06_qa.py','07_report.py']:
 rec=json.loads((R/'logs'/(name+'.exit.json')).read_text());statuses[name]=rec['exit_code'];assert rec['exit_code']==0
for p in (R/'scripts').glob('*.py'):compile(p.read_text(),str(p),'exec')
validation=dict(verified_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),final_analysis_exit_codes=statuses,numerical_QA=qa,report=json.loads((R/'REPORT_VALIDATION.json').read_text()),scripts_compile=True,all_pass=True)
(R/'DELIVERY_VALIDATION.json').write_text(json.dumps(validation,indent=2))
decision=json.loads((R/'INVESTIGATION_DECISION.json').read_text())
(R/'PROGRESS_SNAPSHOT.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2))
(R.parent/'CGNL1_PROGRESS_20261013_NEW_COHORTS.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2))
items=[]
for p in sorted(R.rglob('*')):
 if not p.is_file() or 'raw' in p.relative_to(R).parts or p.suffix=='.part' or p.name in ['SHA256SUMS.txt','BUNDLE_RECEIPT.json','CGNL1_new_cohorts_bundle.zip']:continue
 items.append((p,hashlib.sha256(p.read_bytes()).hexdigest()))
(R/'SHA256SUMS.txt').write_text(''.join(h+'  '+p.relative_to(R).as_posix()+'\n' for p,h in items))
bundle=R/'CGNL1_new_cohorts_bundle.zip'
with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
 for p,_ in items:z.write(p,p.relative_to(R).as_posix())
 z.write(R/'SHA256SUMS.txt','SHA256SUMS.txt')
receipt=dict(path=str(bundle),bytes=bundle.stat().st_size,sha256=hashlib.sha256(bundle.read_bytes()).hexdigest(),files=len(items)+1,raw_count_matrices_not_duplicated_in_delivery=True)
(R/'BUNDLE_RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
