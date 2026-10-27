from pathlib import Path
import json,hashlib,datetime,zipfile
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/increment_context_20261012_v1')
qa=json.loads((R/'FINAL_NUMERICAL_QA.json').read_text());assert qa['all_pass']
rep=json.loads((R/'REPORT_VALIDATION.json').read_text());assert rep['relative_links_valid']
statuses=[]
for p in R.glob('PROCESS_STATUS*.json'):
 x=json.loads(p.read_text());statuses.append((x.get('completed_utc',''),x))
latest={}
for _,rec in sorted(statuses):
 latest.update(rec['stages'])
assert all(latest.get(k)==0 for k in ['02_analyze.py','03_follow_states.py','04_qa.py','05_marker_conditioning.py','06_report.py']),latest
for p in (R/'scripts').glob('*.py'):compile(p.read_text(),str(p),'exec')
validation=dict(verified_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),final_exit_codes=latest,all_stages_success=True,numerical_QA=qa,report=rep,scripts_compile=True,prepare_reconstruction_pass=json.loads((R/'PREPARE_QA.json').read_text())['all_pass'])
(R/'DELIVERY_VALIDATION.json').write_text(json.dumps(validation,indent=2))
decision=json.loads((R/'INVESTIGATION_DECISION.json').read_text())
(R/'PROGRESS_SNAPSHOT.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2))
(R.parent/'CGNL1_PROGRESS_20261012_INCREMENT_CONTEXT.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2))
entries=[]
for p in sorted(R.rglob('*')):
 if p.is_file() and p.name not in ['SHA256SUMS.txt','CGNL1_increment_context_bundle.zip','BUNDLE_RECEIPT.json']:
  entries.append((p,hashlib.sha256(p.read_bytes()).hexdigest()))
(R/'SHA256SUMS.txt').write_text(''.join(h+'  '+p.relative_to(R).as_posix()+'\n' for p,h in entries))
bundle=R/'CGNL1_increment_context_bundle.zip'
with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as z:
 for p,_ in entries:z.write(p,p.relative_to(R).as_posix())
 z.write(R/'SHA256SUMS.txt','SHA256SUMS.txt')
receipt=dict(path=str(bundle),bytes=bundle.stat().st_size,sha256=hashlib.sha256(bundle.read_bytes()).hexdigest(),files=len(entries)+1)
(R/'BUNDLE_RECEIPT.json').write_text(json.dumps(receipt,indent=2))
print(json.dumps(receipt))
