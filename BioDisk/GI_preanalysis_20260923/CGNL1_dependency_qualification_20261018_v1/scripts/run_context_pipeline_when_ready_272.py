from pathlib import Path
import json,subprocess,datetime,hashlib
p=Path.cwd();d=json.loads((p/'metadata/TARGET_RNA_LAUNCH_269.json').read_text())
subprocess.run(['/usr/bin/python3',str(p.parent/'CGNL1_gastric_peritoneal_validation_20261016_v1/scripts/cgnl_wait_pidfd_system.py'),str(d['wrapper_pid'])],check=True)
assert json.loads(Path(d['receipt']).read_text())['exit_code']==0
assert json.loads((p/'metadata/TARGET_RNA_COMPLETE_269.json').read_text())['all_requests_pass']
for name in ['analyze_organoid_context_270.py','independent_organoid_context_qa_271.py']:
 f=p/'scripts'/name;start=datetime.datetime.now(datetime.timezone.utc).isoformat();log=p/'logs'/(f.stem+'.log')
 with log.open('xb') as output:r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(f)],cwd=p,stdout=output,stderr=subprocess.STDOUT)
 receipt={'script':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'start_utc':start,'end_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'exit_code':r.returncode,'log':str(log)}
 (p/'logs'/(f.stem+'.exit.json')).write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt),flush=True);assert r.returncode==0,receipt
assert json.loads((p/'results/context_270/INDEPENDENT_QA_271.json').read_text())['pass']
(p/'metadata/CONTEXT_PIPELINE_COMPLETE_272.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'analysis_and_QA_pass':True,'scientific_interpretation_pending':True,'whole_project_complete':False}),encoding='utf-8')
