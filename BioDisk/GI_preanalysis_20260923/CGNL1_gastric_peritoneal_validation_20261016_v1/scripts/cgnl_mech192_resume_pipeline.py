from pathlib import Path
import subprocess,sys,json,datetime
p=Path.cwd();m=p/'mechanism_investigation_v2'
def run(name):
 f=(p/'logs'/(name+'.resume20261017.log')).open('xb')
 proc=subprocess.Popen([sys.executable,str(p/'scripts'/name)],cwd=p,stdout=f,stderr=subprocess.STDOUT)
 return proc,f
r,rf=run('cgnl_mech191_direct_recovery.py')
q,qf=run('cgnl_mech46c_quantify_available_files.py')
rc=r.wait();rf.close()
receipt=json.loads((m/'metadata/FAILED_FILE_RECOVERY_COMPLETE_191.json').read_text(encoding='utf-8')) if (m/'metadata/FAILED_FILE_RECOVERY_COMPLETE_191.json').exists() else {}
if rc or not receipt.get('all_pass'):
 (m/'metadata/RESUME_RECOVERY_FAILURE_192.json').write_text(json.dumps({'exit_code':rc,'receipt':receipt}),encoding='utf-8')
 raise RuntimeError('Recovery incomplete; quantifier remains independent, inspect before further action')
qc=q.wait();qf.close();assert qc==0
assert json.loads((m/'metadata/DRA_QUANTIFICATION_COMPLETE.json').read_text(encoding='utf-8'))['all_pass']
for name in ['cgnl_mech91c_dra_target_when_ready.py','cgnl_mech94b_transcript_sensitivity_when_ready.py']:
 c,f=run(name);status=c.wait();f.close();assert status==0,(name,status)
(m/'metadata/RESUME_PIPELINE_COMPLETE_192.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':True,'independent_QA_still_required':True}),encoding='utf-8')
