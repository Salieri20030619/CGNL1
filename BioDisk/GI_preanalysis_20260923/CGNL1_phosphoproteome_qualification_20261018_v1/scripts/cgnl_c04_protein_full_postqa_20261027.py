from pathlib import Path
import os,select,json,subprocess,hashlib,datetime,traceback,fcntl,sys
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
W=P/'continuous_exploration_20261023_v1/C04_totalprotein_full_20261026';R=W/'runs/recovery_20261027_v3';Q=W/'postqueue_QA_20261027_v1';Q.mkdir(exist_ok=True)
lock=(Q/'runner.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
expected=sys.argv[1];script=P/'scripts/cgnl_c04_protein_full_q_independent_20261027.R'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def save(name,x):
 x.update(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),whole_project_complete=False)
 (Q/name).write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
try:
 assert not (Q/'EXIT.json').exists() and not (Q/'FAILURE.json').exists()
 assert sha(script)==expected
 try:fd=os.pidfd_open(2150334)
 except ProcessLookupError:assert (R/'EXIT.json').exists()
 else:
  assert 'cgnl_c04_protein_full_recovery_v3_20261027.py' in Path('/proc/2150334/cmdline').read_bytes().decode()
  save('STATUS.json',dict(stage='native_wait_for_fullscope_queue',pid=os.getpid(),queue_pid=2150334,R_script_sha256=expected));select.select([fd],[],[]);os.close(fd)
 terminal=json.loads((R/'EXIT.json').read_text());assert terminal['complete'] and not terminal.get('error')
 manifest=json.loads((W/'joint_identification/INPUT_MANIFEST.json').read_text());assert manifest['complete24scope'] and len(manifest['inputs'])==24
 for x in manifest['inputs']:assert sha(Path(x['path']))==x['sha256']
 assert sha(script)==expected;save('STATUS.json',dict(stage='running_independent_R_QA',pid=os.getpid()))
 r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/R_4.5.2/bin/Rscript','--vanilla',str(script),str(W)],capture_output=True,text=True)
 save('R_EXIT.json',dict(exit=r.returncode,stdout=r.stdout,stderr=r.stderr));assert r.returncode==0
 qa=json.loads((W/'joint_identification/INDEPENDENT_R_PSM_QA.json').read_text());assert qa['pass']
 save('EXIT.json',dict(independent_PSM_QA_complete=True,qa=qa,PSM_table_sha256=sha(W/'joint_identification/PSM_QVALUES.tsv.gz'),protein_quantification_pending=True));save('STATUS.json',dict(stage='PSM_QA_complete_protein_quantification_pending',pid=os.getpid()))
except Exception:
 save('FAILURE.json',dict(error=traceback.format_exc(),independent_PSM_QA_complete=False));raise
