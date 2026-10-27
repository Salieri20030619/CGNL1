from pathlib import Path
import os,json,time,signal,datetime
O=Path.cwd()/'mechanism_investigation_v2'
def processes():
 rows=[]
 for p in Path('/proc').iterdir():
  if p.name.isdigit():
   try:rows.append((int(p.name),(p/'cmdline').read_bytes().split(b'\0')))
   except (FileNotFoundError,PermissionError):pass
 return rows
def assert_no_active_quant():
 for pid,args in processes():
  if b'quant' in args and any(a==str(O/'software/salmon_1_10_3/bin/salmon').encode() for a in args):raise RuntimeError('Active Salmon run '+str(pid)+'; do not interrupt, wait for completion')
assert_no_active_quant();assert not (O/'results/DRA_TARGET_ANALYSIS_COMPLETE.json').exists();assert not (O/'results/DRA_TRANSCRIPT_SENSITIVITY_COMPLETE.json').exists();assert not (O/'metadata/DRA_QUANTIFICATION_COMPLETE.json').exists();before={}
for name in ['DRA_QUANTIFICATION_SPECIFICATION_v2.json','DRA_TARGET_STATISTICAL_SPECIFICATION.json','DRA_TRANSCRIPT_SENSITIVITY_SPECIFICATION.json']:
 import hashlib
 before[name]=hashlib.sha256((O/name).read_bytes()).hexdigest()
stopped=[]
for name in ['cgnl_mech91b_dra_target_when_ready.py','cgnl_mech94_transcript_sensitivity_when_ready.py','cgnl_mech46b_quantify_when_ready.py']:
 for pid,args in processes():
  if str(O.parent/'scripts'/name).encode() in args:os.kill(pid,signal.SIGTERM);stopped.append({'pid':pid,'script':name})
for row in stopped:
 for _ in range(100):
  if not (Path('/proc')/str(row['pid'])).exists():break
  time.sleep(.1)
 else:raise RuntimeError('Prior controller did not exit')
assert_no_active_quant()
r={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'localwholeFASTQpairs alreadyofficialMD5verified, butoldcontroller accepts onlyoriginaldownloadreceipts; remove route-specific scheduling dependency','stopped':stopped,'active_quant_interrupted':False,'frozen_spec_sha256':before,'new_quantifier':'cgnl_mech46c_quantify_available_files.py','new_downstream':['cgnl_mech91c_dra_target_when_ready.py','cgnl_mech94b_transcript_sensitivity_when_ready.py'],'termination':'newquantifier waitsall12runs; ifall input producers become terminal andmissingfiles remain, analyst writes metadata/DRA_INPUTS_TERMINAL.json after verification; original downloader disappearance alone notterminal becauselocalproducers exist'}
(O/'metadata/QUANT_INPUT_UPGRADE_158.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r))
