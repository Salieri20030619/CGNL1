from pathlib import Path
import json,os,signal,datetime
O=Path.cwd()/'mechanism_investigation_v2';assert not (O/'results/DRA_TARGET_ANALYSIS_COMPLETE.json').exists();stopped=[]
for d in Path('/proc').iterdir():
 if d.name.isdigit():
  try:
   if str(O.parent/'scripts/cgnl_mech91_dra_target_when_ready.py').encode() in (d/'cmdline').read_bytes().split(b'\0'):os.kill(int(d.name),signal.SIGTERM);stopped.append(int(d.name))
  except (FileNotFoundError,PermissionError):pass
r={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason':'Beforeanyquantificationoutcomes,implement exactsignedrankwithties/zero-differences correctly; zero-variancepairedt leftunavailableinstead ofinfiniteprecision','frozen_comparisons_unchanged':True,'stopped_waiting_controller_pids':stopped,'replacement':'cgnl_mech91b_dra_target_when_ready.py'};(O/'metadata/DRA_STATISTICAL_IMPLEMENTATION_FIX.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r))
