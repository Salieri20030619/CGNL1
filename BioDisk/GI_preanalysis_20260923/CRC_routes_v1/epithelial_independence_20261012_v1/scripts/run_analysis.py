from pathlib import Path
import os,subprocess,json,datetime,select,traceback
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/epithelial_independence_20261012_v1');P='/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006/env/methods/bin/python'
status={}
def run(name,script):
 print('START',name,datetime.datetime.now().isoformat(),flush=True)
 with (R/'logs'/f'{name}.log').open('w') as log: rc=subprocess.call([P,str(R/'scripts'/script)],cwd=R,stdout=log,stderr=subprocess.STDOUT)
 (R/'logs'/f'{name}.exit').write_text(str(rc)+'\n');status[name]=rc
 print('EXIT',name,rc,flush=True)
 if rc:raise RuntimeError(name+' failed')
try:
 assert (R/'results/B_EXTRACTION_COMPLETE.json').exists()
 run('AB','02_analyze_AB.py')
 if not (R/'logs/extract.exit').exists():
  matches=[]
  for p in Path('/proc').iterdir():
   if not p.name.isdigit():continue
   try:
    cmd=(p/'cmdline').read_bytes().split(b'\0')
    if str(R/'scripts/run_extract.sh').encode() in cmd:matches.append(int(p.name))
   except (OSError,PermissionError):pass
  if not matches:raise RuntimeError('No extraction exit marker or native waiter target')
  print('NATIVE_WAIT_EXTRACTION_WRAPPER',matches,flush=True)
  for pid in matches:
   try:
    fd=os.pidfd_open(pid);poll=select.poll();poll.register(fd,select.POLLIN);poll.poll();os.close(fd)
   except ProcessLookupError:pass
 assert (R/'logs/extract.exit').read_text().strip()=='0'
 assert (R/'results/EXTRACTION_COMPLETE.json').exists()
 run('C','03_analyze_C.py')
 status['all_success']=True
except Exception as e:
 status['all_success']=False;status['error']=repr(e);traceback.print_exc()
finally:
 (R/'PROCESS_STATUS.json').write_text(json.dumps(status,indent=2))
 print('ANALYSIS_COORDINATOR_COMPLETE',json.dumps(status),flush=True)
