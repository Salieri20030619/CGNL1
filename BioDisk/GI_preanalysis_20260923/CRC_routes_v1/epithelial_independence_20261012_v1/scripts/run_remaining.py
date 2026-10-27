from pathlib import Path
import os,subprocess,json,datetime,select,ctypes,traceback
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/epithelial_independence_20261012_v1');P='/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006/env/methods/bin/python'
previous=json.loads((R/'PROCESS_STATUS.json').read_text())
(R/'PROCESS_STATUS_first_waiter.json').write_text(json.dumps(previous,indent=2))
status={'AB':int((R/'logs/AB.exit').read_text()),'waiter_compatibility_fix':'libc pidfd_open for older Python'}
try:
 with (R/'logs/B_subtype.log').open('w') as log:rc=subprocess.call([P,str(R/'scripts/02b_verify_subtype.py')],stdout=log,stderr=subprocess.STDOUT)
 (R/'logs/B_subtype.exit').write_text(str(rc)+'\n');status['B_subtype']=rc;assert rc==0
 if not (R/'logs/extract.exit').exists():
  libc=ctypes.CDLL(None,use_errno=True)
  for d in Path('/proc').iterdir():
   if not d.name.isdigit():continue
   try:
    if str(R/'scripts/run_extract.sh').encode() not in (d/'cmdline').read_bytes().split(b'\0'):continue
   except OSError:continue
   pid=int(d.name);fd=libc.pidfd_open(pid,0)
   if fd<0:
    if ctypes.get_errno()==3:continue
    raise OSError(ctypes.get_errno(),'pidfd_open')
   print('NATIVE_WAIT',pid,flush=True);poll=select.poll();poll.register(fd,select.POLLIN);poll.poll();os.close(fd)
 assert (R/'logs/extract.exit').read_text().strip()=='0'
 assert (R/'results/EXTRACTION_COMPLETE.json').exists()
 with (R/'logs/C.log').open('w') as log:rc=subprocess.call([P,str(R/'scripts/03_analyze_C.py')],stdout=log,stderr=subprocess.STDOUT)
 (R/'logs/C.exit').write_text(str(rc)+'\n');status['C']=rc;status['extract']=0;assert rc==0
 status['all_success']=True
except Exception as e:status['all_success']=False;status['error']=repr(e);traceback.print_exc()
(R/'PROCESS_STATUS.json').write_text(json.dumps(status,indent=2))
print('ANALYSIS_COMPLETE',json.dumps(status),flush=True)
