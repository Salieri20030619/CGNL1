from pathlib import Path
import subprocess,sys
L=Path(__file__).resolve().parent
R='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1'
PY='/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006/env/methods/bin/python'
def run(a,retries=1):
 for i in range(retries):
  p=subprocess.run(a,capture_output=True,text=True,creationflags=subprocess.CREATE_NO_WINDOW)
  if p.returncode==0:
   if p.stdout:print(p.stdout,flush=True)
   return
  print('COMMAND_FAILED',a[0],'attempt',i+1,p.returncode,p.stderr[-1500:],flush=True)
 p.check_returncode()
opts=['-o','BatchMode=yes','-o','ConnectTimeout=20']
run(['ssh','-T']+opts+['sgs10','mkdir -p '+R+'/scripts'],3)
for name in sys.argv[1:]:run(['scp','-q']+opts+[str(L/name),'sgs10:'+R+'/scripts/'+name],3)
run(['ssh','-T']+opts+['sgs10',PY+' -u '+R+'/scripts/'+sys.argv[-1]])
