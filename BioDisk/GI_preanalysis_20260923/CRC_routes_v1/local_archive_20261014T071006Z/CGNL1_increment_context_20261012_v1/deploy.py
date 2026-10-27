from pathlib import Path
import subprocess,sys
L=Path(__file__).resolve().parent
R='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/increment_context_20261012_v1'
PY='/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006/env/methods/bin/python'
def run(args):
 return subprocess.run(args,check=True,creationflags=subprocess.CREATE_NO_WINDOW)
run(['ssh','-T','-o','BatchMode=yes','sgs10','mkdir -p '+R+'/scripts'])
for name in sys.argv[1:]:
 run(['scp','-q',str(L/name),'sgs10:'+R+'/scripts/'+name])
run(['ssh','-T','-o','BatchMode=yes','sgs10',PY+' '+R+'/scripts/'+sys.argv[-1]])
