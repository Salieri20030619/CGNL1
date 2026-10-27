import subprocess,json
from pathlib import Path
code='''
import sys,importlib.util,json
mods=['numpy','pandas','scipy','statsmodels','matplotlib','fast_matrix_market','scanpy']
print(json.dumps({'python':sys.executable,'modules':{m:bool(importlib.util.find_spec(m)) for m in mods}}))
'''
for py in ['/mnt/BioDisk/wangmingyuan/miniconda3/bin/python','/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python','/usr/bin/python3']:
 p=subprocess.run(['ssh','-o','BatchMode=yes','sgs10',py+' -'],input=code.encode('utf-8'),capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 print(p.stdout.decode('utf-8'),p.stderr.decode('utf-8'),flush=True)