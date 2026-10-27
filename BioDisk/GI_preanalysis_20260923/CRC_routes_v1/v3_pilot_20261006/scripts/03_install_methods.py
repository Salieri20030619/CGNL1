from pathlib import Path
import os,subprocess,json,datetime,time
ROOT=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/v3_pilot_20261006');BIG=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
env=os.environ.copy()
for k in list(env):
 if 'proxy' in k.lower():env.pop(k,None)
env.update(NO_PROXY='*',no_proxy='*',PIP_CONFIG_FILE='/dev/null',PIP_DISABLE_PIP_VERSION_CHECK='1',OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8')
base='/mnt/BioDisk/wangmingyuan/miniconda3/envs/cell2loc/bin/python';v=BIG/'env/methods'
subprocess.run([base,'-m','venv','--system-site-packages',str(v)],env=env,check=True)
py=str(v/'bin/python')
cmd=[py,'-m','pip','--isolated','install','--index-url','https://pypi.org/simple','--proxy','','--timeout','30','--retries','1','--no-deps','--no-cache-dir','nichecompass==0.3.3','torch-geometric','liana']
t=time.time()
with (ROOT/'logs/03_pip.log').open('w') as log:
 try:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600);code=r.returncode
 except subprocess.TimeoutExpired:code=124
mods={}
for m in ['torch','torch_geometric','nichecompass','liana']:
 r=subprocess.run([py,'-c',f'import {m};print(getattr({m},"__version__","unknown"))'],env=env,capture_output=True,text=True,timeout=90)
 mods[m]={'exit_code':r.returncode,'stdout':r.stdout[-300:],'stderr':r.stderr[-2000:]}
out={'pip_exit_code':code,'elapsed_s':time.time()-t,'python':py,'network_mode':'direct_only','imports':mods}
(ROOT/'results/ENVIRONMENT_INSTALL.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
