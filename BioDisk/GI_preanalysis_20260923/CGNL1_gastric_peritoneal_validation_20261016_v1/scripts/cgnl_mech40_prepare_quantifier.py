from pathlib import Path
import subprocess,json,datetime,os
O=Path.cwd()/'mechanism_investigation_v2';S=O/'software';S.mkdir(exist_ok=True);env=S/'salmon_1_10_3'
cmd=['/mnt/BioDisk/wangmingyuan/miniconda3/bin/conda','create','-y','-p',str(env),'--override-channels','-c','conda-forge','-c','bioconda','salmon=1.10.3']
assert not env.exists(),'Inspect existing environment rather than overwrite'
r={'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':cmd,'scope':'isolated project environment; existing environments unchanged'}
with (O/'logs/conda_salmon_install.log').open('wb') as f:proc=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=1800)
r['exit_code']=proc.returncode
if proc.returncode==0:
 a=subprocess.run([str(env/'bin/salmon'),'--version'],capture_output=True,text=True);r.update(version=a.stdout.strip(),version_exit=a.returncode)
 a=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/bin/conda','list','-p',str(env),'--explicit'],capture_output=True);(O/'metadata/SALMON_ENV_EXPLICIT.txt').write_bytes(a.stdout)
r['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(O/'metadata/QUANTIFIER_INSTALL.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r))
if proc.returncode:raise RuntimeError('Installation failed; see recorded log')
