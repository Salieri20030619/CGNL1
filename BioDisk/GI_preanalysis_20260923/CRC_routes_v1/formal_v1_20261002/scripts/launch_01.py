from pathlib import Path
import subprocess,json,os,datetime
r=Path(__file__).resolve().parents[1]
env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='2'
p=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(r/'scripts/01_qualify_inputs.py')],env=env)
(r/'PROCESS_EXIT.json').write_text(json.dumps({'exit_code':p.returncode,'finished':datetime.datetime.now().astimezone().isoformat()},indent=2))
raise SystemExit(p.returncode)
