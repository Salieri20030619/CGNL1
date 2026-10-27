
from pathlib import Path
import subprocess,json,datetime,os
r=Path(__file__).resolve().parents[1];env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:env[k]='4'
x=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(r/'scripts/identity.py')],env=env)
(r/'IDENTITY_EXIT.json').write_text(json.dumps({'exit_code':x.returncode,'finished':datetime.datetime.now().astimezone().isoformat()}));raise SystemExit(x.returncode)
