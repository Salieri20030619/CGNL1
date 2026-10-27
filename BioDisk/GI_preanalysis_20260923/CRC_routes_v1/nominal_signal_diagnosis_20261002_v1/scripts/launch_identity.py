from pathlib import Path
import subprocess,json,os
r=Path(__file__).resolve().parents[1];env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='2'
x=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(r/'scripts/identity_audit.py')],env=env)
(r/'IDENTITY_EXIT.json').write_text(json.dumps({'exit_code':x.returncode}))
raise SystemExit(x.returncode)
