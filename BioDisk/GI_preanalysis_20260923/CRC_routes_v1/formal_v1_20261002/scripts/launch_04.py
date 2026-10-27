from pathlib import Path
import subprocess,json,os,datetime
r=Path(__file__).resolve().parents[1];d=r/'04_patient_state';env=os.environ.copy()
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:env[k]='2'
x=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(r/'scripts/04_patient_state.py')],env=env)
(d/'PROCESS_EXIT.json').write_text(json.dumps({'exit_code':x.returncode,'finished':datetime.datetime.now().astimezone().isoformat()}))
raise SystemExit(x.returncode)
