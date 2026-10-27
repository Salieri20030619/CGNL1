from pathlib import Path
import subprocess,json,datetime
r=Path(__file__).resolve().parents[1]
x=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(r/'scripts/association.py')])
(r/'ASSOCIATION_EXIT.json').write_text(json.dumps({'exit_code':x.returncode,'time':datetime.datetime.now().astimezone().isoformat()}))
raise SystemExit(x.returncode)
