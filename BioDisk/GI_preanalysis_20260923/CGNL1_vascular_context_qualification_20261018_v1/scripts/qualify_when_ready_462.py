from pathlib import Path
import json,os,select,subprocess,sys
b=Path(__file__).resolve().parents[1];d=json.loads((b/'metadata/PILOT455_LAUNCH.json').read_text());receipt=Path(d['receipt'])
if not receipt.exists():
 try:
  fd=os.pidfd_open(d['wrapper_pid']);p=select.poll();p.register(fd,select.POLLIN);p.poll();os.close(fd)
 except ProcessLookupError:pass
assert receipt.exists() and json.loads(receipt.read_text())['exit_code']==0,'Quantification not successful'
r=subprocess.run(['/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python',str(b/'scripts/qualify_pilot_cells_457.py')]);sys.exit(r.returncode)
