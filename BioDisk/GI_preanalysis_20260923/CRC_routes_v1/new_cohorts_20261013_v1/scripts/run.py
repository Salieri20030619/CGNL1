from pathlib import Path
import subprocess,json,datetime,sys
R=Path(__file__).resolve().parent.parent
name=sys.argv[1] if len(sys.argv)>1 else '03_extract_existing.py'
status=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=name)
p=R/'logs'/(name+'.log')
if p.exists():p.rename(p.with_name(p.name+'.'+datetime.datetime.now().strftime('%H%M%S')+'.previous'))
with p.open('w') as f:r=subprocess.run([sys.executable,'-u',str(R/'scripts'/name)],stdout=f,stderr=subprocess.STDOUT)
status.update(exit_code=r.returncode,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(R/'logs'/(name+'.exit.json')).write_text(json.dumps(status,indent=2));print(json.dumps(status),flush=True);sys.exit(r.returncode)
