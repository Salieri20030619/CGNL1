from pathlib import Path
import os,signal,subprocess,json
pid=2686413
p=Path('/proc')/str(pid)/'cmdline'
if p.exists():
 cmd=p.read_bytes().replace(b'\0',b' ').decode();assert '/spatial_v3/scripts/download_liver.py' in cmd
 os.kill(pid,signal.SIGSTOP)
 children=subprocess.run(['pgrep','-P',str(pid)],capture_output=True,text=True).stdout.split()
 for child in children:
  try:os.kill(int(child),signal.SIGTERM)
  except ProcessLookupError:pass
 os.kill(pid,signal.SIGKILL);print(json.dumps({'stopped_pid':pid,'children':children}))
