import psutil,json,datetime
from pathlib import Path
L=Path(r'C:\Users\31937\CGNL1_new_cohorts_20261013_v1');stopped=[]
for p in psutil.process_iter(['pid','name','cmdline']):
 try:
  args=p.info['cmdline'] or []
  if p.info['name'].lower().startswith('python') and any(str(a).endswith('CGNL1_new_cohorts_20261013_v1\\download_ranges.py') for a in args):
   children=p.children(recursive=True)
   for c in children:c.kill()
   p.kill();stopped.append(dict(pid=p.pid,children=[c.pid for c in children]))
 except (psutil.NoSuchProcess,psutil.AccessDenied):pass
(L/'HTTP_TRANSFER_INTERRUPTION.json').write_text(json.dumps(dict(stopped=stopped,reason='Official range transfer timed out after180s with only455KB. Preserve parts; use configured km112 route which passed HTTPS probe.'),indent=2))
print(json.dumps(stopped))