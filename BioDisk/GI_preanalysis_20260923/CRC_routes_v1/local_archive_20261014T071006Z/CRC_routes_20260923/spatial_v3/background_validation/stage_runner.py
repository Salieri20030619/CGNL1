from pathlib import Path
import json,hashlib
B=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3\background_validation');j=json.loads((B/'LATEST_QUEUE.json').read_text());Q=Path(j['local'])
for src,dst in [('queue_runner_w2.py','queue_runner.py'),('launch_queue_w2.py','launch_queue.py')]:
 s=(B/src).read_text(encoding='utf-8-sig');compile(s,src,'exec');(Q/'scripts'/dst).write_text(s,encoding='utf-8')
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Q/'scripts').glob('*.py')};(Q/'SCRIPT_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
print('Queue scripts prepared')
