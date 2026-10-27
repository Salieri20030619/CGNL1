import subprocess,json,os,pathlib,datetime
p=pathlib.Path(__file__).resolve().parents[1]
a=pathlib.Path(__file__).resolve().parent
env={k:v for k,v in os.environ.items() if not k.lower().endswith('proxy')}
r=subprocess.run(['python3',str(p/'launcher_direct.py')],env=env)
errors=[]
try:
 c=json.loads((p/'COMPLETE.json').read_text()); e=json.loads((p/'DIRECT_PROCESS_EXIT.json').read_text()); m=json.loads((p/'MANIFEST.json').read_text())
 if r.returncode!=0 or c.get('exit_code')!=0 or e.get('exit_code')!=0: errors.append('nonzero exit or completion status')
 if len(m)!=304 or c.get('verified')!=304: errors.append('not all 304 files verified')
 for x in m:
  f=pathlib.Path(x['path'])
  if x.get('status') not in ['downloaded_verified','verified_existing'] or not f.is_file() or f.stat().st_size!=x['bytes']: errors.append(x['path'])
except Exception as exc: errors.append(str(exc))
result={'status':'verified_complete' if not errors else 'needs_attention','launcher_exit':r.returncode,'errors':errors,'finished':datetime.datetime.now().astimezone().isoformat()}
(a/'RETRY_COMPLETION.json').write_text(json.dumps(result,indent=2))
raise SystemExit(0 if not errors else 1)
