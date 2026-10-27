import pathlib,sys,json,os
R=pathlib.Path(sys.argv[1]);out={'run':R.name,'jobs':{}}
for job in ['V01','DATA_PREP']:
 D=R/job;j={}
 for name in ['STARTED.json','EXIT.json','COMPLETE.json','download_progress.json','RESULT_STATUS.json']:
  p=D/name
  if p.exists():
   try:j[name]=json.loads(p.read_text())
   except json.JSONDecodeError:j[name]='write_in_progress'
 started=j.get('STARTED.json',{});pid=started.get('pid');q=pathlib.Path('/proc')/str(pid)/'cmdline' if pid else None;j['process_matches_run']=bool(q and q.exists() and str(R).encode() in q.read_bytes());out['jobs'][job]=j
print(json.dumps(out,ensure_ascii=False,indent=2))
