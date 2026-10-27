import pathlib,sys,json,subprocess,os,fcntl
R=pathlib.Path(sys.argv[1]);D=R/'V01';D.mkdir(exist_ok=True)
if (D/'STARTED.json').exists():
 j=json.loads((D/'STARTED.json').read_text());pid=j['pid'];p=pathlib.Path('/proc')/str(pid)/'cmdline'
 if p.exists() and str(R).encode() in p.read_bytes():print(json.dumps({'already_running':pid}));sys.exit(0)
 if (D/'COMPLETE.json').exists():print(json.dumps({'already_completed':True}));sys.exit(0)
with (D/'entry.log').open('ab') as log:p=subprocess.Popen([sys.executable,str(R/'scripts/task_entry.py'),str(R),'V01'],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
(R/'V01_PID.json').write_text(json.dumps({'pid':p.pid,'run':str(R)},indent=2));print(json.dumps({'launched_V01_pid':p.pid}),flush=True)
