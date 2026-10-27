import pathlib,sys,json,subprocess,os
Q=pathlib.Path(sys.argv[1]);p=Q/'PID.json'
if p.exists():
 pid=json.loads(p.read_text())['pid'];proc=pathlib.Path('/proc')/str(pid)/'cmdline'
 if proc.exists() and str(Q).encode() in proc.read_bytes():print(json.dumps({'already_running':pid}));sys.exit(0)
 if (Q/'QUEUE_EXIT.json').exists():print(json.dumps({'already_finished':True}));sys.exit(0)
with (Q/'queue.log').open('ab') as log:p=subprocess.Popen([sys.executable,str(Q/'scripts/queue_runner.py'),str(Q)],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
(Q/'LAUNCH.json').write_text(json.dumps({'pid':p.pid,'queue':str(Q),'detached':True},indent=2));print(json.dumps({'launched_queue_pid':p.pid,'queue':str(Q)}))
