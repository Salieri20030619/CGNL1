from pathlib import Path
import subprocess,os,sys,json,datetime
R=Path(__file__).resolve().parent.parent
ps=subprocess.run(['ps','-u',str(os.getuid()),'-o','pid=,args='],capture_output=True,text=True,check=True).stdout
pids=[]
for line in ps.splitlines():
 pid,cmd=line.strip().split(None,1)
 if 'python' in cmd and cmd.endswith('/new_cohorts_20261013_v1/scripts/02_download.py'):pids.append(int(pid))
print('WAIT_DOWNLOADER',pids,flush=True)
for pid in pids:subprocess.run(['tail','--pid='+str(pid),'-f','/dev/null'],check=True)
assert (R/'DOWNLOAD_COMPLETE.json').exists(),'Downloader ended without completion marker'
for name in ['04_extract_primary.py']:
 p=subprocess.run([sys.executable,str(R/'scripts/run.py'),name]);assert p.returncode==0,name
print('NEW_COHORT_PREPARATION_COMPLETE',flush=True)
