from pathlib import Path
import psutil,subprocess,json,datetime,sys
H=Path(r'D:\CGNL1\handoff');W=Path(r'D:\CGNL1\_transfer')
p=psutil.Process(5100)
assert 'migrate.py' in ' '.join(p.cmdline())
scp=psutil.Process(18788)
assert scp.ppid()==5100 and 'SCP2079_counts' in ' '.join(scp.cmdline())
print('Waiting natively for the current raw-file upload before switching archival strategy.',flush=True)
exit_code=scp.wait()
assert exit_code==0,exit_code
# Only stop our own superseded coordinator after the upload has exited successfully.
if p.is_running():
 assert 'migrate.py' in ' '.join(p.cmdline())
 p.terminate();p.wait(timeout=15)
receipt={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'raw_scp_pid':18788,'raw_scp_exit_code':exit_code,'superseded_coordinator_pid':5100,'reason':'Use hash-matched server copies and transfer only missing archival files'}
(H/'TRANSFER_HANDOFF.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt),flush=True)
r=subprocess.run([sys.executable,'-u',str(W/'migrate_reuse.py')],creationflags=subprocess.CREATE_NO_WINDOW)
(H/'MIGRATION_EXIT.json').write_text(json.dumps({'exit_code':r.returncode,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2),encoding='utf-8')
sys.exit(r.returncode)