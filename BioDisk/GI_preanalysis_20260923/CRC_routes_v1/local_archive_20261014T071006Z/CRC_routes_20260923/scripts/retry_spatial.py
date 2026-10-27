import pathlib,json,subprocess,time,hashlib,concurrent.futures
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');jobs=json.loads((O/'spatial_download_plan.json').read_text());bad=[]
for name,url in jobs:
 j=json.loads((O/'raw'/(name+'.complete.json')).read_text())
 if j['exit_code'] or j.get('transfer_exit_code'):bad.append((name,url))
def run(v):
 name,u=v;p=O/'raw'/(name+'.tar.gz');t=time.time()
 with (O/'raw'/(name+'.retry.log')).open('w') as f:r=subprocess.run(['curl.exe','-x','http://127.0.0.1:7897','--http1.1','-L','--fail','--connect-timeout','15','--max-time','180','--speed-time','20','--speed-limit','50000','-C','-','-o',str(p)+'.part',u+'?download=1'],stdout=f,stderr=f)
 j={'sample':name,'url':u,'exit_code':r.returncode,'seconds':time.time()-t}
 if r.returncode==0:
  pathlib.Path(str(p)+'.part').replace(p);j['bytes']=p.stat().st_size;j['sha256']=hashlib.file_digest(p.open('rb'),'sha256').hexdigest();j['transfer_exit_code']=subprocess.run(['scp',str(p),'sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/raw/']).returncode
 (O/'raw'/(name+'.retry.complete.json')).write_text(json.dumps(j,indent=2));print(json.dumps(j),flush=True);return j
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as e:rr=list(e.map(run,bad))
(O/'raw/SPATIAL_RETRY_COMPLETE.json').write_text(json.dumps(rr,indent=2));print('RETRY COMPLETE',flush=True)
