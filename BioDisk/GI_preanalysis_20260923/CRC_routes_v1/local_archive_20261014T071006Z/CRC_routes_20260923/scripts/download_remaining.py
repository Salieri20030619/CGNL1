import pathlib,json,time,hashlib,subprocess,concurrent.futures,urllib.request
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');jobs=[v for v in json.loads((O/'spatial_download_plan.json').read_text()) if v[0] in ['M-ST-33','M-ST-34']]
def run(v):
 name,url=v;p=O/'raw'/(name+'.tar.gz');part=pathlib.Path(str(p)+'.part');records=[]
 opener=urllib.request.build_opener(urllib.request.ProxyHandler({'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}))
 for attempt in range(3):
  t=time.time();offset=part.stat().st_size if part.exists() else 0
  try:
   req=urllib.request.Request(url+'?download=1&attempt='+str(attempt),headers={'Range':'bytes='+str(offset)+'-'})
   with opener.open(req,timeout=25) as response:
    status=response.status;expected=response.headers.get('Content-Length');received=0
    if status==206:assert response.headers.get('Content-Range','').startswith('bytes '+str(offset)+'-')
    else:assert status==200;offset=0
    with part.open('ab' if offset else 'wb') as f:
     while True:
      b=response.read(1024*1024)
      if not b:break
      f.write(b);received+=len(b)
      elapsed=time.time()-t
      if elapsed>180 or (elapsed>25 and received/elapsed<50000):raise TimeoutError('bounded low-throughput stop')
    assert expected is None or received==int(expected)
   part.replace(p);j={'sample':name,'exit_code':0,'attempt':attempt,'seconds':time.time()-t,'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()};j['transfer_exit_code']=subprocess.run(['scp',str(p),'sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/raw/']).returncode;records.append(j);print(json.dumps(j),flush=True);break
  except Exception as e:records.append({'sample':name,'attempt':attempt,'error':repr(e),'seconds':time.time()-t});print(json.dumps(records[-1]),flush=True)
 (O/'raw'/(name+'.openssl.complete.json')).write_text(json.dumps(records,indent=2));return records
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:r=list(ex.map(run,jobs))
(O/'raw/LAST_DOWNLOAD_COMPLETE.json').write_text(json.dumps(r,indent=2))
