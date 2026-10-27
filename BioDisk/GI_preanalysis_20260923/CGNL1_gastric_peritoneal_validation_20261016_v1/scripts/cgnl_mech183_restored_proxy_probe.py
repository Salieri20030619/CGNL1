from pathlib import Path
import requests,subprocess,json,time,hashlib,datetime
O=Path.cwd()/'mechanism_investigation_v2';listeners=subprocess.run(['ss','-ltnp'],capture_output=True,text=True).stdout;hits=[l for l in listeners.splitlines() if ':7897' in l];out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'listener':hits,'purpose':'read-only bounded publicENA transfer via existing loopbackproxy; no tunnel or globalsettings changed'}
if hits:
 task=next(t for t in json.loads((O/'metadata/DRA_DOWNLOAD_PLAN.json').read_text(encoding='utf-8'))['tasks'] if t['name']=='DRR400690_1.fastq.gz');p=O/'raw/DRA014736'/task['name']
 with p.open('rb') as f:expected=f.read(1048576)
 start=time.monotonic();body=bytearray();session=requests.Session();session.trust_env=False
 try:
  with session.get(task['url'],proxies={'https':'http://127.0.0.1:7897'},headers={'Range':'bytes=0-1048575','Accept-Encoding':'identity'},stream=True,timeout=(15,30)) as a:
   out.update(status=a.status_code,content_range=a.headers.get('Content-Range'));assert a.status_code==206 and a.headers['Content-Range']==f'bytes 0-1048575/{task["bytes"]}'
   for b in a.iter_content(65536):
    body.extend(b)
    if len(body)>1048576 or time.monotonic()-start>45:raise TimeoutError('bounded public probe')
  assert bytes(body)==expected;out.update(success=True,identical_to_verified_file_prefix=True)
 except Exception as e:out.update(success=False,error=repr(e))
 out.update(bytes=len(body),seconds=time.monotonic()-start,sha256=hashlib.sha256(body).hexdigest());out['bytes_per_second']=len(body)/out['seconds'];folder=O/'raw/transport_probes';folder.mkdir(exist_ok=True);(folder/'ENA_proxy_183.bin').write_bytes(body)
(O/'metadata/ENA_PROXY_PROBE_183.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out))
