from pathlib import Path
import requests,json,time
O=Path.cwd()/'mechanism_investigation_v2';plan=json.loads((O/'metadata/DRA_DOWNLOAD_PLAN.json').read_text());t=plan['tasks'][0];start=time.monotonic();r={'url':t['url'],'request':'bytes=0-65535','reason':'assessresumeandtransfercapacitywithoutduplicatinglargedownload'}
try:
 with requests.get(t['url'],headers={'Range':'bytes=0-65535','Accept-Encoding':'identity'},stream=True,timeout=(15,30)) as a:
  b=a.raw.read(65536);r.update(status=a.status_code,content_range=a.headers.get('Content-Range'),bytes_read=len(b),elapsed_seconds=time.monotonic()-start)
except Exception as e:r['error']=repr(e)
(O/'metadata/ENA_RANGE_PROBE.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r))
