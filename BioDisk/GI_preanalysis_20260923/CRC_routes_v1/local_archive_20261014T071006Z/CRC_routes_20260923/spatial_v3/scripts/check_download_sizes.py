from pathlib import Path
import json,urllib.request,concurrent.futures,sys
sys.stdout.reconfigure(encoding='utf-8')
O=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3');op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}));jobs=json.loads((O/'metadata/LIVER_DOWNLOAD_PLAN.json').read_text())
def f(j):
 try:
  r=op.open(urllib.request.Request(j['url'],method='HEAD'),timeout=25);return dict(j,bytes=int(r.headers.get('Content-Length',0)))
 except Exception as e:return dict(j,error=str(e))
r=list(concurrent.futures.ThreadPoolExecutor(max_workers=6).map(f,jobs));(O/'metadata/LIVER_FILE_SIZES.json').write_text(json.dumps(r,indent=2));print(json.dumps({'total_bytes_known':sum(x.get('bytes',0) for x in r),'files':r},ensure_ascii=False))
