from pathlib import Path
import json,gzip,hashlib,inspect,urllib.request,urllib.parse
O=Path.cwd()/'mechanism_investigation_v2';r={}
r['urllib_ftp_open_source']=inspect.getsource(urllib.request.FTPHandler.ftp_open)
r['urllib_connect_ftp_source']=inspect.getsource(urllib.request.FTPHandler.connect_ftp)
for p in sorted((O/'raw').glob('binary_retry*')):
 b=p.read_bytes();v={'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'CRLF':b.count(b'\r\n'),'LF':b.count(b'\n'),'head':b[:20].hex()}
 try:gzip.decompress(b);v['gzip_crc']=True
 except Exception as e:v['gzip_error']=repr(e)
 old=O/'raw'/p.name.removeprefix('binary_retry_')
 if old.exists():
  a=old.read_bytes();v['old_bytes']=len(a);v['identical_prefix_bytes']=next((i for i,(x,y) in enumerate(zip(a,b)) if x!=y),min(len(a),len(b)));v['old_CRLF']=a.count(b'\r\n')
 r[p.name]=v;print(p.name,json.dumps(v))
(O/'metadata/TRANSPORT_DIAGNOSIS_49.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(r['urllib_ftp_open_source']);print(r['urllib_connect_ftp_source'])
