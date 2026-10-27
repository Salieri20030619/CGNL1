import requests,concurrent.futures,json,socket
from pathlib import Path
L=Path(__file__).resolve().parent
def test(port):
 try:
  with socket.create_connection(('127.0.0.1',port),timeout=2):pass
 except Exception as e:return dict(port=port,status='not_listening')
 try:
  proxy={'https':'http://127.0.0.1:'+str(port),'http':'http://127.0.0.1:'+str(port)}
  r=requests.get('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200997/suppl/',proxies=proxy,timeout=(8,20));r.raise_for_status()
  (L/'metadata'/('GSE200997_index_proxy'+str(port)+'.txt')).write_text(r.text)
  return dict(port=port,status=r.status_code,bytes=len(r.content),preview=r.text[:2200])
 except Exception as e:return dict(port=port,error=repr(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as p:
 for x in p.map(test,[7897,9876,7890]):print(json.dumps(x),flush=True)
