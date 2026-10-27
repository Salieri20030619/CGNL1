from pathlib import Path
import urllib.request,urllib.parse,os,json,requests,concurrent.futures,gzip,hashlib
O=Path.cwd()/'mechanism_investigation_v2';print('PROXY_HOSTS',{k:urllib.parse.urlparse(v).hostname for k,v in urllib.request.getproxies().items()},flush=True)
base='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223165/suppl/'
urls={'exposure_filelist.txt':base+'filelist.txt','exposure_listing.html':base,'https_GSE223165_fpkm_sample.txt.gz':base+'GSE223165_fpkm_sample.txt.gz'}
def get(item):
 n,u=item;r={'file':n,'url':u}
 try:
  a=requests.get(u,headers={'Accept-Encoding':'identity'},timeout=(15,90));p=O/'raw'/n;p.write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content),sha256=hashlib.sha256(a.content).hexdigest(),headers={k:a.headers.get(k) for k in ['Content-Length','Content-Encoding','Content-Type']})
  if n.endswith('.gz') and a.status_code==200:gzip.decompress(a.content);r['gzip_crc']=True
 except Exception as e:r['error']=str(e)
 print(json.dumps(r),flush=True);return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(get,urls.items()))
(O/'metadata/TRANSPORT_PROBE_29.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
