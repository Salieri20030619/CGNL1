from pathlib import Path
import requests,json,hashlib,gzip,concurrent.futures,threading,time,datetime
p=Path.cwd();url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140973/matrix/GSE140973_series_matrix.txt.gz';size=24671339;folder=p/'raw/matrix297_ranges';folder.mkdir(exist_ok=True);tls=threading.local()
def get_session():
 if not hasattr(tls,'s'):tls.s=requests.Session();tls.s.trust_env=False
 return tls.s
def piece(a):
 z=min(a+1048576,size)-1;f=folder/(str(a)+'.bin')
 for attempt in range(1,4):
  tmp=folder/(str(a)+'.attempt'+str(attempt)+'.partial')
  try:
   with get_session().get(url,headers={'Range':f'bytes={a}-{z}','Accept-Encoding':'identity'},stream=True,timeout=(15,90)) as r:
    r.raise_for_status();assert r.status_code==206 and r.headers['Content-Range']==f'bytes {a}-{z}/{size}'
    with tmp.open('xb') as out:
     for data in r.iter_content(131072):out.write(data)
   assert tmp.stat().st_size==z-a+1;tmp.rename(f);return f
  except Exception:
   if attempt==3:raise
   time.sleep(15*attempt)
with concurrent.futures.ThreadPoolExecutor(max_workers=6) as ex:parts=list(ex.map(piece,range(0,size,1048576)))
f=p/'raw/GSE140973_series_matrix.txt.gz';tmp=f.with_suffix('.assembled.partial');h=hashlib.sha256()
with tmp.open('xb') as out:
 for part in parts:
  data=part.read_bytes();out.write(data);h.update(data)
assert tmp.stat().st_size==size
with gzip.open(tmp,'rb') as inp:
 while inp.read(1048576):pass
tmp.rename(f);receipts=[{'url':url,'path':str(f),'bytes':size,'sha256':h.hexdigest(),'gzip_CRC_pass':True}];(p/'metadata/DOWNLOAD_RECEIPTS_297.json').write_text(json.dumps(receipts,indent=2),encoding='utf-8');print(json.dumps(receipts[-1]),flush=True)
url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10332&targ=self&form=text&view=full';f=p/'raw/GPL10332_self_297.soft';tmp=f.with_suffix('.partial');h=hashlib.sha256();s=requests.Session();s.trust_env=False
with s.get(url,stream=True,timeout=(15,120)) as r:
 r.raise_for_status()
 with tmp.open('xb') as out:
  for data in r.iter_content(131072):out.write(data);h.update(data)
text=tmp.read_text(encoding='utf-8');assert '^PLATFORM = GPL10332' in text and '!platform_table_begin' in text and '!platform_table_end' in text and '^SAMPLE =' not in text;tmp.rename(f);receipts.append({'url':url,'path':str(f),'bytes':f.stat().st_size,'sha256':h.hexdigest(),'platform_self_complete':True});(p/'metadata/DOWNLOAD_RECEIPTS_297.json').write_text(json.dumps(receipts,indent=2),encoding='utf-8');print(json.dumps(receipts[-1]),flush=True)
(p/'metadata/DOWNLOAD_COMPLETE_297.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':True,'files':2}),encoding='utf-8')
