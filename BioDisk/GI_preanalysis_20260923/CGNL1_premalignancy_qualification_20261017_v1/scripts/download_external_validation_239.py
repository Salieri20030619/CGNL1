import pathlib,requests,json,concurrent.futures,gzip,hashlib,datetime
p=pathlib.Path.cwd();out=p/'raw'/'external_239';out.mkdir(exist_ok=True)
urls=['https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78523/matrix/GSE78523_series_matrix.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL18nnn/GPL18990/soft/GPL18990_family.soft.gz']
receipts=[]
for url in urls:
 s=requests.Session();s.trust_env=False
 with s.get(url,headers={'Range':'bytes=0-0','Accept-Encoding':'identity'},timeout=(20,60)) as r:
  r.raise_for_status();assert r.status_code==206;size=int(r.headers['Content-Range'].split('/')[-1]);assert len(r.content)==1
 name=url.rsplit('/',1)[1];folder=out/(name+'.ranges');folder.mkdir(exist_ok=True)
 def piece(a):
  z=min(a+1048576,size)-1;f=folder/(str(a)+'.bin')
  ss=requests.Session();ss.trust_env=False
  with ss.get(url,headers={'Range':f'bytes={a}-{z}','Accept-Encoding':'identity'},timeout=(20,120),stream=True) as rr:
   assert rr.status_code==206 and rr.headers['Content-Range']==f'bytes {a}-{z}/{size}'
   with f.open('xb') as w:
    for chunk in rr.iter_content(131072):w.write(chunk)
  assert f.stat().st_size==z-a+1
  return f
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:parts=list(ex.map(piece,range(0,size,1048576)))
 f=out/name;h=hashlib.sha256()
 with f.open('xb') as w:
  for part in parts:
   data=part.read_bytes();w.write(data);h.update(data)
 assert f.stat().st_size==size
 with gzip.open(f,'rb') as g:
  while g.read(1048576):pass
 receipts.append(dict(url=url,path=str(f),bytes=size,sha256=h.hexdigest(),gzip_CRC_pass=True))
 print(json.dumps(receipts[-1]),flush=True)
 (p/'metadata/EXTERNAL_DOWNLOAD_RECEIPTS_239.json').write_text(json.dumps(receipts,indent=2),encoding='utf-8')
(p/'metadata/EXTERNAL_DOWNLOAD_COMPLETE_239.json').write_text(json.dumps({'all_pass':True,'files':len(receipts)}),encoding='utf-8')
