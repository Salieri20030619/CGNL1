import pathlib,json,requests,concurrent.futures,hashlib,gzip,os,time,datetime
p=pathlib.Path.cwd();old=p/'raw/external_239';out=p/'raw/external_250';out.mkdir(exist_ok=True)
receipt=json.loads((p/'metadata/EXTERNAL_DOWNLOAD_RECEIPTS_239.json').read_text())[0];matrix=pathlib.Path(receipt['path']);assert matrix.stat().st_size==receipt['bytes'];assert hashlib.sha256(matrix.read_bytes()).hexdigest()==receipt['sha256'];os.link(matrix,out/matrix.name)
url='https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL18nnn/GPL18990/soft/GPL18990_family.soft.gz';size=117719368;name='GPL18990_family.soft.gz';folder=out/(name+'.ranges');folder.mkdir(exist_ok=True)
def piece(a):
 z=min(a+1048576,size)-1;f=folder/(str(a)+'.bin');prior=old/(name+'.ranges')/(str(a)+'.bin')
 if prior.exists() and prior.stat().st_size==z-a+1:os.link(prior,f);return f
 for attempt in range(1,4):
  temp=folder/(str(a)+'.attempt'+str(attempt)+'.partial');s=requests.Session();s.trust_env=False
  try:
   with s.get(url,headers={'Range':f'bytes={a}-{z}','Accept-Encoding':'identity'},timeout=(20,120),stream=True) as r:
    assert r.status_code==206 and r.headers['Content-Range']==f'bytes {a}-{z}/{size}'
    with temp.open('xb') as w:
     for chunk in r.iter_content(131072):w.write(chunk)
   assert temp.stat().st_size==z-a+1;temp.rename(f);return f
  except Exception:
   if attempt==3:raise
   time.sleep(15*attempt)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:parts=list(ex.map(piece,range(0,size,1048576)))
f=out/name;temp=out/(name+'.assembled.partial');h=hashlib.sha256()
with temp.open('xb') as w:
 for part in parts:
  data=part.read_bytes();w.write(data);h.update(data)
assert temp.stat().st_size==size
with gzip.open(temp,'rb') as g:
 while g.read(1048576):pass
temp.rename(f);receipt['path']=str(out/matrix.name);receipts=[receipt,dict(url=url,path=str(f),bytes=size,sha256=h.hexdigest(),gzip_CRC_pass=True)]
(p/'metadata/EXTERNAL_DOWNLOAD_RECEIPTS_250.json').write_text(json.dumps(receipts,indent=2),encoding='utf-8')
x={'all_pass':True,'files':2,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_ranges_reused':sum(1 for part in parts if part.stat().st_nlink>1)}
(p/'metadata/EXTERNAL_DOWNLOAD_COMPLETE_250.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x),flush=True)
