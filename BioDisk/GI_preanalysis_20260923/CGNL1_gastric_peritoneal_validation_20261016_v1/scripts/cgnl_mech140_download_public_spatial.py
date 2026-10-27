from pathlib import Path
import requests,json,hashlib,zipfile,datetime,time,concurrent.futures,threading,shutil,os
O=Path('C:/codex_test/output/CGNL1_spatial_public_20261016');q=json.loads((O/'figshare_api.json').read_text(encoding='utf-8'))['data'];assert q['is_public'] and not q['download_disabled'];files=[x for name in ['Human_spatial.zip','Human_scRNA.zip'] for x in q['files'] if x['name']==name];assert len(files)==2 and shutil.disk_usage(O).free>sum(x['size'] for x in files)*2.5
def save(p,v):
 t=p.with_name(p.name+'.writing');t.write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf-8');t.replace(p)
save(O/'DOWNLOAD_PLAN_140.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'article':q['doi'],'files':files,'integrity':'repositorycomputedMD5+completeZIPCRC','workers':4});receipts=[];lock=threading.Lock()
for item in files:
 name=item['name'];folder=O/(name+'.ranges');folder.mkdir(exist_ok=True);manifest={};record=dict(item,start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 def part(iv):
  start,end=iv;dst=folder/f'{start}_{end}.bin'
  for attempt in range(1,4):
   tmp=dst.with_name(dst.name+f'.attempt{attempt}.partial')
   try:
    h=hashlib.sha256();n=0;begin=time.monotonic()
    with requests.get(item['download_url'],headers={'Range':f'bytes={start}-{end-1}','Accept-Encoding':'identity'},stream=True,timeout=(20,90)) as a:
     assert a.status_code==206 and a.headers.get('Content-Range')==f'bytes {start}-{end-1}/{item["size"]}',(a.status_code,a.headers.get('Content-Range'))
     with tmp.open('wb') as f:
      for b in a.iter_content(65536):
       f.write(b);n+=len(b);h.update(b)
       if n>end-start or time.monotonic()-begin>1800:raise TimeoutError('chunk bound')
    assert n==end-start;tmp.replace(dst)
    with lock:manifest[str(start)]={'bytes':n,'sha256':h.hexdigest(),'success':True};save(folder/'manifest.json',manifest)
    return dst
   except Exception as e:err=repr(e)
  with lock:manifest[str(start)]={'success':False,'error':err};save(folder/'manifest.json',manifest)
  raise RuntimeError(err)
 try:
  ivs=[(a,min(a+8388608,item['size'])) for a in range(0,item['size'],8388608)]
  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:parts=list(ex.map(part,ivs))
  p=O/name;tmp=p.with_name(name+'.assembled.partial');h=hashlib.md5();n=0
  with tmp.open('wb') as f:
   for partfile in parts:
    with partfile.open('rb') as src:
     while b:=src.read(4194304):f.write(b);h.update(b);n+=len(b)
  assert n==item['size'] and h.hexdigest()==item['computed_md5']
  with zipfile.ZipFile(tmp) as z:assert z.testzip() is None;members=[{'name':i.filename,'bytes':i.file_size,'compressed_bytes':i.compress_size} for i in z.infolist()]
  assert not p.exists();tmp.rename(p);record.update(success=True,actual_md5=h.hexdigest(),members=members)
 except Exception as e:record.update(success=False,error=repr(e))
 record['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();receipts.append(record);save(O/'DOWNLOAD_RECEIPTS_140.json',receipts);print(json.dumps(record,ensure_ascii=False),flush=True)
save(O/'DOWNLOAD_COMPLETE_140.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':len(receipts)==2 and all(r['success'] for r in receipts),'files':len(receipts)})
