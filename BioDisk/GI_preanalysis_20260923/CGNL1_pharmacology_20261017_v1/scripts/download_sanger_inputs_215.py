
from pathlib import Path
import requests,json,hashlib,zipfile,csv,datetime,concurrent.futures,threading,time
p=Path.cwd();tasks=json.loads((p/'metadata/DOWNLOAD_PLAN_215.json').read_text())['tasks'];tls=threading.local();receipts=[]
def save(f,x):
 tmp=f.with_name(f.name+'.writing');tmp.write_text(json.dumps(x,indent=2),encoding='utf-8');tmp.replace(f)
def sess():
 if not hasattr(tls,'s'):tls.s=requests.Session();tls.s.trust_env=False
 return tls.s
for t in tasks:
 folder=p/'raw'/(t['name']+'.ranges');folder.mkdir(exist_ok=True)
 def piece(iv):
  a,b=iv;f=folder/(str(a)+'_'+str(b)+'.bin')
  for attempt in range(1,4):
   tmp=f.with_name(f.name+'.attempt'+str(attempt)+'.partial')
   try:
    with sess().get(t['url'],headers={'Range':'bytes='+str(a)+'-'+str(b-1),'Accept-Encoding':'identity'},stream=True,timeout=(20,90)) as r:
     assert r.status_code==206 and r.headers.get('Content-Range')=='bytes '+str(a)+'-'+str(b-1)+'/'+str(t['bytes'])
     n=0
     with tmp.open('xb') as out:
      for block in r.iter_content(131072):out.write(block);n+=len(block)
    assert n==b-a
    tmp.rename(f);return f
   except Exception:
    if attempt==3:raise
    time.sleep(30*attempt)
 try:
  intervals=[(a,min(a+2097152,t['bytes'])) for a in range(0,t['bytes'],2097152)]
  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:parts=list(ex.map(piece,intervals))
  dest=p/'raw'/t['name'];tmp=dest.with_name(dest.name+'.assembled.partial');h=hashlib.sha256()
  with tmp.open('xb') as out:
   for f in parts:
    with f.open('rb') as inp:
     while block:=inp.read(1048576):out.write(block);h.update(block)
  assert tmp.stat().st_size==t['bytes']
  if t['name'].endswith('.zip'):
   with zipfile.ZipFile(tmp) as z:assert z.testzip() is None;members=z.namelist()
   info={'zip_members':members,'zip_CRC_pass':True}
  else:
   with tmp.open(encoding='utf-8',newline='') as inp:
    reader=csv.reader(inp);header=next(reader);count=0
    for row in reader:assert len(row)==len(header),(t['name'],count,len(row),len(header));count+=1
   info={'rows':count,'columns':len(header),'header':header}
  tmp.rename(dest);r=dict(t,success=True,sha256=h.hexdigest(),qualification=info)
 except Exception as e:r=dict(t,success=False,error=repr(e))
 receipts.append(r);save(p/'metadata/DOWNLOAD_RECEIPTS_215.json',receipts);print(json.dumps(r),flush=True)
done={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':all(r['success'] for r in receipts) and len(receipts)==4,'files':len(receipts)}
save(p/'metadata/DOWNLOAD_COMPLETE_215.json',done);assert done['all_pass'],done
