from pathlib import Path
import requests,json,csv,re,hashlib,datetime,time,concurrent.futures,threading,shutil,os
O=Path('D:/CGNL1_DRA_public_data_20261016');O.mkdir(exist_ok=True);C=O/'ranges';C.mkdir(exist_ok=True)
meta=Path('C:/codex_test/output/CGNL1_prior_qualification_20261016/ENA_DRA_transport_metadata.tsv');rows=list(csv.DictReader(meta.open(encoding='utf-8'),delimiter='\t'));tasks=[]
for r in rows:
 if r['library_layout']!='PAIRED' or r['run_accession']=='DRR400693':continue
 condition='control' if 'vehicle control' in r['sample_title'] else ('combined' if 'TGF-B1 and IL-1B' in r['sample_title'] else ('TGF' if 'TGF-B1' in r['sample_title'] else 'IL1'))
 for u,n,h in zip(r['fastq_ftp'].split(';'),r['fastq_bytes'].split(';'),r['fastq_md5'].split(';')):tasks.append({'run':r['run_accession'],'donor':re.search(r'patient (GC\d+)',r['sample_title']).group(1),'condition':condition,'url':'https://'+u,'bytes':int(n),'md5':h,'name':u.rsplit('/',1)[1]})
assert len(tasks)==22;tasks.sort(key=lambda x:({'control':0,'combined':1,'TGF':2,'IL1':3}[x['condition']], {'GC42':0,'GC41':1,'GC40':2}[x['donor']], x['name']));assert shutil.disk_usage(O).free>sum(t['bytes'] for t in tasks)*2.2
def save(p,x):
 q=p.with_name(p.name+'.writing');q.write_text(json.dumps(x,indent=2),encoding='utf-8');q.replace(p)
save(O/'PLAN.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'tasks':tasks,'workers':8,'chunk_bytes':8388608,'role':'independent local public-source staging only; do not alter live server writers; import only whole files with officialMD5'})
lock=threading.Lock();results=[]
for task in tasks:
 p=O/task['name'];folder=C/task['name'];folder.mkdir(exist_ok=True);manifest={};start=time.monotonic();record=dict(task)
 def part(iv):
  a,b=iv;dst=folder/f'{a}_{b}.bin'
  for attempt in range(1,4):
   temp=dst.with_name(dst.name+f'.attempt{attempt}.partial')
   try:
    n=0;h=hashlib.sha256();begin=time.monotonic()
    with requests.get(task['url'],headers={'Range':f'bytes={a}-{b-1}','Accept-Encoding':'identity'},stream=True,timeout=(20,120)) as r:
     assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {a}-{b-1}/{task["bytes"]}'
     with temp.open('wb') as f:
      for data in r.iter_content(65536):
       f.write(data);n+=len(data);h.update(data)
       if n>b-a or time.monotonic()-begin>3600:raise TimeoutError('bounded chunk')
    assert n==b-a;temp.replace(dst)
    with lock:manifest[str(a)]={'bytes':n,'sha256':h.hexdigest(),'success':True};save(folder/'manifest.json',manifest)
    return dst
   except Exception as e:
    error=repr(e)
  with lock:manifest[str(a)]={'success':False,'error':error};save(folder/'manifest.json',manifest)
  raise RuntimeError(error)
 try:
  assert not p.exists()
  ivs=[(a,min(a+8388608,task['bytes'])) for a in range(0,task['bytes'],8388608)]
  with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:files=list(ex.map(part,ivs))
  assembled=p.with_name(p.name+'.assembled.partial');h=hashlib.md5();n=0
  with assembled.open('wb') as f:
   for q in files:
    with q.open('rb') as src:
     while data:=src.read(4194304):f.write(data);h.update(data);n+=len(data)
  assert n==task['bytes'] and h.hexdigest()==task['md5'];assembled.replace(p);record.update(success=True,actual_bytes=n,actual_md5=h.hexdigest())
 except Exception as e:record.update(success=False,error=repr(e))
 record.update(seconds=time.monotonic()-start,end_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());results.append(record);save(O/'RECEIPTS.json',results);print(json.dumps(record),flush=True)
save(O/'COMPLETE.json',{'all_pass':len(results)==22 and all(x['success'] for x in results),'files':len(results),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
