from pathlib import Path
import csv,json,requests,hashlib,datetime,concurrent.futures,threading,os,time,re,shutil
O=Path.cwd()/'mechanism_investigation_v2';R=O/'raw/DRA014736';R.mkdir(exist_ok=True)
with (O/'metadata/DRA014736_runs.tsv').open(encoding='utf-8') as f:rows=list(csv.DictReader(f,delimiter='\t'))
rows=[r for r in rows if r['library_layout']=='PAIRED'];assert len(rows)==12
tasks=[]
for r in rows:
 assert r['study_accession']=='PRJDB14249';donor=re.search(r'patient (GC\d+)',r['sample_title']).group(1)
 condition='control' if 'vehicle control' in r['sample_title'] else ('combined' if 'TGF-B1 and IL-1B' in r['sample_title'] else ('TGF' if 'TGF-B1' in r['sample_title'] else 'IL1'))
 for u,n,h in zip(r['fastq_ftp'].split(';'),r['fastq_bytes'].split(';'),r['fastq_md5'].split(';')):tasks.append({'run':r['run_accession'],'donor':donor,'condition':condition,'url':'https://'+u,'bytes':int(n),'md5':h,'name':u.rsplit('/',1)[1]})
tasks.sort(key=lambda r:({'control':0,'combined':1,'TGF':2,'IL1':3}[r['condition']],r['donor'],r['name']))
assert shutil.disk_usage(R).free>sum(t['bytes'] for t in tasks)*2
(O/'metadata/DRA_DOWNLOAD_PLAN.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_donors':3,'n_libraries':12,'fastq_files':24,'total_bytes':sum(t['bytes'] for t in tasks),'exclude':'singleNanoporelibrary not needed for primary expression analysis','tasks':tasks},indent=2),encoding='utf-8')
lock=threading.Lock();out=[]
def get(t):
 p=R/t['name'];temp=p.with_suffix(p.suffix+'.partial');r=dict(t,start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());start=time.monotonic()
 try:
  assert not p.exists() and not temp.exists(),'existing file; inspect before restart'
  h=hashlib.md5();size=0
  with requests.get(t['url'],stream=True,timeout=(20,120)) as a:
   a.raise_for_status()
   with temp.open('wb') as f:
    for chunk in a.iter_content(1024*1024):
     if chunk:f.write(chunk);h.update(chunk);size+=len(chunk)
     if time.monotonic()-start>14400:raise TimeoutError('Four-hour file transfer limit')
  assert size==t['bytes'],f'size mismatch {size}'
  assert h.hexdigest()==t['md5'],'ENA MD5 mismatch'
  temp.rename(p);r.update(success=True,actual_bytes=size,actual_md5=h.hexdigest())
 except Exception as e:r.update(success=False,error=repr(e))
 r['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 with lock:
  out.append(r);(O/'metadata/DRA_DOWNLOAD_RECEIPTS.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(t['name'],r['success'],r.get('error','MD5 verified'),flush=True)
 return r
print('PLAN',len(tasks),'files',sum(t['bytes'] for t in tasks),'bytes',flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:results=list(ex.map(get,tasks))
(O/'metadata/DRA_DOWNLOAD_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':all(r['success'] for r in results),'files':len(results),'passed':sum(r['success'] for r in results)},indent=2),encoding='utf-8')
