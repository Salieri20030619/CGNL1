from pathlib import Path
import os,signal,json,requests,hashlib,time,datetime,concurrent.futures,threading,gzip
O=Path.cwd()/'mechanism_investigation_v2';R=O/'reference/gencode_v47';D=O/'raw/DRA014736';C=O/'raw/verified_ranges';C.mkdir(exist_ok=True);lock=threading.Lock()
def save(p,x):
 tmp=p.with_name(p.name+'.writing');tmp.write_text(json.dumps(x,indent=2),encoding='utf-8');tmp.replace(p)
def hashfile(p,kind='md5'):
 h=hashlib.new(kind)
 with p.open('rb') as f:
  while b:=f.read(4*1024*1024):h.update(b)
 return h.hexdigest()
assert not (O/'metadata/INDEX_COMPLETE.json').exists() and not (R/'gentrome.fa').exists(),'index already started; inspect before intervention'
md={l.split()[1].lstrip('*').removeprefix('./'):l.split()[0] for l in (R/'MD5SUMS').read_text().splitlines() if len(l.split())==2}
jobs=[]
for n in ['gencode.v47.transcripts.fa.gz','GRCh38.primary_assembly.genome.fa.gz']:
 u='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/'+n
 with requests.get(u,headers={'Range':'bytes=0-0','Accept-Encoding':'identity'},stream=True,timeout=(15,30)) as a:
  assert a.status_code==206;total=int(a.headers['Content-Range'].split('/')[1]);assert len(a.raw.read(1))==1
 jobs.append({'name':n,'url':u,'bytes':total,'md5':md[n],'kind':'reference','directory':str(R)})
for t in json.loads((O/'metadata/DRA_DOWNLOAD_PLAN.json').read_text())['tasks']:jobs.append(dict(t,kind='DRA',directory=str(D)))
save(O/'metadata/RESUMABLE_TRANSFER_PLAN.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pid':os.getpid(),'jobs':jobs,'chunk_bytes':8*1024*1024,'workers':8,'reason':'measured~20KB/sperconnection,old4hcapincompatiblewithfilecompletion;verifiedRangeandprefixidentity','success_criterion':'officialMD5ofwholeassembledfile','no_analysis_of_partial':True})
stopped=[]
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  args=(p/'cmdline').read_bytes().split(b'\0')
  for name in ['cgnl_mech38_download_dra.py','cgnl_mech41_acquire_reference.py','cgnl_mech45b_index_when_ready.py','cgnl_mech46_quantify_when_ready.py']:
   if any(a==str(O.parent/'scripts'/name).encode() for a in args):os.kill(int(p.name),signal.SIGTERM);stopped.append({'pid':int(p.name),'script':name})
 except (FileNotFoundError,PermissionError):pass
save(O/'metadata/TRANSFER_UPGRADE_INTERVENTION.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stopped':stopped,'data_preserved':True})
# Wait locally for these specifically identified processes to release their files.
for rec in stopped:
 for _ in range(100):
  p=Path('/proc')/str(rec['pid'])
  if not p.exists():break
  time.sleep(.1)
 else:raise RuntimeError('Prior transfer did not exit; no shared writes allowed')
allresults=[];dr=[];rr=[]
def transfer(t):
 p=Path(t['directory'])/t['name'];prefix=p.with_name(p.name+'.preserved_prefix');part=p.with_name(p.name+'.partial');folder=C/t['name'];folder.mkdir(exist_ok=True);manifest={};result=dict(t,start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 try:
  if p.exists():assert p.stat().st_size==t['bytes'] and hashfile(p)==t['md5'];result.update(success=True,reused_existing=True);return result
  if part.exists():assert not prefix.exists();part.rename(prefix)
  offset=prefix.stat().st_size if prefix.exists() else 0;assert offset<t['bytes'];result['preserved_prefix_bytes']=offset
  intervals=[(s,min(s+8*1024*1024,t['bytes'])) for s in range(offset,t['bytes'],8*1024*1024)]
  def piece(iv):
   start,end=iv;dst=folder/(str(start)+'_'+str(end)+'.bin');item={'start':start,'end_exclusive':end,'bytes':end-start}
   for attempt in range(2):
    temp=dst.with_name(dst.name+f'.attempt{attempt+1}.partial')
    try:
     with requests.get(t['url'],headers={'Range':f'bytes={start}-{end-1}','Accept-Encoding':'identity'},stream=True,timeout=(20,120)) as a:
      assert a.status_code==206 and a.headers.get('Content-Range')==f'bytes {start}-{end-1}/{t["bytes"]}'
      size=0;h=hashlib.sha256();begin=time.monotonic()
      with temp.open('wb') as f:
       for b in a.iter_content(65536):
        if b:f.write(b);h.update(b);size+=len(b)
        if time.monotonic()-begin>3600:raise TimeoutError('Onehourchunkcap')
     assert size==end-start;temp.rename(dst);item.update(success=True,sha256=h.hexdigest(),attempt=attempt+1)
     with lock:manifest[str(start)]=item;save(folder/'range_manifest.json',manifest)
     return dst
    except Exception as e:item['error']=repr(e)
   with lock:manifest[str(start)]=dict(item,success=False);save(folder/'range_manifest.json',manifest)
   raise RuntimeError(item)
  with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:pieces=list(ex.map(piece,intervals))
  assembly=p.with_name(p.name+'.assembled.partial');h=hashlib.md5();size=0
  with assembly.open('wb') as out:
   for source in ([prefix] if offset else [])+pieces:
    with source.open('rb') as f:
     while b:=f.read(4*1024*1024):out.write(b);h.update(b);size+=len(b)
  assert size==t['bytes'] and h.hexdigest()==t['md5'],'WholefileofficialMD5mismatch';assembly.rename(p);result.update(success=True,actual_bytes=size,actual_md5=h.hexdigest(),chunks=len(pieces))
 except Exception as e:result.update(success=False,error=repr(e))
 result['end_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();return result
for i,t in enumerate(jobs):
 r=transfer(t);allresults.append(r);save(O/'metadata/RESUMABLE_TRANSFER_RECEIPTS.json',allresults);print(t['name'],r.get('success'),r.get('error','officialMD5verified'),flush=True)
 if t['kind']=='reference':
  rr.append(r)
  if len(rr)==2:save(O/'metadata/REFERENCE_DOWNLOAD_COMPLETE.json',rr)
 else:dr.append(r);save(O/'metadata/DRA_DOWNLOAD_RECEIPTS.json',dr)
save(O/'metadata/DRA_DOWNLOAD_COMPLETE.json',{'all_pass':len(dr)==24 and all(r['success'] for r in dr),'files':len(dr),'passed':sum(r['success'] for r in dr),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()});save(O/'metadata/RESUMABLE_TRANSFER_COMPLETE.json',{'all_pass':all(r['success'] for r in allresults),'files':len(allresults),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
