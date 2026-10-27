from pathlib import Path
import subprocess,json,hashlib,datetime,tarfile,shlex,base64,mimetypes,shutil,os,traceback
from bs4 import BeautifulSoup
BASE=Path(r'D:\CGNL1');W=BASE/'_transfer';H=BASE/'handoff';U=Path(r'C:\Users\31937')
REMOTE='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1'
STAMP=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
ARCH=json.loads((H/'MIGRATION_STATUS.json').read_text(encoding='utf-8'))['archive']
ANALYSIS=REMOTE+'/SCP2079_paired_20261014_v1'
FLAGS=getattr(subprocess,'CREATE_NO_WINDOW',0)
def run(args,**kw):
 p=subprocess.run(args,capture_output=True,creationflags=FLAGS,**kw)
 if p.returncode:raise RuntimeError(str(args[:3])+': '+p.stderr.decode('utf-8'))
 return p.stdout.decode('utf-8')
def ssh(code):return run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10','python3 -'],input=code.encode('utf-8'))
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(name,data):
 (H/name).write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
def status(stage,**kwargs):
 d={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':stage,'archive':ARCH,'analysis':ANALYSIS,**kwargs};save('MIGRATION_STATUS.json',d);print(json.dumps(d,ensure_ascii=True),flush=True)
try:
 roots=sorted(p for p in U.glob('CGNL1*') if p.is_dir())+[U/'CRC_routes_20260923']
 for p in roots:
  assert p.resolve().parent==U.resolve() and (p.name.startswith('CGNL1') or p.name=='CRC_routes_20260923')
  assert not p.is_symlink()
 files=[]
 for root in roots:
  for p in sorted(root.rglob('*')):
   assert not p.is_symlink(),str(p)
   if p.is_file():files.append((p,p.relative_to(U).as_posix()))
 raw=[]
 expected={'matrix.mtx.gz':1111632586,'barcodes.tsv.gz':429479,'features.tsv.gz':333437}
 for name,size in expected.items():
  p=U/'Downloads'/name
  old={x['local']:x for x in json.loads((H/'RAW_TRANSFER_MANIFEST.json').read_text(encoding='utf-8'))}
  raw.append(old[str(p)])
 save('RAW_TRANSFER_MANIFEST.json',raw)
 ssh('from pathlib import Path\nfor p in '+repr([ARCH+'/Downloads/SCP2079_counts',ANALYSIS+'/raw/expression_epi_counts',ANALYSIS+'/scripts',ANALYSIS+'/logs',ANALYSIS+'/results'])+': Path(p).mkdir(parents=True,exist_ok=True)\n')
 status('uploading_downloaded_matrix',raw_bytes=sum(x['bytes'] for x in raw))
 assert json.loads((H/'TRANSFER_HANDOFF.json').read_text(encoding='utf-8'))['raw_scp_exit_code']==0
 code='''
from pathlib import Path
import hashlib,json,os,datetime
for x in ITEMS:
 p=Path(ARCH)/x['archive_relative']
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 assert p.stat().st_size==x['bytes'] and h.hexdigest()==x['sha256']
 dest=Path(ANALYSIS)/'raw/expression_epi_counts'/p.name
 if dest.exists():assert dest.stat().st_ino==p.stat().st_ino or hashlib.sha256(dest.read_bytes()).hexdigest()==x['sha256']
 else:os.link(p,dest)
print(json.dumps({'all_pass':True,'analysis':ANALYSIS,'raw_files':ITEMS}))
'''
 out=json.loads(ssh('ITEMS='+repr(raw)+'\nARCH='+repr(ARCH)+'\nANALYSIS='+repr(ANALYSIS)+'\n'+code))
 save('RAW_SERVER_READY.json',out);status('raw_verified_preparing_archive')
 # Remove the downloaded raw files only after remote hash verification and local re-check.
 for x in raw:
  p=Path(x['local']).resolve();assert p.parent==(U/'Downloads').resolve() and p.name in expected
  if p.exists():
   assert digest(p)==x['sha256'];p.unlink()
 save('RAW_LOCAL_CLEANUP.json',{'all_verified_before_delete':True,'removed':[x['local'] for x in raw]})
 manifest=[]
 for p,rel in files:manifest.append({'local':str(p),'archive_relative':rel,'bytes':p.stat().st_size,'sha256':digest(p)})
 save('LOCAL_ARCHIVE_MANIFEST.json',manifest)
 run(['scp','-o','BatchMode=yes',str(H/'LOCAL_ARCHIVE_MANIFEST.json'),'sgs10:'+ARCH+'/'])
 code='''
from pathlib import Path
import os,json,hashlib,shutil,collections
R=Path(ARCH);B=R.parent
items=json.loads((R/'LOCAL_ARCHIVE_MANIFEST.json').read_text(encoding='utf-8'))
byname=collections.defaultdict(list);bysize=collections.defaultdict(list)
want={(Path(x['archive_relative']).name,x['bytes']) for x in items};sizes={x['bytes'] for x in items if x['bytes']>1000000}
for root,dirs,files in os.walk(B):
 dirs[:]=[d for d in dirs if not d.startswith('local_archive_') and d!='SCP2079_paired_20261014_v1']
 for name in files:
  p=Path(root)/name
  if p.is_symlink():continue
  size=p.stat().st_size
  if (name,size) in want:byname[name,size].append(p)
  if size in sizes:bysize[size].append(p)
cache={}
def sha(p):
 if p not in cache:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
  cache[p]=h.hexdigest()
 return cache[p]
missing=[];reused=[]
for x in items:
 dst=R/x['archive_relative'];candidates=byname[Path(x['archive_relative']).name,x['bytes']]
 if not candidates and x['bytes']>1000000:candidates=bysize[x['bytes']]
 if dst.exists():
  assert dst.stat().st_size==x['bytes'] and sha(dst)==x['sha256']
  reused.append({'destination':x['archive_relative'],'source':'already_archived'});continue
 found=False
 for p in candidates:
  if sha(p)==x['sha256']:
   dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
   reused.append({'destination':x['archive_relative'],'source':str(p)});found=True;break
 if not found:missing.append(x)
out={'missing':missing,'reused':reused,'reused_files':len(reused),'missing_bytes':sum(x['bytes'] for x in missing)}
(R/'SERVER_REUSE.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out))
'''
 reuse=json.loads(ssh('ARCH='+repr(ARCH)+'\n'+code));save('SERVER_REUSE.json',reuse)
 status('streaming_local_archive',files=len(manifest),reused_files=reuse['reused_files'],bytes_to_upload=reuse['missing_bytes'])
 cmd='tar --keep-old-files -xzf - -C '+shlex.quote(ARCH)
 with (W/'archive_transfer.stderr.log').open('wb') as err:
  proc=subprocess.Popen(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','sgs10',cmd],stdin=subprocess.PIPE,stdout=subprocess.DEVNULL,stderr=err,creationflags=FLAGS)
  with tarfile.open(fileobj=proc.stdin,mode='w|gz',format=tarfile.PAX_FORMAT) as tar:
   for x in reuse['missing']:tar.add(x['local'],arcname=x['archive_relative'],recursive=False)
  proc.stdin.close();rc=proc.wait()
  assert rc==0,(rc,(W/'archive_transfer.stderr.log').read_text(encoding='utf-8'))
 run(['scp','-o','BatchMode=yes',str(H/'LOCAL_ARCHIVE_MANIFEST.json'),str(H/'RAW_TRANSFER_MANIFEST.json'),'sgs10:'+ARCH+'/'])
 status('verifying_server_archive')
 code='''
from pathlib import Path
import hashlib,json,datetime
R=Path(ARCH)
items=json.loads((R/'LOCAL_ARCHIVE_MANIFEST.json').read_text(encoding='utf-8'))
for x in items:
 p=R/x['archive_relative'];h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 assert p.stat().st_size==x['bytes'] and h.hexdigest()==x['sha256'],str(p)
d={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'all_pass':True,'files':len(items),'bytes':sum(x['bytes'] for x in items),'archive':ARCH}
(R/'ARCHIVE_VERIFICATION.json').write_text(json.dumps(d,indent=2),encoding='utf-8');print(json.dumps(d))
'''
 verified=json.loads(ssh('ARCH='+repr(ARCH)+'\n'+code));save('ARCHIVE_VERIFICATION.json',verified)
 # Retain final report reading copies only; embed images and mark server-only attachments.
 retained=[]
 for x in manifest:
  p=Path(x['local']);rel=Path(x['archive_relative'])
  is_report=p.suffix.lower() in ['.html','.pdf','.docx'] and ('report' in p.name.lower() or '报告' in p.name or 'qualification_report' in p.name.lower())
  if not is_report or any(k in p.parts for k in ['sources','metadata','raw','downloads','internal_repair_archive']):continue
  dest=BASE/'final'/rel
  dest.parent.mkdir(parents=True,exist_ok=True)
  if p.suffix.lower()=='.html':
   content=p.read_text(encoding='utf-8');soup=BeautifulSoup(content,'html.parser')
   for tag in soup.find_all(['img','script','link']):
    attr='href' if tag.name=='link' else 'src';s=tag.get(attr,'')
    if not s or s.startswith(('data:','http:','https:','#')):continue
    resource=(p.parent/s).resolve()
    if resource.is_file():
     if tag.name=='img':tag[attr]='data:'+(mimetypes.guess_type(resource.name)[0] or 'application/octet-stream')+';base64,'+base64.b64encode(resource.read_bytes()).decode('ascii')
     elif tag.name=='script':tag.attrs.pop('src',None);tag.string=resource.read_text(encoding='utf-8')
     elif tag.name=='link':tag.name='style';tag.attrs={};tag.string=resource.read_text(encoding='utf-8')
   for a in soup.find_all('a',href=True):
    link=a['href']
    if not link.startswith(('https:','http:','mailto:','data:','#')):
     a['title']='Server archive: '+ARCH+'/'+str(rel.parent/Path(link)).replace('\\','/')
     a['data-server-attachment']=a['title'];del a['href']
   banner=soup.new_tag('p');banner.string='Local reading copy: original report and supporting analysis files are preserved on sgs10. Embedded figures are retained; attachment locations appear in link tooltips.'
   if soup.body:soup.body.insert(0,banner)
   else:soup.insert(0,banner)
   dest.write_text(str(soup),encoding='utf-8')
  else:shutil.copy2(p,dest)
  retained.append({'source':str(p),'local_final':str(dest),'server_original':ARCH+'/'+x['archive_relative'],'sha256':digest(dest)})
 save('FINAL_REPORTS_RETAINED.json',retained)
 status('archive_verified_cleaning_local',retained_reports=len(retained))
 # Revalidate exact bounded paths and file hashes, then delete only verified files.
 for x in manifest:
  p=Path(x['local']).resolve()
  assert any(root.resolve() in p.parents for root in roots)
  assert p.is_file() and digest(p)==x['sha256'],str(p)
 for x in manifest:Path(x['local']).unlink()
 for root in roots:
  resolved=root.resolve();assert resolved.parent==U.resolve() and resolved.name in [r.name for r in roots]
  for d in sorted([p for p in root.rglob('*') if p.is_dir()],key=lambda p:len(p.parts),reverse=True):d.rmdir()
  root.rmdir()
 cleanup={'all_source_files_verified_before_delete':True,'directories_removed':[str(x) for x in roots],'files_removed':len(manifest),'bytes_removed':sum(x['bytes'] for x in manifest),'downloads_removed':True,'unexpected_new_files_would_block_cleanup':True}
 save('LOCAL_CLEANUP_RECEIPT.json',cleanup)
 status('complete',files_archived=len(manifest)+len(raw),bytes_archived=sum(x['bytes'] for x in manifest+raw),retained_reports=len(retained))
except Exception:
 status('failed',traceback=traceback.format_exc());raise