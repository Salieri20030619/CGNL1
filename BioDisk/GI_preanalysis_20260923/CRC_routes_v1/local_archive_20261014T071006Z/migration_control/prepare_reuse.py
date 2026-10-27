from pathlib import Path
R=Path(r'D:\CGNL1\_transfer')
s=(R/'migrate.py').read_text(encoding='utf-8')
s=s.replace("ARCH=REMOTE+'/local_archive_'+STAMP", "ARCH=json.loads((H/'MIGRATION_STATUS.json').read_text(encoding='utf-8'))['archive']")
s=s.replace("assert p.is_file() and p.stat().st_size==size\n  raw.append({'local':str(p),'archive_relative':'Downloads/SCP2079_counts/'+name,'bytes':size,'sha256':digest(p)})", "old={x['local']:x for x in json.loads((H/'RAW_TRANSFER_MANIFEST.json').read_text(encoding='utf-8'))}\n  raw.append(old[str(p)])")
s=s.replace("run(['scp','-o','BatchMode=yes','-o','ConnectTimeout=20']+[x['local'] for x in raw]+['sgs10:'+ARCH+'/Downloads/SCP2079_counts/'])", "assert json.loads((H/'TRANSFER_HANDOFF.json').read_text(encoding='utf-8'))['raw_scp_exit_code']==0")
s=s.replace("assert digest(p)==x['sha256'];p.unlink()", "if p.exists():\n   assert digest(p)==x['sha256'];p.unlink()")
old="status('streaming_local_archive',files=len(manifest),bytes=sum(x['bytes'] for x in manifest))\n cmd='tar --keep-old-files -xf - -C '+shlex.quote(ARCH)"
new="""run(['scp','-o','BatchMode=yes',str(H/'LOCAL_ARCHIVE_MANIFEST.json'),'sgs10:'+ARCH+'/'])
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
 reuse=json.loads(ssh('ARCH='+repr(ARCH)+'\\n'+code));save('SERVER_REUSE.json',reuse)
 status('streaming_local_archive',files=len(manifest),reused_files=reuse['reused_files'],bytes_to_upload=reuse['missing_bytes'])
 cmd='tar --keep-old-files -xzf - -C '+shlex.quote(ARCH)"""
assert old in s;s=s.replace(old,new)
s=s.replace("mode='w|'", "mode='w|gz'")
s=s.replace("for x in manifest:tar.add(x['local'],arcname=x['archive_relative'],recursive=False)", "for x in reuse['missing']:tar.add(x['local'],arcname=x['archive_relative'],recursive=False)")
(R/'migrate_reuse.py').write_text(s,encoding='utf-8')
compile(s,str(R/'migrate_reuse.py'),'exec')
print('Prepared hash-based server reuse migration')