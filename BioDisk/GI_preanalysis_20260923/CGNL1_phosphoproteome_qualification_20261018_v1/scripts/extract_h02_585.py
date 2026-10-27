from pathlib import Path
import os,signal,json,subprocess,struct,zlib,hashlib,datetime,traceback
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
M=B/'metadata/resume_581';D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576/H02/selected_585');D.mkdir(exist_ok=True)
def save(name,x):(M/name).write_text(json.dumps(x,indent=2),encoding='utf-8')
run=json.loads((M/'RUN.json').read_text());pid=run['pid']
proc=Path('/proc')/str(pid)
if proc.exists() and b'resume_h02_581.py' in (proc/'cmdline').read_bytes():
 os.kill(pid,signal.SIGSTOP)
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   if int((p/'stat').read_text().split(') ')[1].split()[1])==pid and b'curl' in (p/'cmdline').read_bytes():os.kill(int(p.name),signal.SIGTERM)
  except (FileNotFoundError,ProcessLookupError):pass
 os.kill(pid,signal.SIGTERM);os.kill(pid,signal.SIGCONT)
save('FULL_ARCHIVE_SUPERSEDED585.json',dict(reason='Verified remote ZIP directory allows retrieval of all four required original quantification/parameter members (81 MB compressed) without unrelated 17.8 GB archive',old_parts_preserved=True,whole_archive_complete=False))
rows=json.loads((M/'REMOTE_ZIP_DIRECTORY584.json').read_text());url=run['url'];records=[]
try:
 for j,item in enumerate(rows):
  if Path(item['name']).name.lower() not in ['phospho (sty)sites.txt','proteingroups.txt','summary.txt','parameters.txt']:continue
  a=item['offset'];b=rows[j+1]['offset']-1;p=D/(Path(item['name']).name+'.compressed_range');h=p.with_suffix('.headers')
  r=subprocess.run(['curl','--proxy','http://127.0.0.1:7897','--noproxy','','--fail','-sS','--retry','2','--connect-timeout','20','--max-time','600','--range',f'{a}-{b}','-D',str(h),'-o',str(p),url],capture_output=True,text=True)
  assert r.returncode==0,(r.returncode,r.stderr)
  assert p.stat().st_size==b-a+1
  assert f'content-range: bytes {a}-{b}/{run["expected_bytes"]}' in h.read_text().lower()
  data=p.read_bytes();header=struct.unpack('<IHHHHHIIIHH',data[:30]);assert header[0]==0x04034b50 and header[3]==8
  n,e=header[-2:];assert data[30:30+n].decode()==item['name'];compressed=data[30+n+e:30+n+e+item['compressed']]
  raw=zlib.decompress(compressed,-15);assert len(raw)==item['bytes'];assert zlib.crc32(raw)==item['crc']
  target=D/Path(item['name']).name;target.write_bytes(raw)
  records.append(dict(file=str(target),bytes=len(raw),crc32=item['crc'],crc_pass=True,sha256=hashlib.sha256(raw).hexdigest(),range=[a,b]))
  save('MEMBERS_PROGRESS585.json',records);print('verified',target.name,len(raw),flush=True)
 save('MEMBERS_COMPLETE585.json',dict(exit_code=0,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source=url,files=records,all_required_members_verified=len(records)==4,whole_archive_downloaded=False))
except Exception:
 save('MEMBERS_EXIT585.json',dict(exit_code=1,error=traceback.format_exc()));raise
