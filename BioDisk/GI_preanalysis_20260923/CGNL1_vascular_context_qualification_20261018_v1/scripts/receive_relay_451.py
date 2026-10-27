from pathlib import Path
import sys,json,os,hashlib,fcntl,datetime
b=Path(__file__).resolve().parents[1];m=json.loads((b/'metadata/PILOT_DOWNLOAD_MANIFEST_440.json').read_text());dest=Path(m['destination']);items={x['name']:x for x in m['items']}
def save(p,d):
 tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(d,indent=2),encoding='utf-8');tmp.replace(p)
def append(p,offset,data,digest,maximum):
 assert hashlib.sha256(data).hexdigest()==digest,'chunk checksum mismatch'
 with p.open('a+b') as f:
  fcntl.flock(f,fcntl.LOCK_EX);f.seek(0,2);size=f.tell()
  if size==offset+len(data):
   f.seek(offset);assert f.read(len(data))==data,'duplicate differs';return size
  assert size==offset,'offset mismatch'
  assert size+len(data)<=maximum,'oversize'
  f.write(data);f.flush();os.fsync(f.fileno());return f.tell()
def digest(p):
 h=hashlib.md5()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(8388608),b''):h.update(x)
 return h.hexdigest()
act=sys.argv[1]
if act=='test':
 p=dest/'RELAY_APPEND_FIXTURE_451.bin';assert not p.exists()
 a=b'CGNL1 relay integrity\x00';z=b'resume suffix'
 assert append(p,0,a,hashlib.sha256(a).hexdigest(),1000)==len(a)
 assert append(p,0,a,hashlib.sha256(a).hexdigest(),1000)==len(a)
 for offset,data,sha in [(0,z,hashlib.sha256(z).hexdigest()),(len(a),z,'bad')]:
  try:append(p,offset,data,sha,1000)
  except AssertionError:pass
  else:raise AssertionError('invalid write accepted')
 assert append(p,len(a),z,hashlib.sha256(z).hexdigest(),1000)==len(a+z)
 assert p.read_bytes()==a+z
 save(b/'metadata/RELAY_APPEND_QA_451.json',{'passed':True,'checks':['append','idempotent replay','wrong offset rejected','bad hash rejected','resume suffix'],'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()});print('PASS')
elif act in ['status','append']:
 name=sys.argv[2];item=items[name];target=dest/name;part=dest/(name+'.part')
 if act=='status':print(json.dumps({'complete':target.exists(),'offset':target.stat().st_size if target.exists() else part.stat().st_size if part.exists() else 0}))
 else:
  old=json.loads((b/'metadata/PILOT441.running.json').read_text());assert not Path('/proc/'+str(old['worker_pid'])).exists(),'old writer active'
  offset=int(sys.argv[3]);sha=sys.argv[4];data=sys.stdin.buffer.read(33554433);assert 0<len(data)<=33554432
  n=append(part,offset,data,sha,item['bytes']);save(b/'metadata'/('RELAY451_'+name+'.json'),{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':n,'total':item['bytes'],'chunk_sha256':sha});print(json.dumps({'offset':n}))
elif act=='finalize':
 receipts=[]
 for item in m['items']:
  target=dest/item['name'];part=dest/(item['name']+'.part');p=target if target.exists() else part
  assert p.stat().st_size==item['bytes'];assert digest(p)==item['md5'],'source MD5 mismatch'
  if p==part:part.rename(target)
  receipts.append(dict(item,verified=True,local_path=str(target)))
 save(b/'metadata/PILOT_DOWNLOAD_COMPLETE_451.json',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'run':m['run'],'files':receipts,'all_md5_verified':True});print('COMPLETE')
elif act=='failure':save(b/'metadata/LOCAL_RELAY_FAILED_451.json',{'error':sys.stdin.read(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
else:raise ValueError(act)
