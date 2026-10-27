from pathlib import Path
import json,sys,hashlib,requests,re,sqlite3
b=Path(__file__).resolve().parents[1];info=json.loads((b/'metadata/APSW_SOURCE_367.json').read_text());sys.path.insert(0,info['installed_to']);import apsw
class ReadFile(apsw.VFSFile):
 def __init__(self,base,reader,size,flags):
  super().__init__('',str(base),[apsw.SQLITE_OPEN_READONLY,0]);self.reader=reader;self.size=size
 def xRead(self,amount,offset):return self.reader(offset,amount)
 def xFileSize(self):return self.size
 def xWrite(self,*a):raise apsw.ReadOnlyError('immutable source')
 def xTruncate(self,*a):raise apsw.ReadOnlyError('immutable source')
 def xDeviceCharacteristics(self):return apsw.SQLITE_IOCAP_IMMUTABLE
class ReadVFS(apsw.VFS):
 def __init__(self,name,base,reader,size):
  super().__init__(name,'');self.base=base;self.reader=reader;self.size=size
 def xOpen(self,name,flags):
  assert flags[0]&apsw.SQLITE_OPEN_MAIN_DB;flags[1]=apsw.SQLITE_OPEN_READONLY;return ReadFile(self.base,self.reader,self.size,flags)
class HTTPReader:
 def __init__(self,url,size,cache,cap=64*1024*1024):
  self.url=url;self.size=size;self.cache=cache;cache.mkdir(exist_ok=True,parents=True);self.cap=cap;self.total=0;self.nreq=0;self.block=65536;self.session=requests.Session();self.etag=None;self.receipts=[]
 def __call__(self,offset,amount):
  assert offset>=0 and offset+amount<=self.size
  chunks=[]
  for start in range(offset//self.block*self.block,offset+amount,self.block):
   end=min(start+self.block,self.size)-1;p=self.cache/(str(start)+'.bin')
   if p.exists():v=p.read_bytes();assert len(v)==end-start+1
   else:
    assert self.total+end-start+1<=self.cap,'pilot byte cap reached'
    r=self.session.get(self.url,headers={'Range':f'bytes={start}-{end}'},timeout=(10,40));r.raise_for_status();assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {start}-{end}/{self.size}';v=r.content;assert len(v)==end-start+1
    etag=r.headers.get('ETag');assert self.etag is None or self.etag==etag;self.etag=etag;p.write_bytes(v);self.nreq+=1;self.total+=len(v);self.receipts.append({'start':start,'end':end,'sha256':hashlib.sha256(v).hexdigest(),'etag':etag})
   chunks.append(v)
  return b''.join(chunks)[offset%self.block:offset%self.block+amount]
def connect(name,base,reader,size):
 vfs=ReadVFS(name,base,reader,size);con=apsw.Connection('file:remote?immutable=1',flags=apsw.SQLITE_OPEN_READONLY|apsw.SQLITE_OPEN_URI,vfs=name);con.execute('PRAGMA query_only=ON');con.execute('PRAGMA trusted_schema=OFF');con.execute('PRAGMA mmap_size=0');return vfs,con
if __name__=='__main__':
 # independent SQLite fixture tests cross-page data and unicode using same VFS interface
 fixture=b/'raw/vfs_fixture_368.sqlite';db=sqlite3.connect(fixture);db.execute('CREATE TABLE IF NOT EXISTS probe(k INTEGER PRIMARY KEY,t TEXT,v REAL)')
 if db.execute('select count(*) from probe').fetchone()[0]==0:db.executemany('insert into probe values(?,?,?)',[(i,'alpha'+chr(945)+(str(i)*200),i/7) for i in range(500)])
 db.commit();expected=db.execute('select k,t,v from probe order by k desc').fetchall();db.close();blob=fixture.read_bytes();v,c=connect('fixture368',fixture,lambda o,n:blob[o:o+n],len(blob));assert list(c.execute('select k,t,v from probe order by k desc'))==expected;c.close();v.unregister()
 pilot=json.loads((b/'metadata/MSF_PILOT_RANGE_366.json').read_text());size=int(pilot['content_range'].split('/')[-1]);reader=HTTPReader(pilot['url'],size,b/'raw/MSF_pilot_range_cache_368');v,c=connect('remote368',b/'raw/MSF_pilot_header_366.bin',reader,size);schema=list(c.execute('SELECT type,name,tbl_name,rootpage,sql FROM sqlite_master ORDER BY name'));c.close();v.unregister()
 out={'fixture_test_pass':True,'apsw':apsw.apswversion(),'source_url':pilot['url'],'source_size':size,'fetched_bytes':reader.total,'requests':reader.nreq,'schema':schema,'range_receipts':reader.receipts,'whole_file_hash_verified':False};(b/'metadata/MSF_PILOT_SCHEMA_368.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print('schema tables',[r[1] for r in schema if r[0]=='table']);print('fetched',reader.total,'requests',reader.nreq)
