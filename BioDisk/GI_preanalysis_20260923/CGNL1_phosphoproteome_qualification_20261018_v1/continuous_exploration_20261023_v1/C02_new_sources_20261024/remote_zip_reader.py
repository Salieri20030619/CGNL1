import io,urllib.request,json,zipfile,pathlib,re,concurrent.futures
P=pathlib.Path(r'C:\Users\31937\cgnl_transfer\restart_20261024')
class Remote(io.RawIOBase):
 def __init__(self,url):
  self.url=url;self.pos=0;self.cache={};self.transferred=0
  r=urllib.request.urlopen(urllib.request.Request(url,headers={'Range':'bytes=0-0'}),timeout=45)
  assert r.status==206,(r.status,dict(r.headers));self.size=int(r.headers['Content-Range'].split('/')[-1]);r.read()
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,off,whence=0):
  self.pos=off if whence==0 else self.pos+off if whence==1 else self.size+off;return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  n=min(n,self.size-self.pos)
  if n==0:return b''
  start=self.pos;end=start+n-1;key=(start,end)
  if key not in self.cache:
   r=urllib.request.urlopen(urllib.request.Request(self.url,headers={'Range':f'bytes={start}-{end}'}),timeout=90)
   assert r.status==206
   b=r.read();assert len(b)==n,(len(b),n);self.cache[key]=b;self.transferred+=len(b)
  self.pos+=n;return self.cache[key]
def inspect(label):
 try:
  url=json.loads((P/(label+'_url.json')).read_text())['data'];remote=Remote(url)
  with zipfile.ZipFile(remote) as z:
   entries=[dict(name=x.filename,size=x.file_size,compressed=x.compress_size,crc=x.CRC) for x in z.infolist()]
  out=dict(label=label,archive_bytes=remote.size,index_bytes_transferred=remote.transferred,entries=entries)
  (P/(label+'_zip_index.json')).write_text(json.dumps(out,indent=2));return {**out,'entries':[x for x in entries if any(t in x['name'].lower() for t in ['pdac','readme'])]}
 except Exception as e:return dict(label=label,error=repr(e))
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
  for x in ex.map(inspect,['BCM','UMich']):print(json.dumps(x))
