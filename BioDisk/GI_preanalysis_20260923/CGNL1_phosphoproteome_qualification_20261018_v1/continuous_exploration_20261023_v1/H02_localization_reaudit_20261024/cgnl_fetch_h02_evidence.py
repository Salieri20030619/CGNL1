from pathlib import Path
import urllib.request,json,struct,hashlib,zlib,time,datetime
P=Path(r'C:\Users\31937\cgnl_transfer\restart_20261024')
index=json.loads((P/'H02_phospho_archive_index.json').read_text());entry=next(x for x in index['entries'] if x['name']=='txt/evidence.txt');url=index['url']
def fetch(a,b):
 r=urllib.request.urlopen(urllib.request.Request(url,headers={'Range':f'bytes={a}-{b}'}),timeout=90)
 assert r.status==206 and r.headers['Content-Range'].startswith(f'bytes {a}-{b}/');return r
with fetch(entry['offset'],entry['offset']+29) as r:h=r.read()
assert h[:4]==b'PK\x03\x04';name,extra=struct.unpack_from('<HH',h,26);start=entry['offset']+30+name+extra
out=P/'H02_evidence.deflate.partial';final=P/'H02_evidence.deflate'
assert not final.exists(),'Already completed;do not overwrite'
crc=0;uncompressed=0;downloaded=0;digest=hashlib.sha256();decomp=zlib.decompressobj(-15);beg=time.time()
if out.exists():
 with out.open('rb') as prior:
  while True:
   block=prior.read(1024*1024)
   if not block:break
   digest.update(block);downloaded+=len(block);plain=decomp.decompress(block);crc=zlib.crc32(plain,crc);uncompressed+=len(plain)
 print(json.dumps({'resuming_from':downloaded,'prefix_replayed':True}),flush=True)
 assert downloaded<entry['compressed']
with out.open('ab') as f:
 for a in range(start+downloaded,start+entry['compressed'],8*1024*1024):
  b=min(a+8*1024*1024,start+entry['compressed'])-1;got=0
  with fetch(a,b) as response:
   while True:
    block=response.read(1024*1024)
    if not block:break
    f.write(block);digest.update(block);got+=len(block);downloaded+=len(block)
    plain=decomp.decompress(block);crc=zlib.crc32(plain,crc);uncompressed+=len(plain)
  assert got==b-a+1
  state=dict(downloaded=downloaded,expected=entry['compressed'],uncompressed=uncompressed,elapsed=time.time()-beg)
  (P/'H02_evidence_download_progress.json').write_text(json.dumps(state));print(json.dumps(state),flush=True)
plain=decomp.flush();crc=zlib.crc32(plain,crc);uncompressed+=len(plain)
assert decomp.eof and uncompressed==entry['bytes'] and crc==entry['crc'],(decomp.eof,uncompressed,crc)
out.rename(final)
receipt=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),url=url,member=entry['name'],compressed_bytes=downloaded,uncompressed_bytes=uncompressed,crc32=crc,crc_verified=True,sha256_compressed=digest.hexdigest(),file=str(final),elapsed=time.time()-beg)
(P/'H02_evidence_download_receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
