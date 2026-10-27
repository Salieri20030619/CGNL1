from pathlib import Path
import io,zlib,csv,json,hashlib,collections,time
P=Path(r'C:\Users\31937\cgnl_transfer\restart_20261024');f=P/'H02_evidence.deflate'
receipt=json.loads((P/'H02_evidence_download_receipt.json').read_text());assert hashlib.sha256(f.read_bytes()).hexdigest()==receipt['sha256_compressed']
class Inflated(io.RawIOBase):
 def __init__(self,path):self.f=path.open('rb');self.z=zlib.decompressobj(-15);self.buf=b'';self.end=False;self.crc=0;self.count=0
 def readable(self):return True
 def readinto(self,b):
  while not self.buf and not self.end:
   c=self.f.read(1024*1024)
   if c:self.buf=self.z.decompress(c)
   else:self.buf=self.z.flush();self.end=True
   self.crc=zlib.crc32(self.buf,self.crc);self.count+=len(self.buf)
  n=min(len(b),len(self.buf));b[:n]=self.buf[:n];self.buf=self.buf[n:];return n
stream=Inflated(f);csv.field_size_limit(20000000);r=csv.DictReader(io.TextIOWrapper(io.BufferedReader(stream),encoding='utf-8'),delimiter='\t')
records=[];total=0;beg=time.time()
for row in r:
 total+=1
 if 'Q0VF96' in row.get('Proteins','').split(';') or 'CGNL1' in row.get('Gene names','').split(';'):records.append(row)
assert stream.crc==receipt['crc32'] and stream.count==receipt['uncompressed_bytes']
out=P/'H02_CGNL1_evidence_rows.json';out.write_text(json.dumps(records,indent=2),encoding='utf-8')
summary=dict(total_evidence_rows=total,target_rows=len(records),types=dict(collections.Counter(x['Type'] for x in records)),crc_verified=True,uncompressed_bytes=stream.count,sha256_target=hashlib.sha256(out.read_bytes()).hexdigest(),elapsed=time.time()-beg)
(P/'H02_target_evidence_receipt.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
