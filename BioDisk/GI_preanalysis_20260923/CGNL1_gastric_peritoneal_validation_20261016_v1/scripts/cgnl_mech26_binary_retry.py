from pathlib import Path
import ftplib,gzip,json,hashlib,datetime
O=Path.cwd()/'mechanism_investigation_v2';names=['GSE223165_fpkm_sample.txt.gz','GSE223165_readcount_genename.txt.gz'];out=[]
for n in names:
 p=O/'raw'/('binary_retry_'+n);r={'file':n,'start':datetime.datetime.now(datetime.timezone.utc).isoformat(),'protocol':'ftplib RETR TYPE I'}
 try:
  with ftplib.FTP('ftp.ncbi.nlm.nih.gov',timeout=60) as ftp:
   ftp.login();ftp.cwd('/geo/series/GSE223nnn/GSE223165/suppl');ftp.voidcmd('TYPE I');r['remote_size']=ftp.size(n)
   with p.open('wb') as f:ftp.retrbinary('RETR '+n,f.write,blocksize=65536)
  b=p.read_bytes();r.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest());assert len(b)==r['remote_size'];raw=gzip.decompress(b);raw.decode('utf-8');r.update(gzip_crc=True,uncompressed_bytes=len(raw),success=True)
 except Exception as e:r.update(success=False,error=str(e))
 print(json.dumps(r),flush=True);out.append(r);(O/'metadata/BINARY_RETRY_26.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
