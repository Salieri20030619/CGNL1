import pathlib,tarfile,gzip,json,urllib.request,subprocess
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2');z=tarfile.open(O/'raw/GSE50760_RAW.tar');print('FILES',z.getnames()[:5])
f=gzip.GzipFile(fileobj=z.extractfile(z.getmembers()[0]));print(f.read(800).decode())
j=json.loads((O/'raw/GSE50760_paper.json').read_text());print('PAPER',j['resultList']['result'][0].get('pmcid'))
u='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz';p=O/'metadata/gencode.v48.annotation.gtf.gz';op=urllib.request.build_opener(urllib.request.ProxyHandler({'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}))
try:
 with op.open(u,timeout=45) as r,p.open('wb') as f:
  while True:
   b=r.read(1024*1024)
   if not b:break
   f.write(b)
 print('GTF_BYTES',p.stat().st_size,flush=True)
 with gzip.open(p,'rt') as f:print(f.readline().strip())
 print('GTF_SCP',subprocess.run(['scp',str(p),'sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2/metadata/']).returncode)
except Exception as e:print(type(e).__name__,str(e),flush=True)
