import pathlib,urllib.request,subprocess,time
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2');op=urllib.request.build_opener(urllib.request.ProxyHandler({'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}));u='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz';p=O/'metadata/GRCh38.p14.ncbi.gtf.gz';t=time.time()
try:
 with op.open(u,timeout=30) as r,p.open('wb') as f:
  while True:
   b=r.read(1024*1024)
   if not b:break
   f.write(b)
   if time.time()-t>180:raise TimeoutError('download time limit')
 print('GTF',p.stat().st_size,flush=True);print('TRANSFER',subprocess.run(['scp',str(p),'sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2/metadata/']).returncode,flush=True)
except Exception as e:print(repr(e),flush=True)
u='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5528589/fullTextXML'
try:
 with op.open(u,timeout=30) as r:(O/'metadata/GSE50760_article.xml').write_bytes(r.read())
 print('ARTICLE_OK',flush=True)
except Exception as e:print(repr(e),flush=True)
