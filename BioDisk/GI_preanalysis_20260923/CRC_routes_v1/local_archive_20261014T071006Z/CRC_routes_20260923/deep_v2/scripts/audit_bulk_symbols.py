import pathlib,tarfile,gzip,pandas as pd,urllib.request,subprocess,json
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2');t=tarfile.open(O/'raw/GSE50760_RAW.tar');f=pd.read_csv(gzip.GzipFile(fileobj=t.extractfile(t.getmembers()[0])),sep='\t');print('DUP',f[f.iloc[:,0].duplicated(keep=False)].to_string(index=False));print('SHAPE',f.shape,'NAN',f.isna().sum().to_dict())
u='https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz';op=urllib.request.build_opener(urllib.request.ProxyHandler({'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}));p=O/'metadata/Homo_sapiens.gene_info.gz'
with op.open(u,timeout=35) as r,p.open('wb') as z:
 while True:
  b=r.read(1024*1024)
  if not b:break
  z.write(b)
print('GENE_INFO',p.stat().st_size,flush=True);print('TRANSFER',subprocess.run(['scp',str(p),'sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2/metadata/']).returncode,flush=True)
