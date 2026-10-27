import pathlib,urllib.request,json,gzip,re
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2');op=urllib.request.build_opener(urllib.request.ProxyHandler({'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}));jobs=[('GSE14297.soft.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14297/soft/GSE14297_family.soft.gz'),('GSE14297.matrix.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14297/matrix/GSE14297_series_matrix.txt.gz')]
for name,url in jobs:
 try:
  with op.open(url,timeout=35) as r:(O/'metadata'/name).write_bytes(r.read())
  print(name,(O/'metadata'/name).stat().st_size,flush=True)
 except Exception as e:print(name,repr(e),flush=True)
s=gzip.open(O/'metadata/GSE14297.soft.gz','rt',encoding='utf-8').read();print('\n'.join(x for x in s.splitlines() if x.startswith(('!Series_title','!Series_platform','!Series_pubmed','!Series_supplementary'))))
for c in s.split('^SAMPLE = ')[1:5]:print(c.splitlines()[0],';'.join(x for x in c.splitlines() if x.startswith(('!Sample_title','!Sample_characteristics'))))
