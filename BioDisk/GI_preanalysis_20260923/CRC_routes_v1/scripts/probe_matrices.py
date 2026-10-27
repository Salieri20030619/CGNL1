import pathlib,gzip,importlib.util,tarfile,json
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
with gzip.open(O/'raw/nonimmune_counts.txt.gz','rt') as f:
 for i in range(3):
  s=f.readline().rstrip('\n').split('\t');print('ROW',i,'FIELDS',len(s),'HEAD',s[:8],flush=True)
print('infercnvpy',bool(importlib.util.find_spec('infercnvpy')))
for p in O.glob('raw/M-ST-*.tar.gz'):
 try:
  with tarfile.open(p) as t:print(p.name,[(x.name,x.size) for x in t.getmembers()][:16])
 except Exception as e:print(p.name,type(e).__name__)
