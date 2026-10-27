import pathlib,json,gzip,urllib.request,concurrent.futures,subprocess,hashlib
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2');j=json.loads((O/'metadata/infercnvpy_pypi.json').read_text());wheel=next(x for x in j['urls'] if x['filename'].endswith('.whl'));jobs=[('GSE50760_RAW.tar','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50760/suppl/GSE50760_RAW.tar'),(wheel['filename'],wheel['url']),('GSE50760_paper.json','https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:25049118&format=json')]
op=urllib.request.build_opener(urllib.request.ProxyHandler({'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}))
def fetch(v):
 n,u=v
 with op.open(u,timeout=40) as r:data=r.read()
 p=O/'raw'/n;p.write_bytes(data);return {'file':n,'url':u,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(fetch,jobs))
(O/'DOWNLOADS.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
p=O/'metadata/refseq_hg38.txt.gz';print('REFSEQ',p.exists(),p.stat().st_size if p.exists() else 0)
if p.exists():
 with gzip.open(p,'rt') as f:print('REFSEQ_HEAD',f.readline()[:200])
