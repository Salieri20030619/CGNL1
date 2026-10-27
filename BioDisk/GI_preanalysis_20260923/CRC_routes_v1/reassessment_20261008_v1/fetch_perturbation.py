from pathlib import Path
import os,json,hashlib,datetime,gzip,zipfile,concurrent.futures
for k in list(os.environ):
 if 'proxy' in k.lower():os.environ.pop(k,None)
os.environ.update(NO_PROXY='*',no_proxy='*')
import requests
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/reassessment_20261008_v1');D=R/'perturbation_inputs';D.mkdir(exist_ok=True)
base='https://ftp.ncbi.nlm.nih.gov/geo/series/'
files={
 'GSE163568_family.soft.gz':base+'GSE163nnn/GSE163568/soft/GSE163568_family.soft.gz',
 'GSE237653_family.soft.gz':base+'GSE237nnn/GSE237653/soft/GSE237653_family.soft.gz',
 'GSE163568_TLA_HUVEC_SE_counts.txt.gz':base+'GSE163nnn/GSE163568/suppl/GSE163568_TLA_HUVEC_SE_counts.txt.gz',
 'GSE163568_TLA-HRECcounts.txt.gz':base+'GSE163nnn/GSE163568/suppl/GSE163568_TLA-HRECcounts.txt.gz',
 'GSE163568_DEG_HUVEC_TLA.xlsx':base+'GSE163nnn/GSE163568/suppl/GSE163568_DEG_HUVEC_TLA.xlsx',
 'GSE163568_DEG_HREC_TLA.xlsx':base+'GSE163nnn/GSE163568/suppl/GSE163568_DEG_HREC_TLA.xlsx',
 'GSE237653_2023_Teichmann_NCoR1_RNAseq_NormCounts.xlsx':base+'GSE237nnn/GSE237653/suppl/GSE237653_2023_Teichmann_NCoR1_RNAseq_NormCounts.xlsx'
}
def fetch(kv):
 name,url=kv;p=D/name;rec={'name':name,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};s=requests.Session();s.trust_env=False
 try:
  q=s.get(url,timeout=(15,60),stream=True);rec['status']=q.status_code;q.raise_for_status();part=p.with_suffix(p.suffix+'.part');n=0;h=hashlib.sha256()
  with part.open('wb') as f:
   for chunk in q.iter_content(1024*1024):
    n+=len(chunk)
    if n>80*1024*1024:raise ValueError('processed_file_exceeds80MiB_safety_limit')
    f.write(chunk);h.update(chunk)
  if name.endswith('.gz'):
   with gzip.open(part,'rb') as f:
    while f.read(1024*1024):pass
  if name.endswith('.xlsx'):
   with zipfile.ZipFile(part) as z:assert z.testzip() is None
  part.rename(p);rec.update(bytes=n,sha256=h.hexdigest(),integrity=True)
 except Exception as e:rec['error']=repr(e)
 return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 out=list(pool.map(fetch,files.items()))
(D/'FETCH_MANIFEST.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
