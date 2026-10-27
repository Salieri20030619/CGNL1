from pathlib import Path
import urllib.request,json,concurrent.futures
r=Path(r'C:/Users/31937/CGNL1_data_prepare_20260929');op=urllib.request.build_opener(urllib.request.ProxyHandler({}))
urls={
 'osf_root':'https://api.osf.io/v2/nodes/c6tfx/files/osfstorage/',
 'gtex_files':'https://storage.googleapis.com/storage/v1/b/adult-gtex/o?prefix=bulk-qtl/v8/single-tissue-cis-qtl/&maxResults=1000',
 'gtex_old_files':'https://storage.googleapis.com/storage/v1/b/gtex-resources/o?prefix=GTEx_Analysis_v8_QTLs/GTEx_Analysis_v8_eQTL_all_associations/&maxResults=1000',
 'atac_ena':'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP157325&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5&format=tsv',
 'biostudies_files':'https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-13652/info'}
def f(k,u):
 try:
  with op.open(u,timeout=20) as h:b=h.read()
  (r/(k+'.json')).write_bytes(b)
  x=json.loads(b)
  if k=='osf_root':print(k,[(z['attributes'].get('name'),z['attributes'].get('kind'),z['links']) for z in x.get('data',[])])
  elif 'gtex' in k:print(k,[(z['name'],z.get('size')) for z in x.get('items',[]) if any(v in z['name'] for v in ['Liver','Artery','Colon','eQTL.tar'])][:40])
  else:print(k,str(x)[:1800])
 except Exception as e:print(k,'ERROR',str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(lambda kv:f(*kv),urls.items()))
