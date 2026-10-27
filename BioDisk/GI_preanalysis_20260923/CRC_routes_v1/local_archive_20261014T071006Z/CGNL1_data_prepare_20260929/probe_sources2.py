import urllib.request,json,concurrent.futures
from pathlib import Path
r=Path(r'C:/Users/31937/CGNL1_data_prepare_20260929');op=urllib.request.build_opener(urllib.request.ProxyHandler({}))
urls={
'gtex_allpairs':'https://storage.googleapis.com/storage/v1/b/adult-gtex/o?prefix=bulk-qtl/v8/single-tissue-cis-qtl/eQTL-Tissue-Specific-All-SNP-Gene-Associations/&maxResults=1000',
'finngen_crc':'https://storage.googleapis.com/storage/v1/b/finngen-public-data-r12/o?prefix=summary_stats/finngen_R12_C3_COLORECTAL&maxResults=50',
'ena_atac_probe':'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERR12397437&result=read_run&fields=run_accession,study_accession,fastq_ftp,fastq_bytes,fastq_md5,submitted_ftp,submitted_bytes,submitted_md5&format=tsv'}
def fetch(k,u):
 try:
  with op.open(u,timeout=30) as f:b=f.read()
  (r/(k+'.json' if k!='ena_atac_probe' else k+'.tsv')).write_bytes(b)
  if k=='ena_atac_probe':print(k,b.decode()[:1600]);return
  x=json.loads(b);print(k,[(v['name'],v.get('size')) for v in x.get('items',[]) if k=='finngen_crc' or any(s in v['name'] for s in ['Liver','Artery','Colon'])])
 except Exception as e:print(k,str(e))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:list(ex.map(lambda kv:fetch(*kv),urls.items()))
