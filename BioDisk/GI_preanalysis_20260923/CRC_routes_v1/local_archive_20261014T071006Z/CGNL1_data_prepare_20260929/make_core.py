from pathlib import Path
import urllib.request,gzip,json,re,traceback
r=Path(r'C:/Users/31937/CGNL1_data_prepare_20260929');m=Path(r'C:/Users/31937/CRC_routes_20260923/expanded_screen_20260925/metadata')
try:
 jobs=[]
 for l in gzip.decompress((m/'GSE234804.soft.gz').read_bytes()).decode().splitlines():
  if l.startswith('!Sample_supplementary_file') and ' = ' in l:
   u=l.split(' = ',1)[1].replace('ftp://','https://')
   if u.startswith('https://'):jobs.append({'dataset':'GSE234804','name':u.rsplit('/',1)[1],'url':u,'role':'human_scRNA_processed'})
 op=urllib.request.build_opener(urllib.request.ProxyHandler({}));u='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131418/matrix/'
 with op.open(u,timeout=30) as h:b=h.read()
 (r/'GSE131418_index.html').write_bytes(b)
 print(b.decode()[:1400])
 for name in sorted(set(re.findall(r'href="([^"]+series_matrix[^"?]*\.gz)"',b.decode()))):
  jobs.append({'dataset':'GSE131418','name':name.rsplit('/',1)[-1],'url':urllib.parse.urljoin(u,name),'role':'bulk_expression'})
 (r/'wave01_core.json').write_text(json.dumps(jobs,indent=2));print('JOBS',len(jobs))
except Exception:traceback.print_exc()
