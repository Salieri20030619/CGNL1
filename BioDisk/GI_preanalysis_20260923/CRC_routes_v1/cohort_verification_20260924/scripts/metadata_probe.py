from pathlib import Path
import gzip,json,pandas as pd,urllib.request,traceback
B=Path(r'C:\Users\31937\CRC_routes_20260923');D=B/'cohort_verification_20260924';op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897','http':'http://127.0.0.1:7897'}))
m=pd.read_csv(B/'metadata/nonimmune_meta.txt.gz',sep='\t',index_col=0)
print(m.groupby(['patients','organs','cluster']).size().reset_index(name='cells').query("organs=='LCT'").to_string(index=False))
s=gzip.open(B/'metadata/GSE178318.soft.gz','rt',encoding='utf-8').read()
lines=[x for x in s.splitlines() if x.startswith(('!Series_overall_design','!Series_pubmed_id','!Series_supplementary_file','^SAMPLE','!Sample_title','!Sample_characteristics_ch1','!Sample_supplementary_file'))]
(D/'GSE178318_metadata_excerpt.txt').write_text('\n'.join(lines),encoding='utf-8')
print('\n'.join(lines))
urls={'GSE159216_matrix_listing.html':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159216/matrix/','GSE73255_matrix_listing.html':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE73nnn/GSE73255/matrix/','GPL17586_annotation_listing.html':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL17nnn/GPL17586/annot/','GPL10558_annotation_listing.html':'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/'}
for name,url in urls.items():
 try:
  text=op.open(url,timeout=45).read().decode('utf-8');(D/name).write_text(text,encoding='utf-8');print(name,text[:7000])
 except Exception as e:print(name,str(e))
