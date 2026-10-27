from pathlib import Path
import gzip,urllib.request,re,json
D=Path(__file__).resolve().parent.parent/'sc_metadata';D.mkdir(exist_ok=True)
op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}))
for acc,prefix in [('GSE245552','GSE245nnn'),('GSE234804','GSE234nnn'),('GSE164522','GSE164nnn')]:
 url=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/soft/{acc}_family.soft.gz';p=D/(acc+'.soft.gz')
 if not p.exists():p.write_bytes(op.open(url,timeout=60).read())
 s=gzip.open(p,'rt',encoding='utf-8').read();lines=[x for x in s.splitlines() if x.startswith(('!Series_title','!Series_summary','!Series_overall_design','!Series_pubmed_id','!Series_supplementary_file','^SAMPLE','!Sample_title','!Sample_source_name','!Sample_characteristics','!Sample_supplementary_file'))]
 (D/(acc+'_metadata_excerpt.txt')).write_text('\n'.join(lines),encoding='utf-8');print(acc,'\n'+'\n'.join(lines)[:20000])
