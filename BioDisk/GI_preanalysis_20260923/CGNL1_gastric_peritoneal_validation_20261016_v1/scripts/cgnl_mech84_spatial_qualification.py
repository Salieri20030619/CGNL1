from pathlib import Path
import requests,json
O=Path.cwd()/'mechanism_investigation_v2';u='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251950&targ=gsm&form=text&view=full';a=requests.get(u,timeout=(15,60));p=O/'metadata/GSE251950_samples.soft';p.write_bytes(a.content);a.raise_for_status();s=a.content.decode('utf-8');rows=[]
for b in s.split('^SAMPLE = ')[1:]:
 ls=b.splitlines();d={'accession':ls[0]}
 for k in ['title','source_name_ch1','characteristics_ch1','description','supplementary_file_1','data_processing']:d[k]=[l.split(' = ',1)[1] for l in ls if l.startswith('!Sample_'+k+' = ')]
 rows.append(d);print(json.dumps({k:v for k,v in d.items() if k!='data_processing'},ensure_ascii=False))
(O/'metadata/SPATIAL_QUALIFICATION_84.json').write_text(json.dumps({'source':u,'samples':rows,'status':'samplemetadataqualification;noexpressionanalysisyet'},ensure_ascii=False,indent=2),encoding='utf-8')
