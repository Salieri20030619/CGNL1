from pathlib import Path
import json,requests,concurrent.futures,xml.etree.ElementTree as ET,re
O=Path.cwd()/'mechanism_investigation_v2'
s=(O/'metadata/GSE308231_family.soft').read_text(encoding='utf-8')
lines=[l for l in s.splitlines() if l.startswith('!Sample_') and any(x in l for x in ['title =','source_name','characteristics','supplementary_file','data_processing'])]
print('TISSUE_METADATA','\n'.join(lines));(O/'metadata/GSE308231_sample_fields.txt').write_text('\n'.join(lines),encoding='utf-8')
p=json.loads((O/'metadata/GSE308231_paper.json').read_text(encoding='utf-8'));res=p['resultList']['result'];print('NEW_TISSUE_PAPER',[(r.get('title'),r.get('pmcid'),r.get('abstractText')) for r in res])
print('CGNL_GEO_QUERY',(O/'metadata/cgnl_geo_search.json').read_text(encoding='utf-8'))
for name in ['ascites_europepmc','cgnl_ko_europepmc']:
 root=ET.fromstring((O/'literature'/(name+'.response')).read_bytes())
 for sec in root.iter():
  if sec.tag in ['notes','sec'] and ('availability' in (sec.attrib.get('sec-type','')+sec.attrib.get('notes-type','')).lower()):print(name,'AVAILABILITY',' '.join(sec.itertext()))
urls={'GSE228955.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228955&targ=self&form=text&view=full','GSE228955_files.html':'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228955/suppl/','HRA002712_retry.html':'https://ngdc.cncb.ac.cn/gsa-human/browse/HRA002712','ascites_supplement_package.zip':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9929081/supplementaryFiles'}
for r in res:
 if r.get('pmcid'):urls[r['pmcid']+'.xml']='https://www.ebi.ac.uk/europepmc/webservices/rest/'+r['pmcid']+'/fullTextXML'
def fetch(k,u):
 r={'file':k,'url':u}
 try:
  a=requests.get(u,timeout=(15,60));(O/'metadata'/k).write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content));print('FETCH',k,a.status_code,len(a.content),flush=True)
 except Exception as e:r['error']=str(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:r=list(ex.map(lambda kv:fetch(*kv),urls.items()))
(O/'metadata/ACQUISITION_06.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
