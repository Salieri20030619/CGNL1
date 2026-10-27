from pathlib import Path
import requests,json,xml.etree.ElementTree as ET
O=Path.cwd()/'mechanism_investigation_v2';ids=['PMC11024677'];a=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',params={'query':'EXT_ID:41882239','format':'json','resultType':'core'},timeout=(15,45));a.raise_for_status();(O/'metadata/GSE314812_publication_record.json').write_bytes(a.content)
for r in a.json()['resultList']['result']:
 print('PUBLICATION',r.get('title'),r.get('pmcid'))
 if r.get('pmcid'):ids.append(r['pmcid'])
out=[]
for i in ids:
 a=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+i+'/fullTextXML',timeout=(15,60));p=O/'literature'/(i+'.xml');p.write_bytes(a.content);a.raise_for_status();root=ET.fromstring(a.content);pars=[]
 for e in root.iter('p'):
  t=' '.join(e.itertext())
  if any(k in t.lower() for k in ['collected','institutional review','enrolled','recruited','consent','patient cohort','26 patients','14 patients']):pars.append(t)
 out.append({'pmcid':i,'clinical_paragraphs':pars});print(i,'\n'.join(pars)[:11500])
(O/'metadata/CLINICAL_PROVENANCE_103.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
