from pathlib import Path
import requests,json,xml.etree.ElementTree as ET
O=Path.cwd()/'mechanism_investigation_v2';u='https://www.ebi.ac.uk/europepmc/webservices/rest/search';a=requests.get(u,params={'query':'EXT_ID:39580151','format':'json','resultType':'core'},timeout=(15,45));a.raise_for_status();(O/'metadata/GSE251950_paper.json').write_bytes(a.content);out=[]
for r in a.json()['resultList']['result']:
 print('PAPER',r.get('title'),r.get('pmcid'),r.get('doi'),r.get('abstractText'))
 if not r.get('pmcid'):continue
 i=r['pmcid'];a=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+i+'/fullTextXML',timeout=(15,50));(O/'literature'/(i+'.xml')).write_bytes(a.content);a.raise_for_status();root=ET.fromstring(a.content);pars=[' '.join(e.itertext()) for e in root.iter('p') if any(k in ' '.join(e.itertext()) for k in ['GC6','peritoneal','Peritoneal','251950'])];out.append({'pmcid':i,'paragraphs':pars});print('SITE_PARAGRAPHS','\n'.join(pars)[:7500])
(O/'metadata/SPATIAL_SITE_QUALIFICATION_86.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
