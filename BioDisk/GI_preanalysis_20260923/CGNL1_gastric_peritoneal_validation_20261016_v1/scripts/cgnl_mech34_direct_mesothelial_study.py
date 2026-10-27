from pathlib import Path
import requests,json,re,xml.etree.ElementTree as ET
O=Path.cwd()/'mechanism_investigation_v2'
q='TITLE:"Mesothelial cells with mesenchymal features"';a=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',params={'query':q,'format':'json','resultType':'core'},timeout=(15,45));(O/'metadata/direct_meso_study_search.json').write_bytes(a.content);a.raise_for_status();r=a.json()['resultList']['result'];print([(x.get('title'),x.get('pmcid'),x.get('doi')) for x in r])
out=[]
for x in r:
 if not x.get('pmcid'):continue
 i=x['pmcid'];b=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+i+'/fullTextXML',timeout=(15,50));(O/'literature'/(i+'.xml')).write_bytes(b.content);b.raise_for_status();root=ET.fromstring(b.content);t=' '.join(root.itertext());ids=sorted(set(re.findall(r'GSE\d+|PRJ\w+\d+|DRA\d+',t)));sections=[]
 for sec in root.iter():
  title=sec.find('title')
  if sec.tag in ['sec','notes'] and title is not None and any(z in ''.join(title.itertext()).lower() for z in ['rna','availability']):sections.append(' '.join(sec.itertext()))
 out.append({'pmcid':i,'ids':ids,'sections':sections});print(i,ids,sections)
(O/'metadata/DIRECT_MESOTHELIAL_STUDY_34.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
