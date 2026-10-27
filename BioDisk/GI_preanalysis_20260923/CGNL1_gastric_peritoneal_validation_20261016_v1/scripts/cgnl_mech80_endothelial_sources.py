from pathlib import Path
import requests,json,xml.etree.ElementTree as ET,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2';ids=['PMC5852532','PMC12238879']
def get(i):
 r={'pmcid':i,'url':'https://www.ebi.ac.uk/europepmc/webservices/rest/'+i+'/fullTextXML'}
 try:
  a=requests.get(r['url'],timeout=(15,60));(O/'literature'/(i+'.xml')).write_bytes(a.content);a.raise_for_status();root=ET.fromstring(a.content);r['abstracts']=[' '.join(x.itertext()) for x in root.iter('abstract')];r['availability']=[' '.join(x.itertext()) for x in root.iter('sec') if x.find('title') is not None and 'availability' in ''.join(x.find('title').itertext()).lower()];r['direct_phenotype']=[' '.join(x.itertext()) for x in root.iter('p') if ('Cgnl1' in ' '.join(x.itertext()) or 'CGNL1' in ' '.join(x.itertext())) and any(q in ' '.join(x.itertext()).lower() for q in ['notch','pericyte','pericytes','rescue','retina'])];print(i,json.dumps({'abstracts':r['abstracts'],'availability':r['availability']},ensure_ascii=False),flush=True)
 except Exception as e:r['error']=repr(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:rs=list(ex.map(get,ids))
(O/'metadata/ENDOTHELIAL_FUNCTIONAL_SOURCES_80.json').write_text(json.dumps(rs,ensure_ascii=False,indent=2),encoding='utf-8')
import pandas as pd
p=O.parent/'results/bulk_patient_deltas.tsv';d=pd.read_csv(p,sep='\t');print('PAIRED_COLUMNS',list(d.columns));print('PAIRED_DESIGN_FIRST_ROWS',d.head(2).to_string(index=False))
