from pathlib import Path
import requests,json,concurrent.futures,xml.etree.ElementTree as ET,re
O=Path.cwd()/'mechanism_investigation_v2'
urls={'GSE53276.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53276&targ=self&form=text&view=full','GSE53276_samples.soft':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53276&targ=gsm&form=text&view=full','DDR2_2016_search.json':'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:10.1038/srep22371&format=json&resultType=core','HPA_CGNL1_gastric.html':'https://www.proteinatlas.org/ENSG00000128849-CGNL1/cancer/stomach+cancer'}
def get(kv):
 n,u=kv;r={'file':n,'url':u}
 try:
  a=requests.get(u,timeout=(15,60));(O/'metadata'/n).write_bytes(a.content);r.update(status=a.status_code,bytes=len(a.content));print(n,a.status_code,len(a.content),flush=True)
 except Exception as e:r['error']=repr(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(get,urls.items()))
(O/'metadata/ACQUISITION_51.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
p=O/'metadata/DDR2_2016_search.json'
if p.exists():
 for v in json.loads(p.read_text())['resultList']['result']:
  if not v.get('pmcid'):continue
  i=v['pmcid'];a=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/'+i+'/fullTextXML',timeout=(15,60));(O/'literature'/(i+'.xml')).write_bytes(a.content);a.raise_for_status();root=ET.fromstring(a.content);t=' '.join(root.itertext());ids=sorted(set(re.findall('GSE[0-9]+',t)));print('DDR2_PAPER',i,ids)
  relevant=[' '.join(e.itertext()) for e in root.iter() if e.tag in ['p','caption'] and any(z in ' '.join(e.itertext()) for z in ['CGNL1','60As6','HSC60','GSE53276'])];(O/'metadata/DDR2_2016_RELEVANT.json').write_text(json.dumps({'pmcid':i,'accessions':ids,'relevant':relevant},indent=2),encoding='utf-8');print('\n'.join(relevant)[:6500])
