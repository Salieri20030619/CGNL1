from pathlib import Path
import requests,json,datetime,xml.etree.ElementTree as ET,re,hashlib,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2'
studies=[('2026_PMGC','DOI:10.1038/s42003-026-10571-8'),('2024_DSP','DOI:10.1053/j.gastro.2024.08.007'),('2025_PIPS','EXT_ID:PMC12518608')]
def get(item):
 name,query=item;out={'name':name,'query':query}
 try:
  if name=='2025_PIPS':pmc='PMC12518608'
  else:
   a=requests.get('https://www.ebi.ac.uk/europepmc/webservices/rest/search',params={'query':query,'format':'json'},timeout=(15,45));a.raise_for_status();data=a.json();(O/'metadata'/(name+'_publication_lookup.json')).write_text(json.dumps(data,indent=2),encoding='utf-8');hits=data['resultList']['result'];assert hits;out['search_hit']=hits[0];pmc=hits[0].get('pmcid');assert pmc,'No accessible PMC fulltext identifier'
  p=O/'literature'/(pmc+'.xml')
  if not p.exists():
   a=requests.get(f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML',timeout=(15,60));a.raise_for_status();ET.fromstring(a.content);p.write_bytes(a.content)
  root=ET.fromstring(p.read_bytes());out.update(pmc=pmc,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),paragraphs=[],links=[])
  for x in root.findall('.//body//p')+root.findall('.//back//p'):
   s=' '.join(''.join(x.itertext()).split())
   if any(k in s.lower() for k in ['spatial','visium','data availab','deposited','geo','ega','zenodo','figshare','github']):out['paragraphs'].append(s)
  for x in root.iter():
   if x.tag in ['ext-link','media']:out['links'].append({'text':''.join(x.itertext()),'attributes':x.attrib})
  out['accessions']=sorted(set(re.findall(r'GSE\d+|PRJ[A-Z]+\d+|EGA[SD]\d+',' '.join(root.itertext()))));out['success']=True
 except Exception as e:out.update(success=False,error=repr(e))
 (O/'metadata'/(name+'_SPATIAL_QUALIFICATION_120.json')).write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');return out
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:results=list(ex.map(get,studies))
for r in results:print(json.dumps(r,ensure_ascii=False)[:23000])
params={'accession':'PRJNA1146495','result':'read_run','fields':'run_accession,sample_accession,sample_title,library_name,library_strategy,library_layout,fastq_ftp,fastq_bytes,fastq_md5','format':'tsv'}
a=requests.get('https://www.ebi.ac.uk/ena/portal/api/filereport',params=params,timeout=(15,60));a.raise_for_status();(O/'metadata/PRJNA1146495_runs.tsv').write_bytes(a.content);print('ENA_PIPS',a.text[:12000])
