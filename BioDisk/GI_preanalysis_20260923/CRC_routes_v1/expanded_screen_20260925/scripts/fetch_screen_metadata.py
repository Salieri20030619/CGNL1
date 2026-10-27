from pathlib import Path
import urllib.request,json,gzip,concurrent.futures,xml.etree.ElementTree as ET
A=Path(__file__).resolve().parent.parent;M=A/'metadata';M.mkdir(exist_ok=True)
op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897','http':'http://127.0.0.1:7897'}))
jobs=[]
for g in ['GSE206552','GSE254364','GSE203612','GSE224235','GSE234804','GSE315534','GSE131418','GSE41258']:
 stem=g[:-3]+'nnn';jobs.append((g+'.soft.gz',f'https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}/{g}/soft/{g}_family.soft.gz'))
for acc in ['E-MTAB-13655','E-MTAB-13651','E-MTAB-13652']:
 jobs.append((acc+'.json','https://www.ebi.ac.uk/biostudies/api/v1/studies/'+acc))
for pmc in ['PMC10352575','PMC12542284','PMC9811165']:
 jobs.append((pmc+'.xml',f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML'))
jobs += [('mendeley_files.json','https://data.mendeley.com/public-api/datasets/yd7bb3shn5/versions/1/files'),('stanford.html','https://dna-discovery.stanford.edu/research/datasets/'),('CRC_metastasis_tree.json','https://api.github.com/repos/EfremovaLab/CRC_metastasis/git/trees/main?recursive=1')]
def run(j):
 name,url=j;p=M/name
 try:
  if not p.exists():
   req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'})
   with op.open(req,timeout=40) as r:data=r.read(60000000)
   assert len(data)<60000000;p.write_bytes(data)
  data=p.read_bytes();rec={'file':name,'url':url,'bytes':len(data),'status':'downloaded'}
  if name.endswith('.soft.gz'):
   lines=gzip.decompress(data).decode('utf8',errors='replace').splitlines();keep=[l for l in lines if l.startswith(('^SAMPLE','!Series_title','!Series_pubmed_id','!Series_summary','!Series_overall_design','!Series_supplementary_file','!Sample_title','!Sample_source_name','!Sample_characteristics','!Sample_supplementary_file'))];(M/(name[:-8]+'_excerpt.txt')).write_text('\n'.join(keep),encoding='utf8')
  elif name.endswith('.xml'):
   root=ET.fromstring(data);(M/(name[:-4]+'.txt')).write_text('\n'.join(' '.join(n.itertext()) for n in root.iter() if n.tag in ['p','title','table']),encoding='utf8')
  print(name,len(data),'OK',flush=True);return rec
 except Exception as e:print(name,str(e),flush=True);return {'file':name,'url':url,'status':'error','error':str(e)}
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:res=list(ex.map(run,jobs))
(A/'METADATA_FETCH.json').write_text(json.dumps(res,indent=2));print('COMPLETE',flush=True)
