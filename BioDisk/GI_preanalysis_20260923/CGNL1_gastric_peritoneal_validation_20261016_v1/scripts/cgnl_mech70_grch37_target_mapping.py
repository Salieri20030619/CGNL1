from pathlib import Path
import requests,json,datetime,hashlib
O=Path.cwd()/'mechanism_investigation_v2';out=[]
for gene in ['CGNL1','CGN']:
 u='https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/'+gene
 r={'gene':gene,'url':u,'assembly_expected':'GRCh37','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  a=requests.get(u,params={'expand':1,'content-type':'application/json'},headers={'Accept':'application/json'},timeout=(15,60));p=O/'metadata'/('GRCh37_'+gene+'_expanded.json');p.write_bytes(a.content);a.raise_for_status();d=a.json();assert d['assembly_name']=='GRCh37' and d['display_name']==gene;r.update(gene_id=d['id'],transcripts=[{'id':t['id'],'version':t.get('version'),'biotype':t.get('biotype'),'start':t['start'],'end':t['end']} for t in d.get('Transcript',[])],sha256=hashlib.sha256(a.content).hexdigest(),success=True)
 except Exception as e:r.update(success=False,error=repr(e))
 out.append(r);print(json.dumps(r),flush=True)
(O/'metadata/GRCH37_MAPPING_RECEIPTS.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
