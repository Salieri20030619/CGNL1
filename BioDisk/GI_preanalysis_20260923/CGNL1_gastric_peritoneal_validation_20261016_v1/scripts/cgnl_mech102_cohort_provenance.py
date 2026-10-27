from pathlib import Path
import gzip,json,re,requests
O=Path.cwd()/'mechanism_investigation_v2';B=O.parent.parent;out={}
for name,p in [('GSE314812',B/'01_metadata/GSE314812.soft.gz'),('GSE237876',O/'metadata/GSE237876.soft')]:
 b=p.read_bytes();s=(gzip.decompress(b) if b[:2]==b'\x1f\x8b' else b).decode('utf-8');keys=['title','summary','overall_design','contributor','contact_institute','contact_name','pubmed_id','relation','submission_date'];fields=[l for l in s.splitlines() if any(l.startswith('!Series_'+k+' = ') for k in keys)];out[name]={'file':str(p),'series_fields':fields};print(name,'\n'.join(fields)[:6200])
u='https://www.ebi.ac.uk/europepmc/webservices/rest/search';a=requests.get(u,params={'query':'EXT_ID:38498378','format':'json','resultType':'core'},timeout=(15,45));a.raise_for_status();(O/'metadata/GSE237876_paper.json').write_bytes(a.content);print('ROUTE_PAPER',[(r.get('title'),r.get('pmcid'),r.get('doi')) for r in a.json()['resultList']['result']]);out['route_paper_records']=a.json()['resultList']['result'];(O/'metadata/COHORT_PROVENANCE_102.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
