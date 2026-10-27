import json, urllib.request, urllib.parse
from pathlib import Path
from datetime import datetime, timezone
urllib.request.install_opener(urllib.request.build_opener(urllib.request.ProxyHandler({})))
ROOT=Path(__file__).resolve().parent
queries={
 'target_all':'(CGNL1 OR paracingulin OR "cingulin-like 1" OR "junction-associated coiled-coil protein")',
 'target_cancer':'(CGNL1 OR paracingulin OR JACOP) AND (cancer OR carcinoma OR metastasis OR tumor)',
 'axis_crc':'(CAMSAP3 OR PLEKHA7 OR TIAM1 OR ARHGEF2) AND (colorectal OR colon) AND (metastasis OR polarity OR adhesion)',
 'target_geo':'(CGNL1 OR paracingulin) AND (GEO OR RNA-seq OR transcriptome OR proteomics)',
}
records={};audit=[]
for key,query in queries.items():
 url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':'('+query+') AND FIRST_PDATE:[1900-01-01 TO 2026-10-12]','format':'json','resultType':'core','pageSize':1000})
 try:
  with urllib.request.urlopen(url,timeout=50) as r:data=json.load(r)
  (ROOT/(key+'.json')).write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
  audit.append({'query':query,'url':url,'hitCount':data.get('hitCount'),'retrieved':len(data['resultList']['result']),'utc':datetime.now(timezone.utc).isoformat()})
  for row in data['resultList']['result']:
   ident=row.get('source','')+':'+row.get('id','')
   if ident not in records:records[ident]={**row,'searches':[]}
   records[ident]['searches'].append(key)
  print(key,data.get('hitCount'),flush=True)
 except Exception as e:audit.append({'query':query,'error':str(e)});print(key,str(e),flush=True)
(ROOT/'search_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
(ROOT/'all_records.json').write_text(json.dumps(list(records.values()),ensure_ascii=False,indent=2),encoding='utf-8')
for r in records.values():
 title=r.get('title','')
 if any(s in title.lower() for s in ['cgnl1','paracingulin','cingulin-like','camsap3']):
  print(r.get('id'),r.get('pmcid'),r.get('pubYear'),title,flush=True)
