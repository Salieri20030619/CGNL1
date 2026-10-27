from pathlib import Path
import urllib.request,json,xml.etree.ElementTree as ET,concurrent.futures,re
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/mainline_search_20261012_v1'); F=R/'fulltext';F.mkdir(exist_ok=True)
ids=['PMC2555940','PMC3514035','PMC10184829','PMC10202830','PMC11513637','PMC12238879','PMC10417749','PMC8227492','PMC11697589','PMC12683414','PMC13163822','PMC4101344','PMC4159652','PMC5852532','PMC12278653','PMC13390483','PMC7351134','PMC13380819','PMC3089516','PMC4072573']
def fetch(i):
 u='https://www.ebi.ac.uk/europepmc/webservices/rest/'+i+'/fullTextXML'
 try:
  op=urllib.request.build_opener(urllib.request.ProxyHandler({}))
  with op.open(u,timeout=35) as r: raw=r.read()
  tree=ET.fromstring(raw);(F/(i+'.xml')).write_bytes(raw)
  paras=[' '.join(''.join(p.itertext()).split()) for p in tree.findall('.//body//p')]
  (F/(i+'.txt')).write_text('\n\n'.join(paras),encoding='utf-8')
  row={'pmcid':i,'url':u,'bytes':len(raw),'paragraphs':len(paras),'status':'ok'}
 except Exception as e:row={'pmcid':i,'url':u,'status':'failed','error':str(e)}
 print(i,row['status'],flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:out=list(pool.map(fetch,ids))
(R/'fulltext_fetch_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print('FULLTEXT_COMPLETE',len(out),sum(x['status']=='ok' for x in out),flush=True)
