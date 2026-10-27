from pathlib import Path
import json,requests,concurrent.futures
O=Path.cwd()/'mechanism_investigation_v2';ids=['GSE305534','GSE331489','GSE331488','GSE316970','GSE325526','GSE325523','GSE288730','GSE202503','GSE283522','GSE221512','GSE193566','GSE222859','GSE203611','GSE210132','GSE203612','GSE198832','GSE153562','GSE184462','GSE165659','GSE140182']
def get(i):
 r={'id':i}
 try:
  a=requests.get('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi',params={'acc':i,'targ':'self','form':'text','view':'full'},timeout=(15,50));(O/'metadata'/(i+'_gastric_candidate.soft')).write_bytes(a.content);a.raise_for_status();s=a.content.decode('utf-8');r['fields']=[l for l in s.splitlines() if l.startswith(('!Series_title','!Series_summary','!Series_overall_design','!Series_pubmed_id','!Series_supplementary_file'))];print(i,'\n'.join(l for l in r['fields'] if not l.startswith('!Series_summary'))[:2000],flush=True)
 except Exception as e:r['error']=str(e)
 return r
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(get,ids))
(O/'metadata/GASTRIC_CANDIDATE_CATALOG_33.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8')
