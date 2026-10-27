from pathlib import Path
import json,re,requests,concurrent.futures,datetime
O=Path.cwd()/'mechanism_investigation_v2';ids=set()
for p in sorted((O/'metadata').glob('search_*.json')):
 try:
  d=json.loads(p.read_text());q=d.get('esearchresult',{});series=['GSE'+i[3:] for i in q.get('idlist',[]) if i.startswith('200')];ids.update(series);print(p.name,'series',series)
 except Exception as e:print(p.name,str(e))
for n in ['GSE237876.soft','GSE165211.soft']:
 s=(O/'metadata'/n).read_text();print(n,'\n'.join(l for l in s.splitlines() if l.startswith(('!Series_title','!Series_overall_design','!Series_supplementary_file','!Series_pubmed_id'))))
s=(O/'metadata/GSE237876_samples.soft').read_text();samples=[]
for b in s.split('^SAMPLE = ')[1:]:
 lines=b.splitlines();d={'accession':lines[0]}
 for k in ['title','characteristics_ch1','source_name_ch1']:
  d[k]=[l.split(' = ',1)[1] for l in lines if l.startswith('!Sample_'+k+' = ')]
 samples.append(d)
(O/'metadata/GSE237876_SAMPLE_MAP.json').write_text(json.dumps(samples,indent=2),encoding='utf-8');print('ROUTE_SAMPLES',json.dumps(samples)[:7000])
def get(i):
 p=O/'metadata'/(i+'_candidate.soft')
 try:
  a=requests.get('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi',params={'acc':i,'targ':'self','form':'text','view':'full'},timeout=(15,40));p.write_bytes(a.content);a.raise_for_status();s=a.text;v=[l for l in s.splitlines() if l.startswith(('!Series_title','!Series_overall_design','!Series_supplementary_file'))];print(i,'\n'.join(v),flush=True);return {'id':i,'fields':v}
 except Exception as e:return {'id':i,'error':str(e)}
# Restrict to focused mesothelial queries, not broad gastric all-series matches.
focus=set()
for key in ['mesothelial_gastric','mesothelial_perturbation']:
 p=O/'metadata'/('search_'+key+'.json')
 if p.exists():focus.update('GSE'+i[3:] for i in json.loads(p.read_text()).get('esearchresult',{}).get('idlist',[]) if i.startswith('200'))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:r=list(ex.map(get,sorted(focus)))
(O/'metadata/FOCUSED_CANDIDATES_15.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
