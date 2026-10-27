from pathlib import Path
import os, json, datetime, hashlib, gzip, re, concurrent.futures
for key in list(os.environ):
 if 'proxy' in key.lower(): os.environ.pop(key,None)
os.environ.update(NO_PROXY='*',no_proxy='*')
import requests
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
R=B/'feasibility_audit_20261008_v1';D=R/'metadata';D.mkdir(exist_ok=True)
urls={
 'GSE294385':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294385&targ=all&form=text&view=full',
 'GSE217414':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217414&targ=all&form=text&view=full',
 'GSE335552':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335552&targ=all&form=text&view=full',
 'Zenodo19043057':'https://zenodo.org/api/records/19043057',
 'Zenodo21431802':'https://zenodo.org/api/records/21431802',
 '10x_spatial_outputs':'https://www.10xgenomics.com/support/software/space-ranger/2.1/analysis/outputs/spatial-outputs'
}
def fetch(item):
 key,url=item; sess=requests.Session();sess.trust_env=False
 record={'key':key,'url':url,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  q=sess.get(url,timeout=(15,35));record.update(status=q.status_code,final_url=q.url,bytes=len(q.content),sha256=hashlib.sha256(q.content).hexdigest());(D/(key+'.response')).write_bytes(q.content)
  if key.startswith('Zenodo') and q.ok:
   j=q.json();m=j.get('metadata',{});record.update(title=m.get('title'),access=m.get('access_right'),files=[{'key':f.get('key'),'size':f.get('size')} for f in j.get('files',[])],latest=j.get('links',{}).get('latest_html'),notes=m.get('notes'))
  elif key.startswith('GSE') and q.ok:
   good=q.text.startswith('^SERIES') or '^SAMPLE' in q.text;record['valid_soft']=good
   if good:
    lines=[s for s in q.text.splitlines() if s.startswith(('^SAMPLE','!Sample_title','!Sample_source','!Sample_characteristics','!Sample_supplementary_file','!Sample_data_processing','!Series_title','!Series_overall_design','!Series_supplementary_file'))]
    (D/(key+'_relevant.txt')).write_text('\n'.join(lines));record['samples']=q.text.count('^SAMPLE');record['supplementary_lines']=sum('supplementary_file' in s for s in lines)
   else:record['response_excerpt']=q.text[:300]
  elif key=='10x_spatial_outputs':
   record['mentions_theoretical65']=bool(re.search(r'theoretical.{0,80}65|65.{0,80}spot',q.text,re.I));i=q.text.find('spot_diameter_fullres');record['definition_excerpt']=q.text[i:i+1500]
 except Exception as e:record['error']=repr(e)
 return record
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(fetch,urls.items()))
(D/'DIRECT_METADATA_CHECK.json').write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
