from pathlib import Path
import json,gzip,re,urllib.request,concurrent.futures,sys
sys.stdout.reconfigure(encoding='utf-8')
O=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3')
op=urllib.request.build_opener(urllib.request.ProxyHandler({'https':'http://127.0.0.1:7897'}))
jobs=[('figshare_'+str(i)+'.json','https://api.figshare.com/v2/articles/'+str(i)) for i in [24803352,24803571,24803616]]
jobs += [('mosaic_liver_zenodo.json','https://zenodo.org/api/records/7737590')]
def run(job):
 name,url=job
 try:
  b=op.open(url,timeout=30).read();(O/'metadata'/name).write_bytes(b);j=json.loads(b)
  files=[{'name':x.get('name',x.get('key')),'bytes':x.get('size'),'url':x.get('download_url',x.get('links',{}).get('self'))} for x in j.get('files',[])]
  return {'id':name,'title':j.get('title',j.get('metadata',{}).get('title')),'files':files,'ok':True}
 except Exception as e:return {'id':name,'ok':False,'error':str(e)}
r=list(concurrent.futures.ThreadPoolExecutor(max_workers=4).map(run,jobs));(O/'results/external_data_availability.json').write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(r,ensure_ascii=False,indent=2))
