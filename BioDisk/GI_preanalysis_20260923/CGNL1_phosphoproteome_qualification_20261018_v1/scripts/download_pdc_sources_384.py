from pathlib import Path
import requests,json,hashlib,concurrent.futures
b=Path(__file__).resolve().parents[1];d=b/'raw/PDC_original_384';d.mkdir(exist_ok=True);items=json.loads((b/'metadata/PDC_SOURCE_DOWNLOAD_MANIFEST_384.json').read_text())
def run(it):
 q='{ filesPerStudy(pdc_study_id:'+json.dumps(it['study'])+',file_name:'+json.dumps(it['file_name'])+',offset:0,limit:10){file_id signedUrl{url}} }';r=requests.post('https://pdc.cancer.gov/graphql',json={'query':q},timeout=(15,45));r.raise_for_status();j=r.json();assert not j.get('errors'),j.get('errors');rows=j['data']['filesPerStudy'];row=next(x for x in rows if x['file_id']==it['file_id']);url=row['signedUrl'][0]['url'];p=d/it['file_name'];md=hashlib.md5();sha=hashlib.sha256();n=0
 with requests.get(url,stream=True,timeout=(15,60)) as r:
  r.raise_for_status()
  with p.with_suffix(p.suffix+'.partial').open('wb') as f:
   for chunk in r.iter_content(1048576):f.write(chunk);md.update(chunk);sha.update(chunk);n+=len(chunk)
 assert n==int(it['file_size']) and md.hexdigest()==it['md5sum'];p.with_suffix(p.suffix+'.partial').rename(p);out={'file_id':it['file_id'],'file_name':it['file_name'],'bytes':n,'MD5_verified':True,'sha256':sha.hexdigest()};print(json.dumps(out),flush=True);return out
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:out=list(ex.map(run,items))
(b/'metadata/PDC_SOURCE_DOWNLOAD_COMPLETE_384.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
