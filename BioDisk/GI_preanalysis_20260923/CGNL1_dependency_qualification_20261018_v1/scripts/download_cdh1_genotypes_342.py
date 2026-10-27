from pathlib import Path
import requests,json,hashlib,zipfile,datetime,concurrent.futures
b=Path(__file__).resolve().parents[1];d=b/'raw/cdh1_genotype_341';d.mkdir(exist_ok=True)
items=[('mutations_summary_20260724.zip',584662),('cnv_summary_20260804.zip',9098915)]
def fetch(item):
 name,size=item;url='https://cog.sanger.ac.uk/cmp/download/'+name;r=requests.get(url,timeout=(15,60));r.raise_for_status();assert len(r.content)==size;(d/name).write_bytes(r.content)
 with zipfile.ZipFile(d/name) as z:assert z.testzip() is None;members=[{'name':i.filename,'size':i.file_size} for i in z.infolist()]
 return {'name':name,'url':url,'bytes':size,'sha256':hashlib.sha256(r.content).hexdigest(),'members':members,'ZIP_CRC_pass':True}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:out=list(ex.map(fetch,items))
(b/'metadata/CDH1_GENOTYPE_DOWNLOAD_COMPLETE_342.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out),flush=True)
