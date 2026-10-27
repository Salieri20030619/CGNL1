from pathlib import Path
import json,subprocess,concurrent.futures,hashlib,datetime,time,os
b=Path(__file__).resolve().parents[1];m=json.loads((b/'metadata/PILOT_DOWNLOAD_MANIFEST_440.json').read_text());dest=Path(m['destination']);root=dest.parent;stage=root/'aspera_probe_480';exe=json.loads((b/'metadata/ASPERA_CLIENT_475.json').read_text())['executable'];key=root/'software/aspera_474/ibm_public_client_rsa.pem'
def digest(p):
 h=hashlib.md5()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(8388608),b''):h.update(x)
 return h.hexdigest()
def run(item):
 target=dest/item['name'];received=stage/item['name'];logs=stage/('production_'+item['name']);logs.mkdir(exist_ok=True)
 if target.exists():assert target.stat().st_size==item['bytes'] and digest(target)==item['md5'];return dict(item,verified=True,local_path=str(target))
 for attempt in range(1,4):
  cmd=[exe,'-Q','-k','1','-l','20M','-P','33001','-i',str(key),'-L',str(logs),'era-fasp@fasp.sra.ebi.ac.uk:/vol1/run/ERR100/ERR10026294/'+item['name'],str(stage)]
  with (b/'logs'/('ASPERA484_'+item['name']+'.log')).open('ab') as f:r=subprocess.run(cmd,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT)
  if r.returncode==0:break
  print(item['name'],'attempt',attempt,'exit',r.returncode,flush=True)
  if attempt==3:raise RuntimeError('Aspera failed; original HTTPS partial preserved')
  time.sleep(5)
 assert received.stat().st_size==item['bytes'];actual=digest(received);assert actual==item['md5'],'Source MD5 mismatch; staged file preserved'
 assert not target.exists();received.rename(target)
 receipt=dict(item,verified=True,local_path=str(target),transfer='server-direct Aspera',attempts=attempt,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(b/'metadata'/('ASPERA484_'+item['name']+'.json')).write_text(json.dumps(receipt,indent=2),encoding='utf-8');print('VERIFIED',item['name'],flush=True);return receipt
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:files=list(pool.map(run,m['items']))
(b/'metadata/PILOT_DOWNLOAD_COMPLETE_484.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':files,'run':m['run'],'all_md5_verified':True,'server_direct':True},indent=2),encoding='utf-8')
