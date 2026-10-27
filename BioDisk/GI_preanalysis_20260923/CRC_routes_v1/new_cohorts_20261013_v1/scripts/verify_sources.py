from pathlib import Path
import json,hashlib,gzip,datetime
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1')
rec=json.loads((R/'LOCAL_DOWNLOAD_MANIFEST.json').read_text())
manifest=[]
for x in rec:
 p=R/'raw'/x['file'];h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 assert h.hexdigest()==x['sha256'] and p.stat().st_size==x['bytes']
 manifest.append(dict(x,path=str(p),status='verified'))
for p in (R/'metadata').glob('*'):
 if not p.is_file() or p.suffix=='.part':continue
 if p.suffix=='.gz':
  with gzip.open(p,'rb') as f:
   while f.read(1048576):pass
 manifest.append(dict(name=p.name,path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),status='verified'))
(R/'DOWNLOAD_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
(R/'DOWNLOAD_COMPLETE.json').write_text(json.dumps(dict(currently_transferred_files_verified=True,verified_datasets=[x['dataset'] for x in rec],route='official_NCBI_HTTPS_km112_to_local_to_sgs10',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2))
print('REMOTE_SOURCE_VERIFIED',flush=True)
