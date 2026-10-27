from pathlib import Path
import subprocess,json,hashlib,sys,base64
L=Path(__file__).resolve().parent;D=L/'raw';D.mkdir(exist_ok=True)
K='/data/nas1/wangmingyuan_OD/CGNL1_public_download_20261013_v1'
R='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1'
def run(a,**kw):return subprocess.run(a,check=True,creationflags=subprocess.CREATE_NO_WINDOW,**kw)
records=json.loads((L/'LOCAL_DOWNLOAD_MANIFEST.json').read_text()) if (L/'LOCAL_DOWNLOAD_MANIFEST.json').exists() else []
for ds in sys.argv[1:]:
 source='SEGMENT132_COMPLETE.json' if ds=='GSE132465' else 'MANIFEST.json'
 response=run(['ssh','-T','-o','BatchMode=yes','km112','cat '+K+'/'+source],capture_output=True,text=True)
 meta=json.loads(response.stdout);rec=meta if isinstance(meta,dict) else next(x for x in meta if x['dataset']==ds)
 assert rec['status']=='verified';dest=D/rec['file']
 run(['scp','-q','km112:'+K+'/'+rec['file'],str(dest)])
 assert dest.stat().st_size==rec['bytes'] and hashlib.sha256(dest.read_bytes()).hexdigest()==rec['sha256']
 run(['scp','-q',str(dest),'sgs10:'+R+'/raw/'+rec['file']])
 records=[r for r in records if r['dataset']!=ds]+[dict(rec,local=str(dest),transport='km112_to_local_to_sgs10')]
(L/'LOCAL_DOWNLOAD_MANIFEST.json').write_text(json.dumps(records,indent=2))
run(['scp','-q',str(L/'LOCAL_DOWNLOAD_MANIFEST.json'),'sgs10:'+R+'/LOCAL_DOWNLOAD_MANIFEST.json'])
print('CROSS_SERVER_TRANSFER_VERIFIED',[x['dataset'] for x in records],flush=True)
