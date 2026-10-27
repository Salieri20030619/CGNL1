from pathlib import Path
import subprocess,requests,concurrent.futures,json,hashlib,gzip,time,datetime,re,os
L=Path(__file__).resolve().parent;D=L/'raw';D.mkdir(exist_ok=True)
curl=Path(os.environ['SystemRoot'])/'System32/curl.exe'
records=[]
for ds,ext in [('GSE200997','csv'),('GSE132465','txt')]:
 name=ds+'_GEO_processed_CRC_10X_raw_UMI_count_matrix.'+ext+'.gz'
 url='http://ftp.ncbi.nlm.nih.gov/geo/series/'+ds[:-3]+'nnn/'+ds+'/suppl/'+name
 s=requests.Session();s.trust_env=False
 with s.get(url,headers={'Range':'bytes=0-0'},stream=True,timeout=(12,20)) as r:
  r.raise_for_status();cr=r.headers['Content-Range'];assert cr.startswith('bytes 0-0/'),cr;total=int(cr.split('/')[-1]);assert len(r.content)==1
 parts=D/(ds+'_parts');parts.mkdir(exist_ok=True)
 ranges=[(i,start,min(start+4*1024**2,total)-1) for i,start in enumerate(range(0,total,4*1024**2))]
 def fetch(item):
  i,start,end=item;p=parts/(str(i).zfill(4)+'.part');h=parts/(str(i).zfill(4)+'.headers');log=parts/(str(i).zfill(4)+'.log')
  for attempt in range(2):
   with log.open('ab') as err:
    proc=subprocess.run([str(curl),'--noproxy','*','--fail','--location','--silent','--show-error','--connect-timeout','12','--max-time','180','--speed-limit','1024','--speed-time','30','--range',str(start)+'-'+str(end),'--dump-header',str(h),'--output',str(p),url],stderr=err,creationflags=subprocess.CREATE_NO_WINDOW)
   if proc.returncode==0 and p.stat().st_size==end-start+1 and ('content-range: bytes '+str(start)+'-'+str(end)+'/'+str(total)) in h.read_text().lower():
    print(ds,'PART',i+1,'/',len(ranges),'OK',flush=True);return p
  raise RuntimeError('download failed '+str(p)+' '+log.read_text())
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:files=list(pool.map(fetch,ranges))
 dest=D/name;sha=hashlib.sha256()
 with dest.open('wb') as out:
  for p in files:
   b=p.read_bytes();out.write(b);sha.update(b)
 assert dest.stat().st_size==total
 with gzip.open(dest,'rb') as f:
  while f.read(1048576):pass
 rec=dict(dataset=ds,file=name,url=url,bytes=total,sha256=sha.hexdigest(),status='verified',range_parts=len(ranges),local=str(dest))
 records.append(rec);(L/'LOCAL_DOWNLOAD_MANIFEST.json').write_text(json.dumps(records,indent=2))
 print('FILE_VERIFIED',json.dumps(rec),flush=True)
for rec in records:
 subprocess.run(['scp','-q',rec['local'],'sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1/raw/'+rec['file']],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
subprocess.run(['scp','-q',str(L/'LOCAL_DOWNLOAD_MANIFEST.json'),'sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1/LOCAL_DOWNLOAD_MANIFEST.json'],check=True,creationflags=subprocess.CREATE_NO_WINDOW)
print('DOWNLOAD_AND_UPLOAD_COMPLETE',flush=True)
