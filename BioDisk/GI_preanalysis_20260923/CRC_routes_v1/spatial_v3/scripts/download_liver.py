import pathlib,json,subprocess,concurrent.futures,hashlib,time,tarfile,sys
O=pathlib.Path(sys.argv[1]);D=O/'raw/liver';D.mkdir(exist_ok=True,parents=True);jobs=json.loads((O/'metadata/LIVER_DOWNLOAD_PLAN.json').read_text())
if len(sys.argv)>3:jobs=[j for j in jobs if j['sample'] in sys.argv[3].split(',')]
def fetch(j):
 p=D/(j['sample']+'.tar.gz');start=time.time();record=dict(j)
 try:
  if not p.exists():
   with (D/(j['sample']+'.download.log')).open('w') as f:r=subprocess.run(['curl']+(['-x',sys.argv[2]] if len(sys.argv)>2 else [])+['-L','--fail','--retry','3','--retry-all-errors','--retry-delay','2','--max-time','1200','-C','-','-o',str(p)+'.part',j['url']],stdout=f,stderr=f)
   if r.returncode:raise RuntimeError('curl exit '+str(r.returncode))
   pathlib.Path(str(p)+'.part').replace(p)
  with tarfile.open(p) as t:
   names=t.getnames();assert any(n.endswith('matrix.mtx.gz') for n in names);assert any('tissue_positions' in n for n in names)
  with p.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
  record.update(ok=True,bytes=p.stat().st_size,sha256=h,members=len(names))
 except Exception as e:record.update(ok=False,error=str(e))
 record['seconds']=time.time()-start;(D/(j['sample']+'.complete.json')).write_text(json.dumps(record,indent=2));print(json.dumps(record),flush=True);return record
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:r=list(ex.map(fetch,jobs))
(O/'results/LIVER_DOWNLOAD_COMPLETE.json').write_text(json.dumps(r,indent=2));assert all(j['ok'] for j in r),'At least one input failed';print('VERIFIED '+str(len(jobs))+' DOWNLOADS',flush=True)
