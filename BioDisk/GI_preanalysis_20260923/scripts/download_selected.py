import pathlib,subprocess,hashlib,json,concurrent.futures,os,time,gzip,xml.etree.ElementTree as ET
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
def run(item):
 name,url=item; p=R/'02_raw'/name; log=R/'logs'/(name+'.fetch.log'); t=time.time()
 rec={'name':name,'url':url,'relative_path':str(p.relative_to(R)),'start':time.strftime('%FT%T%z')}
 try:
  with log.open('w') as f:
   rc=subprocess.run(['curl','-L','--fail','--connect-timeout','20','--max-time','900','--retry','2','--retry-delay','2','--speed-limit','1024','--speed-time','90','-o',str(p)+'.part',url],stdout=f,stderr=f).returncode
  rec['exit_code']=rc
  if rc:raise ValueError('download failed')
  tmp=pathlib.Path(str(p)+'.part')
  if tmp.stat().st_size==0:raise ValueError('empty response')
  if name.endswith('.gz'):
   with gzip.open(tmp,'rb') as f:
    while f.read(8*1024*1024):pass
  if name.endswith('.h5'):
   with tmp.open('rb') as f:assert f.read(8)==b'\x89HDF\r\n\x1a\n'
  with tmp.open('rb') as f:
   h=hashlib.sha256()
   while True:
    b=f.read(8*1024*1024)
    if not b:break
    h.update(b)
  rec.update(bytes=tmp.stat().st_size,sha256=h.hexdigest(),status='verified_download')
  os.replace(tmp,p)
 except Exception as e:rec.update(status='failed',error=str(e))
 rec['seconds']=round(time.time()-t,2); (R/'logs'/(name+'.fetch.json')).write_text(json.dumps(rec,indent=2));print(json.dumps(rec),flush=True);return rec
items=[]
for suffix in ['barcodes.tsv.gz','features.tsv.gz','meta.tsv.gz','count_matrix.mtx.gz']:
 name='GSE234129_'+suffix;items.append((name,'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234129/suppl/'+name))
items += [('GSE145356_counts.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145356/suppl/GSE145356_merged_HTSeq_readCounts_V2.txt.gz'),('GSE200697_family.soft.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200697/soft/GSE200697_family.soft.gz'),('PXD064711_protein.raw.intensity.txt','https://download.iprox.org/IPX0012159000/IPX0012159002/protein.raw.intensity.txt')]
for s in ['TP001','TO001','TP002','TO002']:
 name='GSE294559_'+s+'.h5';items.append((name,'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294559/suppl/'+name))
(R/'01_metadata'/'selected_downloads.json').write_text(json.dumps(items,indent=2))
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:res=list(ex.map(run,items))
(R/'logs'/'selected_downloads.complete.json').write_text(json.dumps({'exit_code':0 if all(x['status']=='verified_download' for x in res) else 1,'files':res},indent=2))
