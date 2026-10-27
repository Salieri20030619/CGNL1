import pathlib,subprocess,json,time,hashlib
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');D=O/'raw';D.mkdir(exist_ok=True);u='https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7058nnn/GSM7058755/suppl/GSM7058755_non_immune_counts.txt.gz';p=D/'nonimmune_counts.txt.gz';t=time.time()
with (D/'download_counts.log').open('w') as f:r=subprocess.run(['curl.exe','-x','http://127.0.0.1:7897','-L','--fail','--max-time','1800','--retry','2','-C','-','-o',str(p)+'.part',u],stdout=f,stderr=f)
j={'url':u,'exit_code':r.returncode,'elapsed_seconds':time.time()-t}
if r.returncode==0:
 pathlib.Path(str(p)+'.part').replace(p);j['bytes']=p.stat().st_size;j['sha256']=hashlib.file_digest(p.open('rb'),'sha256').hexdigest();j['transfer_exit_code']=subprocess.run(['scp',str(p),'sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/raw/']).returncode
(D/'COUNTS_DOWNLOAD.json').write_text(json.dumps(j,indent=2));print(json.dumps(j),flush=True)
