from pathlib import Path
import requests,re,concurrent.futures,json,hashlib,subprocess,os
b=Path(r'C:\Users\31937\CGNL1_v3_pilot_20261006\direct_local');tasks=[]
for rel in re.findall(r'href="([^"]*cp311[^"\s]*linux_x86_64.whl)"',(b/'pyg_cpu.html').read_text(encoding='utf-8')):
 tasks.append((rel.split('/')[-1].replace('%2B','+'),'https://data.pyg.org/whl/'+rel,b/'wheels'))
tasks += [('GSM9815454_cells.parquet.gz','https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9815nnn/GSM9815454/suppl/GSM9815454_tumor1_M242T_cells.parquet.gz',b)]
def f(it):
 n,u,d=it;p=d/n;part=p.with_suffix(p.suffix+'.part')
 env=os.environ.copy()
 for k in list(env):
  if 'proxy' in k.lower():env.pop(k,None)
 env.update(NO_PROXY='*',no_proxy='*')
 cmd=[r'C:\Windows\System32\curl.exe','-q','-4','--proxy','','--noproxy','*','-L','--fail','--connect-timeout','12','--max-time','180','--retry','1','--retry-delay','1','-o',str(part),u]
 r=subprocess.run(cmd,env=env,creationflags=subprocess.CREATE_NO_WINDOW,capture_output=True,text=True,timeout=390)
 rec={'file':n,'url':u,'exit_code':r.returncode,'network_mode':'direct_only'}
 if r.returncode==0:
  part.replace(p);rec.update(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 else:rec['error']=r.stderr[-300:]
 return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rs=list(ex.map(f,tasks))
(b/'additional_download_log.json').write_text(json.dumps(rs,indent=2),encoding='utf-8');print(json.dumps(rs))
