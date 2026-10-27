from pathlib import Path
import tarfile,json,subprocess,re,datetime
O=Path.cwd();D=O/'downloads';out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'archive_entries':[]}
try:
 with tarfile.open(D/'HGP_archive_range.response','r|*') as tf:
  for m in tf:
   out['archive_entries'].append({'name':m.name,'size':m.size,'offset_data':m.offset_data})
   if len(out['archive_entries'])>=8:break
except Exception as e:out['partial_archive_status']=str(e)
u='https://raw.githubusercontent.com/TissUUmaps/TissUUmaps/master/tissuumaps/views.py';p=D/'TissUUmaps_views.py';r=subprocess.run(['curl','--noproxy','*','-L','--connect-timeout','10','--max-time','30','-sS','-w','%{http_code}','-o',str(p),u],capture_output=True,timeout=35);out['route_source']={'url':u,'exit':r.returncode,'http':r.stdout.decode('utf-8')}
if r.returncode==0 and r.stdout==b'200':
 lines=p.read_text(encoding='utf-8').splitlines();out['h5_routes']=[]
 for i,line in enumerate(lines):
  if 'h5ad' in line or 'HTTP_RANGE' in line or 'X-Requested-With' in line:out['h5_routes'].append('\n'.join(lines[max(0,i-5):i+12]))
(O/'metadata/ARCHIVE_ROUTE_CHECK.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False))
