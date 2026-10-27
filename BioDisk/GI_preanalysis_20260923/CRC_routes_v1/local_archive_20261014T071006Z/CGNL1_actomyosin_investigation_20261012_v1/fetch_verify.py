from pathlib import Path
import subprocess,time,hashlib,zipfile,json
from PIL import Image
r=Path(r'C:\Users\31937\CGNL1_actomyosin_investigation_20261012_v1')
p=r/'CGNL1_actomyosin_investigation_bundle.zip'
remote='sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/actomyosin_investigation_20261012_v1/CGNL1_actomyosin_investigation_bundle.zip'
attempts=[]
for i in range(3):
 proc=subprocess.run(['scp','-q','-o','BatchMode=yes','-o','ConnectTimeout=20',remote,str(p)],capture_output=True,text=True,creationflags=subprocess.CREATE_NO_WINDOW)
 attempts.append({'attempt':i+1,'exit_code':proc.returncode,'stderr':proc.stderr})
 if proc.returncode==0:break
 time.sleep(2)
(r/'transfer_log.json').write_text(json.dumps(attempts,indent=2),encoding='utf-8')
assert proc.returncode==0,proc.stderr
digest=hashlib.sha256(p.read_bytes()).hexdigest()
assert digest=='d61cb3b8ef440229c03642cf54ce4b69a21d0d8852aa73c45ecfe3e80d1431b7',digest
dest=(r/'delivery').resolve();dest.mkdir(exist_ok=True)
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None
 checked=0
 for line in z.read('SHA256SUMS.txt').decode('utf-8').splitlines():
  h,n=line.split(maxsplit=1);assert hashlib.sha256(z.read(n)).hexdigest()==h,n;checked+=1
 for info in z.infolist():
  out=(dest/info.filename).resolve();assert out==dest or dest in out.parents
  if info.is_dir():out.mkdir(parents=True,exist_ok=True)
  else:out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(z.read(info.filename))
sheet=Image.new('RGB',(1800,1300),'white')
for i,name in enumerate(['components.png','robustness.png','external_replication.png','knockout.png']):
 im=Image.open(dest/'figures'/name).convert('RGB');im.thumbnail((870,620));sheet.paste(im,((i%2)*900+(900-im.width)//2,(i//2)*650+(650-im.height)//2))
sheet.save(r/'visual_qa.jpg',quality=94)
result={'sha256':digest,'files_verified':checked,'failed':0,'report_utf8_ok':'CGNL1' in (dest/'CGNL1_actomyosin_investigation_report.html').read_text(encoding='utf-8'),'transfer_attempts':len(attempts)}
(r/'LOCAL_DELIVERY_CHECK.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))
