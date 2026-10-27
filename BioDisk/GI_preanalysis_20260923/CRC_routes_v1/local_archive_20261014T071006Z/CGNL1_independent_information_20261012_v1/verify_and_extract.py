from pathlib import Path
import zipfile,hashlib,json
from PIL import Image,ImageOps,ImageDraw
r=Path(r'C:\Users\31937\CGNL1_independent_information_20261012_v1')
p=r/'CGNL1_independent_information_bundle.zip';digest=hashlib.sha256(p.read_bytes()).hexdigest()
assert digest=='98ca7e29ce91e5e5978c37f01624a57d96fb90f23f7f16fa6300eafb253e6d28',digest
dest=r/'delivery';dest.mkdir(exist_ok=True);root=dest.resolve()
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None
 n=0
 for line in z.read('SHA256SUMS.txt').decode('utf-8').splitlines():
  h,name=line.split(maxsplit=1);assert hashlib.sha256(z.read(name)).hexdigest()==h,name;n+=1
 for item in z.infolist():
  out=(dest/item.filename).resolve();assert out==root or root in out.parents
  if item.is_dir():out.mkdir(parents=True,exist_ok=True)
  else:out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(z.read(item.filename))
report=(dest/'CGNL1_independent_information_report.html').read_text(encoding='utf-8')
assert '<html lang="zh-CN">' in report and len(report)>10000
files=['A_conditional_associations.png','A_increment_marker_comparison.png','B_paired_specificity.png','C_spatial_attenuation.png']
sheet=Image.new('RGB',(1800,1300),'white')
for i,name in enumerate(files):
 im=Image.open(dest/'figures'/name).convert('RGB');im.thumbnail((870,620))
 x=(i%2)*900+(900-im.width)//2;y=(i//2)*650+(650-im.height)//2;sheet.paste(im,(x,y))
sheet.save(r/'visual_qa.jpg',quality=94)
check={'zip_sha256':digest,'manifest_files_checked':n,'failed':0,'safe_extraction':True,'report_utf8_ok':True}
(r/'LOCAL_DELIVERY_CHECK.json').write_text(json.dumps(check,indent=2),encoding='utf-8')
print(json.dumps(check))
