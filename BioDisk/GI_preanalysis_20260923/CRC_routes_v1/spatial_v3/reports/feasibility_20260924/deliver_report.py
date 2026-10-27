from pathlib import Path
import json,shutil,hashlib
O=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3\reports\feasibility_20260924');D=Path(r'C:\Users\31937\Downloads\CGNL1项目可行性报告_20260924');D.mkdir(exist_ok=True)
status=json.loads((O/'REPORT_STATUS.json').read_text());status['render_review']='all8 pages passed';(O/'REPORT_STATUS.json').write_text(json.dumps(status,indent=2))
for p in O.iterdir():
 if p.is_file() and p.suffix in ['.pdf','.md','.png','.json']:
  shutil.copy2(p,D/p.name)
  assert hashlib.sha256(p.read_bytes()).hexdigest()==hashlib.sha256((D/p.name).read_bytes()).hexdigest()
print('Feasibility report copied to Downloads')
