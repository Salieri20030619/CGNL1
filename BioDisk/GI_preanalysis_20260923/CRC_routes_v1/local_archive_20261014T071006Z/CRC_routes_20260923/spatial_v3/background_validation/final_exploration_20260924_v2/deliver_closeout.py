from pathlib import Path
import zipfile,hashlib,json,shutil
Q=Path(__file__).parent
D=Path(r'C:\Users\31937\Downloads\结肠癌肝转移_最后一轮探索_20260924');D.mkdir(exist_ok=True)
paths=[p for p in Q.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix.lower() in ['.py','.tsv','.json','.md','.png','.log'] and p.name not in ['deliver_closeout.py','stop_superseded.py']]
manifest=[{'path':p.relative_to(Q).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]
(Q/'CLOSEOUT/DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
for p in (Q/'CLOSEOUT').iterdir():
 if p.is_file():shutil.copy2(p,D/p.name)
Z=D/'最后一轮探索_可复核证据包.zip'
with zipfile.ZipFile(Z,'w',zipfile.ZIP_DEFLATED) as z:
 for p in paths:z.write(p,p.relative_to(Q))
 z.write(Q/'CLOSEOUT/DELIVERY_MANIFEST.json','CLOSEOUT/DELIVERY_MANIFEST.json')
with zipfile.ZipFile(Z) as z:assert z.testzip() is None
result={'zip_bytes':Z.stat().st_size,'zip_sha256':hashlib.sha256(Z.read_bytes()).hexdigest(),'files':len(paths),'status':'final_exploration_delivered','next_phase':'formal_study_plan','jobs_running':0}
(D/'DELIVERY_STATUS.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result))
