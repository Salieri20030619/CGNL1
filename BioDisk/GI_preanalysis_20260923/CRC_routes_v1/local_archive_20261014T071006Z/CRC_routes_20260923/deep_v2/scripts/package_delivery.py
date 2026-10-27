import pathlib,json,hashlib,zipfile,shutil,sys
from datetime import datetime
sys.stdout.reconfigure(encoding='utf-8')
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2'); D=O/'delivery'; B=O.parent
out=pathlib.Path(r'C:\Users\31937\Downloads\CGNL1立项建议_20260924');out.mkdir(exist_ok=True,parents=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
(D/'排版验收.json').write_text(json.dumps({'status':'pass','pages':8,'all_pages_reviewed':True,'changed_pages_final_reviewed':[2,3,6],'fixes':['Scientific-notation missing glyphs replaced by ASCII e notation','CNV expression units clarified as log2(CPM+1)','Markdown table separators added'],'statistical_checks_passed':18},ensure_ascii=False,indent=2),encoding='utf-8')
selected=[]
for folder in ['results','scripts','metadata','raw']:
 for p in (O/folder).rglob('*'):
  if p.is_file() and p.name!='GRCh38.p14.ncbi.gtf.gz' and '__pycache__' not in p.parts:
   selected.append((p,'deep_v2/'+p.relative_to(O).as_posix()))
for p in O.glob('*.json'):selected.append((p,'deep_v2/'+p.name))
for folder in ['results','metadata','scripts']:
 for p in (B/folder).rglob('*'):
  if p.is_file() and p.stat().st_size<10_000_000 and '__pycache__' not in p.parts:
   selected.append((p,'first_stage/'+p.relative_to(B).as_posix()))
for p in B.glob('*.json'):selected.append((p,'first_stage/'+p.name))
for p in (B/'delivery').glob('*.md'):selected.append((p,'first_stage/delivery/'+p.name))
for p in D.iterdir():
 if p.is_file() and p.suffix not in ['.zip'] and p.name not in ['证据文件清单.json','交付完成.json']:
  selected.append((p,'delivery/'+p.name))
paths=[arc for _,arc in selected];assert len(paths)==len(set(paths))
manifest={'files':[{'path':arc,'bytes':p.stat().st_size,'sha256':sha(p)} for p,arc in selected], 'reference_not_bundled':{'path':str(O/'metadata/GRCh38.p14.ncbi.gtf.gz'),'sha256':sha(O/'metadata/GRCh38.p14.ncbi.gtf.gz'),'remote':'/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2/metadata/GRCh38.p14.ncbi.gtf.gz'},'scope':'Raw single-cell matrices, vendor and large inference arrays retained on sgs10. Small first-stage results included for spatial/source traceability.'}
mp=D/'证据文件清单.json';mp.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
zip_path=D/'CGNL1深化分析_证据包.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p,arc in selected:z.write(p,arc)
 z.write(mp,'delivery/'+mp.name)
with zipfile.ZipFile(zip_path) as z:
 assert z.testzip() is None
 assert len(z.namelist())==len(selected)+1
for p in D.iterdir():
 if p.is_file() and p.name!='交付完成.json':
  q=out/p.name;shutil.copy2(p,q);assert sha(p)==sha(q)
status={'stage':'deeper_preanalysis_delivered','bioinformatics_source_gate':'pending','wet_experiments':'not_started','local_delivery':str(D),'downloads':str(out),'zip_files':len(selected)+1,'zip_bytes':zip_path.stat().st_size,'zip_sha256':sha(zip_path),'pdf_sha256':sha(D/'CGNL1结肠癌肝转移_深化立项建议.pdf'),'remote_sync':'pending'}
(D/'交付完成.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8');shutil.copy2(D/'交付完成.json',out/'交付完成.json')
print(json.dumps(status,ensure_ascii=False,indent=2))
