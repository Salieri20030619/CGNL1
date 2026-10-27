import pathlib,json,hashlib,zipfile,shutil,sys
sys.stdout.reconfigure(encoding='utf-8');O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');D=O/'delivery';m=json.loads((O/'FINAL_INPUT_MANIFEST.json').read_text());q=json.loads((O/'results/SERVER_INPUT_INTEGRITY.json').read_text());assert q['exit_code']==0 and len(q['files'])==9 and all(x['matches_local'] for x in q['files']);assert len(json.loads((O/'QA.json').read_text(encoding='utf-8'))['checks'])==183
sources=['clinical_eligibility.json','spatial_patient_map.tsv','spatial_table2.xlsx','spatial_table2.json','nonimmune_meta.txt.gz','visium_meta_after_qc.rds','zenodo19043057.json','zenodo_latest.json','visium_Readme.md','GSE225857.soft.gz','GSE294385.soft.gz']
zpath=D/'结肠癌两方向预分析_完整证据包.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
 for p in [D/'结肠癌两方向首轮预分析报告.md',D/'README.txt']:z.write(p,p.name)
 for root,prefix in [(O/'results','results'),(O/'scripts','scripts')]:
  for p in root.rglob('*'):
   if p.is_file():z.write(p,prefix+'/'+p.relative_to(root).as_posix())
 for p in (D/'figures').iterdir():z.write(p,'figures/'+p.name)
 for name in ['PLAN.json','ANALYSIS_PLAN.json','FINAL_INPUT_MANIFEST.json','QA.json','metadata_downloads.json','spatial_download_plan.json']:z.write(O/name,name)
 for name in sources:z.write(O/'metadata'/name,'metadata/'+name)
 # Preserve both failed attempts and eventual successes for download auditing.
 for p in (O/'raw').glob('*.json'):z.write(p,'download_records/'+p.name)
with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
out=pathlib.Path(r'C:\Users\31937\Downloads')/'结肠癌两方向预分析_20260924';out.mkdir(exist_ok=True);shutil.copy2(D/'结肠癌两方向首轮预分析报告.md',out);shutil.copytree(D/'figures',out/'figures',dirs_exist_ok=True);target=out.parent/zpath.name;shutil.copy2(zpath,target)
h=lambda p:hashlib.file_digest(p.open('rb'),'sha256').hexdigest();assert h(target)==h(zpath);assert h(out/'结肠癌两方向首轮预分析报告.md')==h(D/'结肠癌两方向首轮预分析报告.md');manifest={p.name:{'bytes':p.stat().st_size,'sha256':h(p)} for p in [D/'结肠癌两方向首轮预分析报告.md',zpath]};(D/'SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({'report':str(out/'结肠癌两方向首轮预分析报告.md'),'package':str(target),'manifest':manifest,'archive_crc':'passed'},ensure_ascii=False,indent=2))
