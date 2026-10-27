import pathlib,json,hashlib,zipfile,shutil,sys
sys.stdout.reconfigure(encoding='utf-8')
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3');D=O/'delivery';D.mkdir(exist_ok=True);out=pathlib.Path(r'C:\Users\31937\Downloads\CGNL1空间初筛_20260924');out.mkdir(exist_ok=True)
report=O/'空间初筛结果.md';shutil.copy2(O/'background_validation/后台验证计划.md',D/'后台验证计划.md')
qa={'initial_screen':'complete','sections':24,'mapped_patients':11,'QC_spots':219762,'primary_independent_patients':11,'primary_raw_count_reconstruction':'pass','FDR_recheck':'7 families passed','visual_review':'all24 maps reviewed; shared title clipping repaired; final summary and representative full-resolution maps checked','background_validation':'planned_not_launched','remaining':'cell-source identification, deconvolution, external spatial replication and perturbation response'}
(D/'完成状态.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'RUNNING_STATE.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
readme='''先阅读空间初筛结果.md，再阅读后台验证计划.md。后台验证当前仅完成规划，尚未启动。

原始空间包24份、合计5,777,998,995字节，保留sgs10：/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/spatial_v3/raw/liver。证据包不重复打包5.78GB原始文件，包含原始文件URL/SHA256、逐spot候选计数与坐标、全部患者统计、图件、方案版本和脚本。

完整复现入口（sgs10既有项目）：
/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python /mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/spatial_v3/scripts/screen_liver.py /mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/spatial_v3

运行脚本依赖原项目父目录的作者注释（results/spatial_author_metadata.tsv.gz，或metadata下同名文件）。证据包provenance保留该注释及患者映射副本。不要将旧screen_pilot结果当作正式结果；正式结果均在results/screen。初筛报告中的相对图链接在此目录中可用。
'''
(D/'交付说明.txt').write_text(readme,encoding='utf-8-sig')
payload=[]
for p in (O/'results/screen').rglob('*'):
 if p.is_file():payload.append((p,p.relative_to(O).as_posix()))
for name in ['screen_liver.py','summarize_screen.py','finalize_screen.py','redraw_spatial_maps.py','download_liver.py','run_on_input_completion.py','coordinate_calibration.py','probe_liver_features.py','plan_background_validation.py']:
 payload.append((O/'scripts'/name,'scripts/'+name))
for p in O.glob('SCREEN_PLAN*.json'):payload.append((p,p.name))
for name in ['LIVER_DOWNLOAD_PLAN.json','LIVER_FILE_SIZES.json']:
 payload.append((O/'metadata'/name,'metadata/'+name))
for name in ['LIVER_INPUTS_VERIFIED.json','SCREEN_RUN_EXIT.json','full_screen.log','coordinate_calibration_audit.json','IMPLEMENTATION_AMENDMENT.json']:
 payload.append((O/'results'/name,'results/'+name))
for name in ['spatial_author_metadata.tsv.gz','spatial_patient_map.tsv','spatial_table2.json']:
 payload.append((O.parent/'metadata'/name,'provenance/'+name))
payload.extend([(report,report.name),(O/'background_validation/tasks.json','background_validation/tasks.json')])
for p in D.iterdir():
 if p.is_file() and p.suffix!='.zip' and p.name!='文件校验清单.json':payload.append((p,p.name))
assert len(set(x[1] for x in payload))==len(payload)
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
manifest=[{'path':arc,'bytes':p.stat().st_size,'sha256':sha(p)} for p,arc in payload];mp=D/'文件校验清单.json';mp.write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8');payload.append((mp,mp.name))
zp=D/'CGNL1空间初筛_证据包.zip'
with zipfile.ZipFile(zp,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p,arc in payload:z.write(p,arc)
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
for p,arc in payload:
 dest=out/arc;dest.parent.mkdir(exist_ok=True,parents=True);shutil.copy2(p,dest);assert sha(p)==sha(dest)
shutil.copy2(zp,out/zp.name);assert sha(zp)==sha(out/zp.name)
status={'downloads':str(out),'zip_bytes':zp.stat().st_size,'zip_sha256':sha(zp),'report_sha256':sha(report),'background':'planned_not_launched','files':len(payload)};(D/'交付校验.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8');shutil.copy2(D/'交付校验.json',out/'交付校验.json')
print(json.dumps(status,ensure_ascii=False,indent=2))
