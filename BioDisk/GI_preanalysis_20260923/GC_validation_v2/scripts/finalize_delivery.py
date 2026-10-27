import pathlib,json,csv,zipfile,hashlib
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2');D=R/'05_reports'
old='修订针对细胞身份，在检查外部目标基因显著性之外完成；初始与修订结果分目录保留。';new='修订依据参考标签的一致性和谱系特异性，不以外部终点显著性作为阈值选择标准；初始外部结果已查看，初始与修订结果分目录保留。'
for p in [D/'internal_evidence_zh.md',R/'scripts'/'internal_report.py']:
 s=p.read_text(encoding='utf-8-sig');p.write_text(s.replace(old,new),encoding='utf-8')
(D/'client_message_zh.txt').write_text('我们已完成胃癌腹膜转移方向的进一步预验证，建议启动“脂质相关髓系微环境”的正式候选机制筛选。现已有17对患者组织、4对同患者巨噬细胞及独立28例腹腔细胞的数据基础。结果将重点指向腹膜脂质环境与髓系细胞状态，因此正式分析会先区分细胞比例与细胞状态，再筛选具有创新增量的可干预节点。后续以靶向脂质组、组织定位和关键功能实验支持主线。建议先开展有明确验收标准的生信发现阶段，候选通过重复性和查新后再集中投入湿实验。\n',encoding='utf-8')
checks=list(csv.DictReader((R/'logs'/'validation_checks.tsv').open(),delimiter='\t'));assert all(z['pass']=='True' for z in checks)
qa=json.loads((D/'build_artifacts'/'document_QA.json').read_text());assert qa['status']=='pass' and qa['pages']==6
for ext,key in [('docx','docx_sha256'),('pdf','pdf_sha256')]:
 p=D/('胃癌腹膜转移研究建议书.'+ext);assert hashlib.sha256(p.read_bytes()).hexdigest()==qa[key]
status={'status':'complete','date':'2026-09-23','direction1':'paused by user; no new CRC analysis','direction2':'prevalidation completed; client proposal supports gated formal discovery','scientific_checks_passed':len(checks),'client_pages':6,'all_pages_visually_verified':True,'patient_units':{'tissue_pairs':17,'paired_macrophage_patients':4,'external_localization_patients':28},'final_external_gate':'C1Q macrophage-enriched, endpoint genes excluded','no_confirmed_new_mechanism':True,'no_unfinished_analysis_jobs':True}
(R/'logs'/'TASK_COMPLETE.json').write_text(json.dumps(status,ensure_ascii=False,indent=2));(R/'00_docs'/'DECISION.md').write_text('# 第二轮决策\n\n方向1按用户要求暂缓。方向2已完成本轮验证与客户版交付，支持有验收门槛的正式候选机制筛选。\n\n保留组织组成解释、FABP4依赖、4对探索性样本和外部定位边界；最终候选和机制尚未确认。\n',encoding='utf-8')
# Refresh final input/script checksums after all text/code edits.
manifest=[]
for top in ['01_metadata','02_raw','scripts']:
 for f in sorted((R/top).glob('*')):
  if not f.is_file() or f.name=='integrity_manifest.tsv' or f.name.endswith('.part'):continue
  h=hashlib.sha256()
  with f.open('rb') as stream:
   while True:
    chunk=stream.read(8*1024*1024)
    if not chunk:break
    h.update(chunk)
  manifest.append({'path':str(f.relative_to(R)),'bytes':f.stat().st_size,'sha256':h.hexdigest()})
with (R/'01_metadata'/'integrity_manifest.tsv').open('w',newline='') as stream:
 writer=csv.DictWriter(stream,fieldnames=['path','bytes','sha256'],delimiter='\t');writer.writeheader();writer.writerows(manifest)

pack=D/'GC_validation_evidence_package.zip'
with zipfile.ZipFile(pack,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for top in ['00_docs','04_results','scripts','logs','05_reports']:
  for p in (R/top).rglob('*'):
   if p.is_file() and p!=pack and p.suffix!='.zip' and '__pycache__' not in str(p):z.write(p,str(p.relative_to(R)))
 for p in [R/'01_metadata'/'integrity_manifest.tsv',R/'01_metadata'/'novelty_review.json']:
  z.write(p,str(p.relative_to(R)))
with zipfile.ZipFile(pack) as z:assert z.testzip() is None
j={'package':pack.name,'bytes':pack.stat().st_size,'sha256':hashlib.sha256(pack.read_bytes()).hexdigest(),'status':status};(D/'delivery_manifest.json').write_text(json.dumps(j,ensure_ascii=False,indent=2));print(json.dumps(j,ensure_ascii=False))
