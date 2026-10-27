from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,zipfile,html,base64,io
from bs4 import BeautifulSoup
from PIL import Image
import pandas as pd
O=Path.cwd();R=O/'report';now=datetime.now(timezone.utc).isoformat()
source=R/'CGNL1_spatial_mainline_review.html';b=source.read_bytes();text=b.decode('utf-8',errors='strict');assert not b.startswith(b'\xef\xbb\xbf');assert '\ufffd' not in text
soup=BeautifulSoup(text,'html.parser');images=[]
for tag in soup.find_all('img'):
 src=tag['src'];assert src.startswith('data:image/png;base64,');im=Image.open(io.BytesIO(base64.b64decode(src.split(',',1)[1])));im.verify();images.append({'alt':tag.get('alt'),'valid_png':True})
assert len(images)==3
render=json.loads((O/'metadata/LOCAL_RENDER_CHECK.json').read_text());assert render['exit_code']==0 and render['screenshot_exists']
assert 'CC' not in soup.title.get_text()
tables=pd.read_html(io.StringIO(text));assert len(tables)>=10
manifest=json.loads((O/'baseline/SOURCE_MANIFEST.json').read_text());assert all(hashlib.sha256(Path(x['source']).read_bytes()).hexdigest()==x['sha256'] for x in manifest)
searches=['"GSE294385" spatial colorectal liver metastasis','"CGNL1" "spatial" metastasis','colorectal liver metastasis tumor liver interface spatial transcriptomics replacement desmoplastic 2025 2026','"Spatial Multi-Omics Landscape of Colorectal Cancer Macro- and Micrometastases"','"CGNL1" "interface" cancer','"CGNL1" "spatial transcriptomics"','"Evidence of endothelial anergy" "GSE"','"CGNL1" "colorectal" "spatial"','"paracingulin" "colorectal" metastasis','"GSE206552" liver spatial','"GSE254364" liver spatial','"Prognostic value and spatial transcriptomic landscape" histopathological growth patterns 2026','"Spatial Multi-Omics Landscape" "CGNL1"','"liver metastases" "CGNL1" -site:genecards.org -site:editxor.com','"37157887" spatial senescence colorectal liver four','"36652202" spatial Visium two patients','"TissUUmaps" "h5ad" "h5" send_file']
(O/'literature/WEB_SEARCH_AUDIT.json').write_text(json.dumps({'recorded_utc':now,'execution_context':'assistant web search during current authorized turn; exact query strings recorded; repository/API downloads have separate receipts','queries':searches,'scope':'focused scientific novelty, available material and public viewer route qualification; not systematic literature review'},ensure_ascii=False,indent=2),encoding='utf-8')
qa={'utc':now,'UTF8_no_BOM':True,'embedded_figures_valid':images,'HTML_tables':len(tables),'render_exit_zero':True,'visual_check':'Chinese title, opening decision and table rendered in headless Edge; first two scientific figures inspected separately; all-section map included but not manual histopathology validation','original_snapshots_unchanged':True,'primary_numeric_qa':json.loads((O/'NUMERICAL_QA.json').read_text())['all_pass'],'report_status':'ready_for_user_review'}
(O/'REPORT_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
readme=f'''CGNL1 空间主线核定：审核包

先打开 report/CGNL1_spatial_mainline_review.html；图件已嵌入，可离线阅读。
简要结论见 report/CGNL1_空间主线核定_审核摘要.txt。
本轮保留一个候选问题，但未核定为已获充分数据支持的正式主线。

results/：完整主/敏感性汇总、逐患者效应、资格及诊断表（不包含原始表达矩阵）。
figures/：最终图件。
方案、决策、QA和文献登记保存于同名路径。
全量逐spot测量、权重、原始响应、脚本、运行日志、旧结果快照及下载未完成记录留在服务器：
{O}

结果版本：spatial_mainline_validation_20261016_v1。
UTC创建于2026-10-16，对应纽约当地2026-10-15。
原结果未覆盖；研究项目整体未标记完成；当前报告等待用户审核。
'''
(O/'REVIEW_PACKAGE_README.txt').write_text(readme,encoding='utf-8')
exclude_results={'spatial_spot_features.tsv.gz','overlap_weights.tsv.gz'}
files=[p for p in R.iterdir() if p.is_file()]+[p for p in (O/'figures').iterdir() if p.is_file()]
files += [p for p in (O/'results').iterdir() if p.is_file() and p.name not in exclude_results and not p.name.startswith('matched_pairs_')]
files += [O/n for n in ['BASELINE_AND_AUTHORIZATION.json','ANALYSIS_SPECIFICATION.json','SPEC_FREEZE.json','IMPLEMENTATION_FREEZE.json','POST_PRIMARY_DIAGNOSTIC_SPEC.json','OVERLAP_VALIDATION_SPEC.json','GEOMETRY_DIAGNOSTIC_SPEC.json','GEODESIC_VALIDATION_SPEC.json','MAINLINE_DECISION.json','NUMERICAL_QA.json','FINAL_SUMMARY_QA.json','REPORT_QA.json','REVIEW_PACKAGE_README.txt']]
files += [O/'literature'/n for n in ['LITERATURE_REGISTER.json','WEB_SEARCH_AUDIT.json','FETCH_RECEIPTS.json','FETCH_REPAIR_RECEIPTS.json','PRIMARY_SOURCE_INSPECTION.json']]
files += [O/'metadata'/n for n in ['RAW_SOURCE_MANIFEST.json','interface_geometry_qualification.tsv','panel_coverage.tsv','ACCESS_QUALIFICATION.json','PUBLIC_FILE_PROBES.json','ARCHIVE_ROUTE_CHECK.json']]
records=[{'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files]
(O/'DELIVERY_MANIFEST.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8');files.append(O/'DELIVERY_MANIFEST.json')
zipname=O/'CGNL1_spatial_mainline_review_20261016_v1.zip'
with zipfile.ZipFile(zipname,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for p in files:z.write(p,str(p.relative_to(O)))
with zipfile.ZipFile(zipname) as z:assert z.testzip() is None
status={'utc':now,'stage':'validation_and_report_complete_pending_user_review','exit_status':'all_required_final_stages_success','mainline_decision':'not_qualified_as_established_primary_story','candidate_question':json.loads((O/'MAINLINE_DECISION.json').read_text())['candidate_question_zh'],'report':str(source),'review_zip':str(zipname),'zip_sha256':hashlib.sha256(zipname.read_bytes()).hexdigest(),'report_sha256':hashlib.sha256(b).hexdigest(),'raw_data_scripts_logs':str(O),'original_project_outputs_unchanged':True,'whole_project_complete':False,'external_messages_sent':False}
(O/'STATUS.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8');(O/'VALIDATION_COMPLETE.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8');(O.parent/'CGNL1_CURRENT_STATE_20261016_SPATIAL_REVIEW.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'report_bytes':len(b),'package_files':len(files),'zip_bytes':zipname.stat().st_size,'zip_sha256':status['zip_sha256'],'report_sha256':status['report_sha256'],'status':status['stage']},ensure_ascii=False))
