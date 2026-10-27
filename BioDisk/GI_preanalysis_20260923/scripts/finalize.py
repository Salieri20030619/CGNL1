import pathlib,json,csv,hashlib,zipfile,importlib.util,re,html
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
p=R/'05_reports'/'feasibility_report_zh.md';s=p.read_text(encoding='utf-8-sig');s=s.replace('WGS原始数据公开不等于本轮已完成基因组整合。','WGS原始数据公开不等于本轮已完成基因组整合。可用临床表也没有完整的逐样本用药史、取材时间及实测肿瘤纯度，原发/转移差异仍需结合这些因素解释，不能直接写成时间转换或因果机制。');p.write_text(s,encoding='utf-8')
rows=[
 ['CRC','GSE145356','25 libraries / one donor for principal contrasts','downloaded_and_analyzed','CPT DTP/regrowth descriptive RNA effects','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145356'],
 ['CRC','GSE294559','9 patients in study; P001 TP/TO pilot analyzed','2_h5_downloaded_and_analyzed','RNA/ATAC background, no treatment timecourse','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294559'],
 ['CRC','Zenodo11058999','two principal PDO models; 37.4GB package','metadata_only','not downloaded; no independent human recovery validation','https://zenodo.org/records/11058999'],
 ['CRC','GSE163171;GSE187650;GSE163035;GSE187512','mouse/model perturbation','metadata_and_some_SOFT_tables_downloaded','not matched human chemotherapy recovery timecourse','https://www.nature.com/articles/s43018-022-00402-0'],
 ['GC','GSE234129','19488 cells; 6 patients; 17 samples; zero malignant/epithelial annotations','downloaded_and_analyzed','immune/stromal reference only','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234129'],
 ['GC','GSE200697','3 pairs / 6 arrays','downloaded_and_analyzed','probe level only; platform lacks symbols','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200697'],
 ['GC','PXD064711','6847 protein rows; 25 coded channels + 2 refs','protein_table_and_metadata_downloaded','group and replicate mapping not resolved','https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD064711'],
 ['GC','GSE314812','26 patients / 55 RNA libraries; 17 paired patients selected','downloaded_and_analyzed','paired primary vs peritoneal tissue RNA, not ascites transition','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314812'],
 ['GC','HRA002712','35 patients in paper','controlled_raw_not_obtained; supplementary/source_sheets_checked','no complete gene-cell matrix in obtained source sheets','https://www.nature.com/articles/s41467-023-36310-9'],
 ['GC','GSE163558','one peritoneal sample','metadata_only','insufficient peritoneal replication','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163558'],
 ['GC','GSE228598','28 fluid samples mixing ascites/washing','metadata_only','not a paired three-site transition design','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228598'],
 ['GC','EGAS00001004443;EGAS00001006397;EGAS50000000501','controlled genomic/transcriptomic studies','literature_data_availability_review','not counted as downloaded expression matrices','https://www.nature.com/articles/s41591-020-1125-8']]
with (R/'01_metadata'/'dataset_decision_manifest.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['direction','accession','design','acquisition_status','role_or_limitation','source_url']);w.writerows(rows)
(R/'00_docs'/'DECISION.md').write_text('# Decision, 2026-09-23\n\nCRC: original public-multiomics-led regrowth regulatory project is not ready; model exploration is feasible.\nGC: conditionally proceed after narrowing to lipid-uptake-associated peritoneal microenvironment.\nPriority: GC before CRC. No target or publication outcome is confirmed.\n\nGSE314812 added after initial candidate audit: 17 independent paired patients, primary model adjusted for preservation; composition, frozen-only and serial-patient sensitivity recorded.\n')
# Compact navigable HTML overview with full report link, results and figures.
body='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>GI 多组学预分析</title><style>body{max-width:1050px;margin:40px auto;padding:0 24px;font:17px/1.7 system-ui,"Microsoft YaHei",sans-serif;color:#243746}h1{color:#12344c}img{max-width:100%;border:1px solid #ddd;margin:16px 0}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ddd;padding:10px;text-align:left}a{color:#086788}.note{background:#eef6f7;padding:20px}</style><h1>消化系统肿瘤方向1、2：预分析</h1><p>2026-09-23 · 已完成下载、设计审查、实际矩阵分析与必要验收</p><div class="note"><b>更新优先级：胃癌方向2 ＞ 结直肠癌方向1。</b><p>方向2有17对患者的组织转录组基础，脂质摄取相关信号受组织组成影响，应收窄为微环境问题；方向1缺少多患者治疗恢复期调控组学，原题暂缓作为首选。</p></div><p><a href="feasibility_report_zh.md">阅读完整中文报告（包括阴性结果、方法与限制）</a></p><table><tr><th>方向</th><th>已获得证据</th><th>立项条件</th></tr><tr><td>CRC残存细胞再生</td><td>同患者PDX模型中增殖抑制与恢复；RNA/ATAC参考文件可分析</td><td>增加独立模型和治疗恢复期调控组学；本轮没有患者层面的推断性统计</td></tr><tr><td>GC腹膜转移</td><td>17对患者，脂质摄取面板13/17同向，P=0.0010，FDR=0.0081</td><td>组成校正后FDR=0.180；需要确定细胞来源及完成其他组学/功能验证</td></tr></table>'''
for fn,caption in [('GC_03_panel_sensitivity.png','胃癌：主分析与组成敏感性；置信区间按模型自由度计算。'),('GC_04_paired_lipid_uptake.png','胃癌：患者配对方向。组织信号不等同于癌细胞内机制。'),('CRC_02_panel_effects.png','结直肠癌：单患者模型描述性效应，不作患者推断。')]:body+=f'<p>{caption}</p><img src="../04_results/figures/{fn}">'
body+='</html>';(R/'05_reports'/'index.html').write_text(body,encoding='utf-8')
checks=list(csv.DictReader((R/'logs'/'validation_checks.tsv').open(),delimiter='\t'));assert all(x['pass']=='True' for x in checks)
status={'status':'complete','scope':'preanalysis and feasibility, not full mechanism validation','date':'2026-09-23','server_root':str(R),'CRC_decision':'defer original project; conditional model exploration','GC_decision':'conditional proceed with narrowed microenvironment hypothesis','priority':['GC','CRC'],'checks_passed':len(checks),'analysis_exit_codes':{'CRC_R':0,'preanalysis_python':0,'GC_bulk_R':0,'figures_verification':0},'unfinished_analysis_jobs':False}
(R/'logs'/'TASK_COMPLETE.json').write_text(json.dumps(status,indent=2))
# Exclude raw/matrices and temporary partials from review package. Preserve full server inputs.
package=R/'05_reports'/'GI_preanalysis_review_package.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for top in ['00_docs','04_results','scripts','logs']:
  for f in (R/top).rglob('*'):
   if f.is_file() and '__pycache__' not in str(f) and not f.name.endswith('.part'):z.write(f,str(f.relative_to(R)))
 for f in [p,R/'05_reports'/'index.html',R/'01_metadata'/'dataset_decision_manifest.tsv',R/'01_metadata'/'file_integrity_manifest.tsv',R/'01_metadata'/'GSE314812_sample_design.tsv',R/'01_metadata'/'GSE314812_paired_design.tsv']:
  z.write(f,str(f.relative_to(R)))
print(json.dumps(status));print('PACKAGE_BYTES',package.stat().st_size);print('PACKAGE_SHA256',hashlib.sha256(package.read_bytes()).hexdigest())
