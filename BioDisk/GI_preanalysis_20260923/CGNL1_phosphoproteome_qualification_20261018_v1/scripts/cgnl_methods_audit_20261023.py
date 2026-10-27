from pathlib import Path
import json,hashlib,datetime,csv
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
O=P/'continuous_exploration_20261023_v1';D=O/'source_methods'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,z):
 p.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
 assert json.loads(p.read_text(encoding='utf-8'))==z
t=(D/'Data_S1.txt').read_text();pages=t.split('\f')
assert (D/'Data_S1.pdf').read_bytes().startswith(b'%PDF')
assert sha(D/'Data_S1.pdf')=='02f3dcc06106de18d830f5006906a81f1a4ab8fc003c72eb007d435ed77e4143'
extract=[]
for i,page in enumerate(pages):
 if any(k in page for k in ['PTMProphet','Abundance table preprocessing','To formally impute','Gene and phosphosite intensities']):
  extract.append(dict(pdf_page_1based=i+1,text=page))
save(D/'METHOD_PAGE_EXTRACTS.json',extract)
coverage=[]
for p in sorted((P/'raw').glob('PDAC_phospho_site_abundance*')):
 with p.open() as f:
  reader=csv.reader(f,delimiter='\t');h=next(reader)
  for r in reader:
   if r[0].startswith('ENSG00000128849.'):
    assert len(r)==len(h)
    coverage.append(dict(file=p.name,site=r[0].split('|')[2],samples=len(h)-1,missing=sum(x in ['NA','NaN','nan',''] for x in r[1:]),sha256=sha(p)))
z=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),doi='10.1016/j.ccell.2023.06.009',
 source_url='https://ars.els-cdn.com/content/image/1-s2.0-S1535610823002192-mmc1.pdf',source_sha256=sha(D/'Data_S1.pdf'),
 classification='partial_resolution_of_source_qualification_not_observation_level_validation',
 established=['Michigan pipeline uses PTMProphet; site-level reports define confidently localized sites by probability >=0.75',
 'BCM harmonization uses Michigan gene and phosphosite intensities, reference-intensity median centering, and isoform-based site reannotation',
 'BCM composite identifier final flag 1 denotes primary selected protein isoform, not localization probability',
 'Methods describe reference-normalized log2 PSM ratios and weighted rollup',
 'PDAC table lists 36 proteome and 11 phosphosite outlier protein/site-TMTplex sets; no ComBat correction for either PDAC matrix',
 'Methods also describe DreamAI imputation on >=50% covered features; workflow presence alone does not identify the exact exported matrix as imputed',
 'Current CGNL1 S198 and S297 both retain 9 missing values among 44 normals; tumor missing counts 22/105 and 27/105 respectively'],
 inference='Coverage is inconsistent with a fully filled DreamAI version of these high-coverage rows, but retained NA alone does not prove every reported value was directly observed.',
 unresolved=['Exact export-stage provenance and observation mask for data_freeze_v1.2_reorganized',
 'Sample/PSM-level localization probabilities for C01 S198/S297',
 'Independent biological PDAC cohort confirmation'],
 coverage=coverage,new_tests=0,original_results_unchanged=True)
save(D/'HARMONIZED_METHODS_AUDIT.json',z)
report=(O/'CGNL1_CONTINUOUS_PROGRESS_20261023_v2.txt').read_text(encoding='utf-8')
report+='''

第五轮证据补充：统一重处理方法已取得（2026-10-23）
本补充更新前文“定位/处理证据仍需补齐”的状态，不改变任何模型结果。
原文：Proteogenomic Data and Resources for Pan-Cancer Analysis，Cancer Cell 2023，DOI 10.1016/j.ccell.2023.06.009。
出版社补充方法 Data S1 已下载并通过 PDF 格式及 SHA256 校验。Michigan 流程用 PTMProphet 定位，位点级报告以定位概率≥0.75为可靠定位标准；BCM统一流程使用其结果，进行参考强度中心化及蛋白异构体坐标重注释。因此，C01已有流程层面的定位依据；CDAP缺少S198/S297不能等同于统一流程结果错误。两个导出仍来自同一生物学队列，不能计为独立验证。
复合位点ID末尾的“1”表示选定的主要蛋白异构体，不是定位概率1，也不是100%定位置信度。
补充方法同时描述异常位点-TMTplex数据置NA、部分癌种ComBat和DreamAI插补。PDAC表中列出36个蛋白及11个位点层面的异常TMTplex集合，蛋白和磷酸化矩阵均未列ComBat。不能把跨癌种的方法步骤全部套到PDAC当前文件。
当前使用的data_freeze_v1.2_reorganized矩阵中，S198/S297在44例正常组织中各缺失9例，在105例肿瘤中分别缺失22和27例。高覆盖行仍保留缺失，不符合这些行已被完整填补的版本；但仅凭NA不能证明每个现有数值都是直接观测。精确导出阶段与逐观测掩码仍未核实，不能把“未见完整填补”写成“已证明完全无插补”。
候选排序暂不变：C01胰腺癌S198−S297与C02 S486相对模式继续保留；C03保留HCC异质性和关键阴性；C04仍缺可分析的目标位点材料。此次没有新增检验，没有更换主终点，没有按接近原论文或更小P值重新选择模型。
下一步：核对精确矩阵的导出说明/观测掩码与可获取的位点级定位文件，再决定C01是否具备进入疾病相关功能验证的材料条件。Goal保持active，机制与最终主线均未宣告完成。
'''
f=O/'CGNL1_CONTINUOUS_PROGRESS_20261023_v3.txt';f.write_text(report,encoding='utf-8')
assert f.read_text(encoding='utf-8')==report
s=json.loads((O/'STATUS.json').read_text());s.update(utc=z['utc'],source_methods='PTMProphet>=0.75 method recovered;exact export imputation status and per-observation localization pending',latest_report=str(f),next_task='Exact harmonized export-stage provenance and observation mask; C01 targeted localization evidence')
save(O/'STATUS.json',s)
save(D/'ACCEPTANCE.json',dict(utc=z['utc'],source_pdf_valid=True,source_sha256=z['source_sha256'],coverage_rows=len(coverage),new_tests=0,report=str(f),report_sha256=sha(f),goal_status='active',whole_project_complete=False))
c=P/'metadata/PROJECT_CHECKPOINT_CONTINUOUS_20261023.json'
x=json.loads(c.read_text());x.update(utc=z['utc'],latest_report=str(f),source_methods_audit=str(D/'HARMONIZED_METHODS_AUDIT.json'),goal_status='active',whole_project_complete=False);save(c,x)
print(json.dumps(dict(report=str(f),sha256=sha(f),coverage_rows=len(coverage),method_pages=[x['pdf_page_1based'] for x in extract])))
