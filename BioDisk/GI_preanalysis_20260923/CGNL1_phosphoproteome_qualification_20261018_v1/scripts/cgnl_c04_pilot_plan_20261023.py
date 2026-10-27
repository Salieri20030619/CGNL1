from pathlib import Path
import json,re,hashlib,datetime,xml.etree.ElementTree as E
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'C04_sources'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,z):
 p.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert json.loads(p.read_text())==z
utc=datetime.datetime.now(datetime.timezone.utc).isoformat()
params={x.attrib['name']:x.text for x in E.parse(D/'PAK_params.xml').getroot() if 'name' in x.attrib}
comments=params['dataset.comments'];part=comments.split('One TMT-set was run for phosphopeptide mapping.')[1]
channels=dict(re.findall(r'sample (\d+) \((\d+[nc])\)',part))
records=[]
for block in (D/'GSE95705_samples_quick.soft').read_text().split('^SAMPLE = ')[1:]:
 gsm=block.splitlines()[0];title=re.search(r'!Sample_title = (.+)',block).group(1)
 ids=set(re.findall(r'_JL(\d+)\.(?:CEL|chp)\.gz',block));assert len(ids)==1
 sid=ids.pop();group='control' if 'negative control' in title else 'inflamed' if 'severely inflamed' in title else None
 assert group and sid in channels
 records.append(dict(sample=sid,gsm=gsm,title=title,group=group,phospho_TMT_channel=channels[sid],mapping_evidence='GEO supplementary CEL/CHP filename JL identifier joined to MassIVE sample identifier'))
assert len(records)==8 and len(set(x['sample'] for x in records))==8
assert sum(x['group']=='control' for x in records)==5 and sum(x['group']=='inflamed' for x in records)==3
mapping=dict(utc=utc,records=records,sources=[dict(path=str(D/n),sha256=sha(D/n)) for n in ['PAK_params.xml','GSE95705_samples_quick.soft']],cross_source_link_inferred_by_exact_sample_ID=True,biological_sample_count=8,not_independent_of_MK2_TCT_reanalysis=True)
save(D/'PAK_BIOLOGICAL_SAMPLE_MAPPING.json',mapping)
download=json.loads((D/'PAK_RAW_PILOT_DOWNLOAD.json').read_text());assert download['http']=='429' and not download['size_matches']
plan=dict(utc=utc,status='design_qualified_download_rate_limited_not_started',accession='MSV000081198',
 question='Does the archived colitis phosphoproteome contain reliably localized Cgnl1 sites with usable MS3 reporter measurements?',
 stage1='One preselected file phospho_fraction_1.raw; conversion/acquisition/scan-linkage feasibility only. Absence of Cgnl1 in one fraction is not a cohort negative.',
 stage1_gates=['download HTTP success and expected size 221702430 plus hash and RAW reader acceptance','isolated pinned ThermoRawFileParser >=2 self-contained Linux build; archive/license/version receipt','preserve MS1, MS2 and MS3, verify high/low-resolution MS2 and MS3 linkage and reporter channels','predefine full mouse proteome target-decoy identification and site localization; do not search only Cgnl1 and treat resulting score as proteome-wide FDR'],
 stage2='Only after stage1 passes: all phospho fractions for target detection; matched totalprotein data only if Cgnl1 quantification and sample linkage are usable.',
 decision_rules=['No group-effect testing before sample identity, localization, target detection and source scale qualify','n=8 biological samples: exploratory direction/effect/individual consistency; no new confirmatory claim','No CGNL1 mechanism claim from pathway enrichment or MS3 technical spectra','Do not claim 0.75 isolation specificity in original Methods is a phosphosite localization probability','No full raw archive download until single-file feasibility gates pass'],
 current_blockers=['single RAW download returned HTTP429,227byte HTML not221702430byteRAW','server PATH lacks dotnet/mono/raw converter/search engine;self-contained converter is available upstream but not yet installed'],
 software_primary_documentation='https://github.com/CompOmics/ThermoRawFileParser',new_tests=0)
save(D/'PAK_RAW_PILOT_SPEC.json',plan)
comparison=[
 dict(id='C01',question='PDAC S198 versus S297 relative remodeling',priority='conditional co-lead',strength='within-patient pair pattern and robust mean direction; purity proxies null; method-level localization recovered',gap='no independent PDAC confirmation; exact observation/localization provenance unresolved',next='measurement validation then same-disease independent targeted phosphoproteomics'),
 dict(id='C02',question='S486 regulation relative to protein and S283 across liver/biliary disease',priority='conditional co-lead',strength='source peptide reconstruction and60paired iCCA samples;relative direction also supported inH01',gap='not iCCA-specific;source-level localization probabilities and causal functional evidence absent',next='confirm S486 localization then compare iCCA/HCC within compatible prospective design'),
 dict(id='C03',question='HCC site pattern consistency and heterogeneity',priority='supporting comparative branch',strength='relativeS283-S297 negative mean acrossthree sources',gap='H02primary test null;H03imputation;unknown donoroverlap;not3independentreplications',next='donor/source clarification before pooled inference'),
 dict(id='C04',question='nonmalignant intestinal inflammation and barrier-associated Cgnl1 phosphorylation',priority='material qualification branch',strength='new8sampleidentity/channelmapping;matchedproteinandphospho design',gap='target detection unknown;processedtable inaccessible;rawpilotHTTP429',next='bounded raw conversion feasibility pilot before biological analysis')]
save(O/'CANDIDATE_COMPARISON_v1.json',dict(utc=utc,selection_basis=['measurement reliability','independent replication','disease interpretation','novelty and functional feasibility'],candidates=comparison,final_unique_mainline_selected=False,not_ranked_by_minimum_P=True))
report=(O/'CGNL1_CONTINUOUS_PROGRESS_20261023_v3.txt').read_text()
report+='''

新增：非肿瘤炎症分支样本身份闭合与候选比较（2026-10-23）
C04取得MassIVE原始params.xml（7589字节）及GEO GSE95705的完整quick元数据（25997字节）。按GEO CEL/CHP文件中的JL编号与MassIVE样本号精确连接，8例无冲突：
对照410/411/426/431/432，对应磷酸化TMT 127n/127c/129n/130n/130c；炎症421/422/430，对应128n/128c/129c。原文5对照、3炎症，与映射一致。这是跨库编号连接证据，尚不等于已确认Cgnl1被测到。
原始实验是MS2鉴定、SPS-MS3定量，并为磷酸化样本采集高/低分辨率MS2。原文0.75指标是前体分离特异性，不是磷酸位点定位概率。重处理必须保留MS3关联，不能套用只有MS2报告离子的定量流程。
已冻结单文件试验：仅用预先选定phospho_fraction_1.raw检查下载、转换与扫描关联；通过后才考虑全部磷酸化分馏。单个分馏没检出Cgnl1不能当作队列阴性。8个生物学样本保持探索性，技术扫描不能充当重复。
当前下载返回HTTP429，内容只有227字节HTML，已记录为失败而非RAW；服务器尚无转换/搜库工具。已确认上游提供Linux自包含转换器，暂未安装。没有对失败链接连续重试，没有启动全量重处理。

当前候选判断：
C01 胰腺癌S198/S297：条件性共同主候选。发现效应与稳健性有依据，仍缺同疾病独立验证和逐观测测量证据。
C02 S486及其相对S283模式：条件性共同主候选。60对iCCA和可核对源肽段是优势，HCC同向支持更适合跨疾病问题；不能称iCCA特异，也未证明功能因果。
C03 HCC位点模式：保留为比较分支，不能隐藏H02阴性、H03插补及队列重叠问题，不合并成三次独立验证。
C04 炎症/屏障：样本与实验设计资格前进一层，目标蛋白/位点检出仍未知；需要先完成原始数据试验。
这些排序依据测量可靠性、独立重复、疾病解释和功能验证可行性，没有按最小P值排序，也没有选定唯一主线。下一阶段优先完成C04小规模可行性验收，并保留C01/C02的测量与独立验证缺口。Goal仍active。
'''
f=O/'CGNL1_CONTINUOUS_PROGRESS_20261023_v4.txt';f.write_text(report,encoding='utf-8');assert f.read_text()==report
s=json.loads((O/'STATUS.json').read_text());s.update(utc=utc,C04='8 biological samples linked to group and MS3 TMT channel;target unknown;singleRAW pilot HTTP429;spec frozen',latest_report=str(f),next_task='C04 isolated RAW converter and bounded single-file pilot when rate limit permits;C01 exact matrix provenance remains unresolved');save(O/'STATUS.json',s)
g=json.loads((O/'GOAL_PLAN.json').read_text());g['utc']=utc
for task in g['tasks']:
 if task['id']=='T04':task['progress']='8sampleidentity/channelmapping complete;singlefile pilot spec frozen,download HTTP429'
 if task['id']=='T06':task['status']='interim_comparison_complete';task['artifact']=str(O/'CANDIDATE_COMPARISON_v1.json')
save(O/'GOAL_PLAN.json',g)
save(D/'PILOT_PLANNING_ACCEPTANCE.json',dict(utc=utc,mapping_unique=True,n_controls=5,n_inflamed=3,source_ids_match=True,new_tests=0,raw_download_accepted=False,report=str(f),report_sha256=sha(f),whole_project_complete=False))
c=P/'metadata/PROJECT_CHECKPOINT_CONTINUOUS_20261023.json';z=json.loads(c.read_text());z.update(utc=utc,latest_report=str(f),C04_pilot_spec=str(D/'PAK_RAW_PILOT_SPEC.json'));save(c,z)
print(json.dumps(dict(report=str(f),sha256=sha(f),samples=records)))
