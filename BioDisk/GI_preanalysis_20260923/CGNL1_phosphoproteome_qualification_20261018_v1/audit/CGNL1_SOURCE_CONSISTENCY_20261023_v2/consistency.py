from pathlib import Path
import json,re,hashlib,datetime,zipfile,xml.etree.ElementTree as ET
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
O=P/'audit/CGNL1_SOURCE_CONSISTENCY_20261023_v2'
S=O/'sources'
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576')
ns={}
exec((P/'audit/CGNL1_FOLLOWUP_20261023_v1/followup.py').read_text(encoding='utf-8').split('baseline=read')[0],ns)
rows,sheets=ns['rows'],ns['sheets']
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def read(f):return json.loads(f.read_text(encoding='utf-8'))
def save(name,d):
 f=O/name;f.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert read(f)==d
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
old=P/'audit/CGNL1_FOLLOWUP_20261023_v1'
before={str(f):sha(f) for f in old.rglob('*') if f.is_file()}
for path,h in read(old/'SOURCE_HASHES_BEFORE.json').items():assert sha(Path(path))==h
save('BASELINE_HASHES_BEFORE.json',before)

# Three independently available layers: original design, repository samples, clinical table.
design=read(P/'audit/CGNL1_FIRST_REVIEW_20261022_v1/B02_TMT_CHANNEL_MAP.json')
sam=read(S/'NODE_samples.json')['data']['sapTableData']
assert sam['totalElements']==len(sam['content'])==530 and sam['totalPages']==1
clinical_rows=list(rows(D/'source_search_593/B02_TableS1.xlsx','S1A. Clinical info'))
head=clinical_rows[1][1]
clinical={str(c.get('A'+str(n))):{head.get(re.sub(r'\d+$','',k)+'2',k):v for k,v in c.items()} for n,c in clinical_rows[2:] if c.get('A'+str(n)) is not None}
assert len(clinical)==262
figure_pairs={'587':['iCCA_16','iCCA_22'],'647':['iCCA_02','iCCA_22'],'321':['iCCA_06','iCCA_22'],'133':['iCCA_10','iCCA_22'],'823':['iCCA_18','iCCA_22']}
matches=[]
for case,batches in figure_pairs.items():
 found=[x for x in sam['content'] if x['name']==case]
 assert len(found)==1
 x=found[0];a=x['attributes'];c=clinical[case]
 size=float(re.search(r'Tumor_size_diameter:([^;]+)',a['disease_phenotype']).group(1))
 checks={'age':float(a['age'].removesuffix('y'))==c['Age'],'sex':a['sex'].lower()==c['Sex'].lower(),
 'tumor_size':abs(size-c['Tumor_size_diameter (cm)'])<1e-10,'TNM_stage':a['development_stage']==c['TNM stage']}
 assert all(checks.values())
 actual=sorted(x['batch'] for x in design if x['case']==case)
 assert actual==sorted(batches)
 matches.append(dict(case=case,batches=actual,repository_sample=x['sapNo'],repository_records=1,
  clinical_checks=checks,figure_S1C_pair_matches_design=True,
  date_present=bool(a.get('sample_collection_date')),site_detail_present=bool(a.get('sampling_site') or a.get('anatomical_site')),
  tissue=x['tissue'],pathological_nature=a['pathological_nature']))
assert '227' not in clinical
s227=next(x for x in sam['content'] if x['name']=='227')
save('B02_REPEAT_IDENTITY_CONSISTENCY.json',dict(utc=now,issue='CGNL-A04',records=matches,matched_clinical_fields=20,
 conclusion='同一登记患者/肿瘤样本跨批次复测用于平台QC，受到原文、S1C、通道表、NODE样本和临床字段一致支持。',
 evidence_level='strong_source_concordance_for_same_patient_sample_remeasurement',
 not_claimed=['未确定是独立提取、独立酶解、重复标记或其后步骤','未取得同一物理组织块/提取物的唯一标识','未得到最终重复值合并/选择规则','没有原始通道值，不声称重新算出复测相关性或一致性系数'],
 longitudinal_or_multiregion='当前公开材料未显示这些重复来自不同时间点或不同位置；不能将字段缺失当作绝对排除。主文初次未治疗切除和粉碎分装设计支持技术复测解释。',
 case227=dict(repository_sample=s227['sapNo'],in_repository=True,in_design=True,in_clinical_registry=False,in_final_quantification=False,
  repository_sex=s227['attributes'].get('sex'),repository_age=s227['attributes'].get('age'),exclusion_reason='未公开，不能由缺失反推具体原因'),
 phenotype_note='主文另述五对肝胆癌及邻近非肿瘤组织的benchmark，设计/通道不同，不能与这5个重复患者自动视为同一组。',
 modeling_remains_paused=True))
experiments=read(S/'NODE_proteomic_experiments.json')['data']['expTableData']['content']
save('SOURCE_METADATA_DISCREPANCIES.json',dict(issue='CGNL-A09',severity='medium',classification='source_metadata_inconsistency',
 findings=[dict(experiment=x['expNo'],name=x['name'],quantification=x['attributes']['quantification'],instrument_type=x['attributes']['instrument_type'],platform=x['attributes']['platform']) for x in experiments],
 conflict='NODE写TMT 10plex及ESI-QUAD-TOF；主文为11plex/Q Exactive HF-X，实际通道表包含10个肿瘤通道和131C混合参照。',
 impact='仓库概括字段不能覆盖主文/通道级实验设计；不足以单独推翻已验证的最终矩阵。',
 resolution='后续以主文和实测通道表建立设计，保留仓库原字段。若重处理原始谱图，须先核对实际仪器和搜索参数。'))

# H01 reports are genuine historical source-linked QC, not a missing version certificate.
qc=[]
audit=read(P/'audit/CGNL1_FIRST_REVIEW_20261022_v1/H01_NORMALIZATION_AUDIT.json')
for study,kind,date in [('PDC000198','Proteome','2019-07-31'),('PDC000199','Phosphoproteome','2019-07-15')]:
 receipt=read(S/(study+'_qc_receipt.json'));assert receipt['md5_pass'] and receipt['exit']==0
 f=S/receipt['name'];assert sha(f)==receipt['sha256']
 t=(S/(f.stem+'.visible.txt')).read_text(encoding='utf-8')
 expected='Zhou_Liver_Cancer_'+kind+('.tmt11.tsv' if kind=='Proteome' else '.phosphosite.tmt11.tsv')
 assert expected in t and 'Biological Samples 330' in t and 'Analytical Samples 33' in t
 qc.append(dict(study=study,report=receipt['name'],date=date,input_quantification=expected,biological_sample_columns=330,
  TMT_plexes=33,official_MD5_verified=True,exact_pipeline_commit_in_visible_report=False,
  explicit_normalization_recipe_in_visible_report=False))
save('H01_EXPORT_LINEAGE_CONSISTENCY.json',dict(utc=now,issue='CGNL-A01',reports=qc,
 judgment='QC报告明确指向当前使用的导出矩阵，330样本列和33批设计一致；增加历史来源证据，但未给出确切归一化版本或参数。',
 warning='报告中的330个Biological Samples是组织/通道列，不是330名独立患者；本项目配对分析仍为159名患者。',
 normalization_decision='沿用前轮双估计量并列约束；不因找到QC报告而取消归一化敏感性。'))

# Cohort chronology: verify actual Excel date system, distinguish acquisition-file dates.
h01=next((D/'H01').glob('*.xlsx'))
removal=list(rows(h01,'Table 4. PHI removed'))
assert any('Date of surgery removed' in str(c) for _,c in removal)
h03=P/'raw/intermittent_564/HCC_verified/mmc2.xlsx'
with zipfile.ZipFile(h03) as z:
 wp=ET.fromstring(z.read('xl/workbook.xml')).find('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}workbookPr')
 assert wp is None or wp.attrib.get('date1904','0') in ['0','false']
rr=list(rows(h03,'2. Detialed Clinicopathological'))
assert rr[1][1]['AD2']=='Date of surgery'
dates=[c['AD'+str(n)] for n,c in rr[2:] if isinstance(c.get('AD'+str(n)),(int,float))]
assert len(dates)==160
dt=lambda x:(datetime.datetime(1899,12,30)+datetime.timedelta(days=x)).date().isoformat()
f02=next((D/'H02/global_589').rglob('*.xlsx'))
h02_headers=dict(rows(f02,'Quantitative analysis (222)',{1}))[1]
assert 'Raw_file' in h02_headers.values() and not any('surgery' in str(x).lower() for x in h02_headers.values())
overlap=dict(utc=now,issue='CGNL-A06',status='source_evidence_extended_patient_overlap_unresolved',
 H01=dict(file=str(h01),sha256=sha(h01),surgery_dates_explicitly_removed=True,hospital_accession_masked=True,PHI_removal_sheet=[c for _,c in removal]),
 H02=dict(file=str(f02),sha256=sha(f02),columns=h02_headers,rule='Raw_file日期不解释为取材日期；当前跨表映射仍不包含跨研究稳定患者标识。'),
 H03=dict(file=str(h03),sha256=sha(h03),date_field='2. Detialed Clinicopathological!AD2',date_system='Excel1900',dated_patients=160,
  earliest_surgery=dt(min(dates)),latest_surgery=dt(max(dates)),scope='整套160例临床队列，不替代118例位点完整子集的样本量'),
 decision='H03手术日期已补齐；H01日期被明确移除。因此无法按当前公开表证明H01/H02/H03无重叠，也没有证据可断言已存在重叠。',
 next_required='作者或数据管理方的去标识化跨研究映射/重叠声明；不用年龄性别相似度猜测患者身份。')
save('HCC_COHORT_CHRONOLOGY_CONSISTENCY.json',overlap)
report='''CGNL来源与一致性核查 v2（2026-10-23）

本轮仅来源核查，B02继续暂停建模。此前结果、统计口径、原始数据及v1审核均保留。

一、B02重复性质：已有较强证据支持同一患者样本的跨批次复测，用于定量平台QC。
新取得的原文补充图S1C红字标明5组重复：587（16/22批）、647（02/22批）、321（06/22批）、133（10/22批）、823（18/22批）。图示每组相邻聚类，作者将其用于平台定量稳定性评价。5组与原始通道表完全对应。
NODE公开样本元数据完整取得530条，5个ID各只有一个登记的肿瘤样本记录。其年龄、性别、肿瘤大小、TNM分期与原文临床表逐项一致，共20/20项。主文说明患者来自首次、治疗前切除，组织粉碎后分装。
因此当前最符合证据的是同一登记患者/肿瘤样本的技术或制备复测；没有找到这些5组来自不同时间点或不同取材位置的标注。但尚不能细分重复发生在组织块、提取、酶解、标记或上机哪一层，也没有最终患者列合并或选择规则。不能简单称为“同一提取物再次上机”。图中相邻聚类是作者提供的全局蛋白QC证据，不是本轮重新计算的通道相关性，更不证明S283无批次偏差。
病例227：确实存在于NODE和设计表，但不在262例最终临床表及214例最终矩阵；NODE年龄、性别为not collected。具体排除原因仍不明。
主文另有五对肝胆癌/邻近非肿瘤benchmark，其通道设计不同，不能与这5个重复患者自动混为一组。

二、发现的来源差异（CGNL-A09）。
NODE实验元数据写TMT 10plex、ESI-QUAD-TOF；主文及实际设计表为11plex、Q Exactive HF-X，包含10个肿瘤通道和131C混合内参。该差异保留为仓库元数据不一致；后续设计解释采用主文与通道表，不静默改原字段，也不据此宣布最终矩阵无效。

三、H01来源一致性。
新获取PDC000198和PDC000199官方QC报告，大小、官方MD5及SHA256验证通过。分别为2019-07-31和2019-07-15历史报告，明确引用本项目使用的蛋白/磷酸位点矩阵，330个样本列、33个TMT批次一致。报告不提供精确处理commit或完整归一化配方，A01仍部分解决。330是组织/通道列数，不是330位患者；配对分析仍为159位患者。保留源口径与额外中心化口径，不选有利方向。

四、HCC队列独立性。
H03原始Table S1的AD列明确为Date of surgery，核实Excel日期系统后，160例手术日期范围为2014-02-18至2019-06-14。该信息补充了前轮“日期未取得”的状态。
H01原始工作簿Table 4明确说明手术日期等信息因PHI保护被移除，医院编号也已遮蔽。H02可用MaxQuant映射表提供Raw_file、Experiment、BCLC_stage及Sample_ID；不能把质谱文件名前缀日期当作手术日期。
目前无可靠的跨研究患者映射。因此既不能认证三队列完全独立，也不能断言已经发现患者重叠。若要彻底解决，需数据管理方提供去标识化跨研究映射或明确重叠说明；不尝试用临床特征拼接真实身份。

结论：B02“重复用途和同一患者样本对应关系”已有多来源一致证据；“重复发生层级及最终数值合并规则”仍未解决。H01历史导出来源得到加强，归一化口径问题保留。H03日期得到补充，但三队列患者级排重仍无法由现有公开字段完成。本轮未拟合临床模型、未增加推断检验。

主要来源：
B02补充图：https://ars.els-cdn.com/content/image/1-s2.0-S1535610821006590-mmc1.pdf（PDF第3页S1C）
B02主文：https://ueaeprints.uea.ac.uk/id/eprint/102354/1/1-s2.0-S1535610821006590-main.pdf（PDF第24-28页）
B02仓库：https://www.biosino.org/node/project/detail/OEP001105
H01：https://proteomic.datacommons.cancer.gov/pdc/study/PDC000198 及PDC000199
H03：https://pmc.ncbi.nlm.nih.gov/articles/PMC10783603/ 及原始Table S1
'''
f=O/'CGNL1_SOURCE_CONSISTENCY_v2.txt';f.write_text(report,encoding='utf-8');assert '\ufffd' not in f.read_text(encoding='utf-8')
assert all(sha(Path(k))==v for k,v in before.items())
assert all(sha(Path(k))==v for k,v in read(old/'SOURCE_HASHES_BEFORE.json').items())
save('CONSISTENCY_ACCEPTANCE.json',dict(utc=now,previous_files_unchanged=len(before),previous_source_workbooks_unchanged=4,
 source_samples_complete=530,repeat_cases=5,clinical_fields_agree=20,figure_design_agree=5,
 H03_surgery_dates_verified=160,clinical_models=0,new_inferential_tests=0,B02_modeling_paused=True,
 all_issues_resolved=False,whole_project_complete=False,
 report_sha256=sha(f),artifacts=[dict(path=str(x),sha256=sha(x)) for x in O.rglob('*') if x.is_file() and x.name!='CONSISTENCY_ACCEPTANCE.json']))
checkpoint=dict(utc=now,baseline_checkpoint='PROJECT_CHECKPOINT_595.json',latest_followup=str(O),acceptance=str(O/'CONSISTENCY_ACCEPTANCE.json'),
 status='source_consistency_followup_complete_remaining_provenance_gaps',B02_modeling_paused=True,whole_project_complete=False)
(P/'metadata/PROJECT_CHECKPOINT_CONSISTENCY_20261023_v2.json').write_text(json.dumps(checkpoint,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'report':str(f),'sha256':sha(f),'repeat_cases':5,'clinical_checks':20,'H03_dates':[dt(min(dates)),dt(max(dates))],'baseline_files_unchanged':len(before)},ensure_ascii=False))
