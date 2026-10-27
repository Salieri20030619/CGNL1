from pathlib import Path
import json,datetime,hashlib,shutil
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
O=P/'continuous_exploration_20261023_v1'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,d):
 f=O/name;f.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert read(f)==d
for w in ['wave1','wave2']:
 a=read(O/w/'ACCEPTANCE.json')
 for x in a['artifacts']:assert sha(Path(x['path']))==x['sha256']
w1=read(O/'wave1/RESULTS.json');w2=read(O/'wave2/RESULTS.json')
sources=[
 dict(id='C04_TCT',title='Integrated in vivo multiomics identifies p21-activated kinase signaling as a driver of colitis',pmcid='PMC6719711',dataset='MSV000081198',
  design='Rag1-null mouse T-cell transfer colitis; matching same-animal proteome/phosphoproteome;8 deposited sample IDs;group assignment still needs source-table verification',
  needed=['NIHMS1046950-supplement-Table_2.xlsx','NIHMS1046950-supplement-Table_3.xlsx'],
  status='processed_tables_exist_but_download_challenge_not_resolved',target_Cgnl1_measured=None,
  evidence='Valid full article HTML/BioC and MassIVE complete file inventory obtained. PMC attachment endpoints return HTML instead of XLSX; isolated headless browser remained at reCAPTCHA. EPMC supplementary API explicitly says article not open access. MassIVE inventory38files contains37raw plus params.xml, no processed abundance matrix.',
  action='Seek legitimate accessible publisher/author copy or assess raw-reprocessing feasibility; do not repeat identical blocked URLs or infer target absence.'),
 dict(id='C04_MK2',title='Substrate-based kinase activity inference identifies MK2 as driver of colitis',pmcid='PMC7208439',dataset='MSV000082217',
  design='TCT and TNF-deltaARE intestinal inflammation models; TCT overlaps earlier study, not a second independent TCT cohort',
  needed=['NIHMS1553501-supplement-Suppl_Table_S2.txt','NIHMS1553501-supplement-Suppl_Table_S3.txt','S4-S7 as needed after mapping'],
  status='processed_tables_listed_download_challenge',target_Cgnl1_measured=None,
  evidence='MassIVE inventory39files:37raw,params.xml,zero-byte statistics.tsv; published processed tables not present in that inventory. PMC txt downloads are HTML, not valid tables.'),
 dict(id='C04_PCDH20',title='Protocadherin20 maintains intestinal barrier function to protect against Crohn disease by targeting ATF6',pmcid='PMC10320870',
  design='WT vsPcdh20CKO mice,n3/group; paired omics stated in Methods,localization>0.75 and protein normalization described',
  status='assay_relevant_public_patient_or_sample_matrix_not_located',target_Cgnl1_measured=None,
  evidence='FullXML acquired and parsed. Listed XLSX supplements describe interactors/enrichment rather than sample phosphosite intensities. Data availability points to RNA and microscopy/blot data, no identified phosphosite matrix accession.'),
 dict(id='C04_UC_HUMAN',title='Phosphoproteomics Reveals Novel Insights into the Pathogenesis and Identifies New Therapeutic Kinase Targets of Ulcerative Colitis',doi='10.1093/ibd/izad291',pmid='38085663',
  design='Human UC rectal tissue and healthy controls,proteome plus class-I phosphosite profiling reported',
  status='relevant_paper_confirmed_processed_matrix_accession_pending',target_CGNL1_measured=None,
  evidence='EPMC core bibliographic record obtained; publisher access returned403. No repository accession confirmed in this round. Do not label target missing or substitute totalproteome-only studies.')]
save('C04_SOURCE_QUALIFICATION.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),records=sources,
 PRIDE_search='colitis query page0=100/page1=36 records inspected for phosphoproteomic candidates; not an exhaustive assay-specific or all-repository search',
 clinical_or_model_tests_performed=0,no_target_absence_claim=True))
plan=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),goal_status='active',goal_user_authorized=True,
 scope='CGNL1 phosphorylation across digestive diseases;C01/C02/C03/C04 candidates,not onlyS283',
 tasks=[dict(id='T01',status='complete',task='Wave1 seven locked contrasts/recalculations and independent QA'),
 dict(id='T02',status='complete',task='Wave2 locked cross-disease checks;7validtests and5non-evaluable/descriptive scenarios'),
 dict(id='T03',status='complete_for_available_inputs',task='Requalify Fudan,Britton,COAD for denominator-free exact pairs;none supplies current locked pairs'),
 dict(id='T04',status='active',task='Obtain and qualify C04 processed site+protein tables;four relevant sources tracked'),
 dict(id='T05',status='pending',task='Audit candidate-site localization/processing and assess composition alternatives using source-appropriate,predeclared checks'),
 dict(id='T06',status='pending',task='Compare candidates after independent-source and function/novelty assessment;do not select by smallestP')],
 rules=['Preserve all primary nulls and alternativeprocessing sensitivity;do not replace fullsample result with favorable sharedsubset',
 'No pooled HCC inference until overlap and measurement compatibility resolved',
 'No raw archive counts or running monitors as proof of scientific completion',
 'No author contacts sent without explicit authorization',
 'No identical access retries without new route/evidence;retain valid partials and distinguish access failure from biological absence'],
 next_actions=['Resolve C04 processed-matrix access using public alternative sources or a concrete raw-analysis feasibility pilot',
 'Qualify CGNL1 measured sites,biological sample mapping,localization and totalprotein before testing C04',
 'Perform targeted source/functional qualification for C01 S198-S297,C02 S486-S283,C03 S283-S297;no automatic clinicalendpoint sweep'],
 active_analysis_processes=0,whole_project_complete=False)
save('GOAL_PLAN.json',plan)
save('STATUS.json',dict(utc=plan['utc'],goal_status='active',wave1='complete_QA_pass',wave2='complete_QA_pass',
 C01='S198-S297 supported in discovery;S283-S708 primary null;no eligible independentPDAC exactpair source',
 C02='S486 effect andrelativeS486-S283 supported inB03;H01relativepattern also supported;notiCCA-specific proof',
 C03='relativeS283-S297 negative in3HCC cohorts;H01/H03 supported,H02P0.1089;GRR/overlap limits persist',
 C04='four relevant sources identified;sampleleveltarget qualification pending on processedtable access',
 next_task='T04/T05',whole_project_complete=False))
report='''CGNL1持续探究：两轮阶段结果与接续检查点（2026-10-23）

Goal状态：active。四个候选继续保留，未确定唯一主线，未把项目标记完成。本报告是阶段检查点，不是机制结题。

一、已执行并验收
第一轮7项预设对比/复核；第二轮7项有效跨疾病对比，另记录5项低覆盖或不可评估场景。两轮均以独立患者为单位，Python与R独立数值核验通过。所有原始和既有结果保持不变。新对比是在已知单个位点结果后提出，属于明确标注的探索性分析，不是盲法验证。

二、C01：PDAC多位点模式
S198-S297：28对，平均差-0.832048，95%患者bootstrap均值区间[-1.181905,-0.502391]，21降7升，精确符号秩P=0.000105701，C01家族FDR=0.000211403。逐例剔除均值保持负向，P范围0.0000474602至0.000208884。共同总蛋白分母代数抵消验证通过。
S283-S708：21对，平均差-0.191461，P=0.320457，FDR=0.320457，不支持预设位点间差异。可在替代CDAP处理中对应的14例子集呈负向及P<0.05，但原处理同一14例也更强；这是覆盖/子集敏感性，不能替换21例主要阴性结果，也不是独立生物学验证。
跨疾病：胃癌G02的S198-S297只有14对，平均-0.031944，P=0.714844，未验证该模式；另一胃癌无对应共同位点。不能由此宣布PDAC特异性。
独立PDAC再资格：Fudan公开补充表只有目标磷酸位点汇总FC/P，缺个体位点矩阵；Britton只含S283及多磷酸化组合，缺固定的S198/S297或S283/S708精确配对。解除总蛋白分母限制仍无法填补第二位点缺失，未改用有利替代位点。

三、C02：iCCA S486
已知S486蛋白校正效应在60例直接源行重建通过：均值-0.319315，P=0.000120399，FDR=0.000240798；不额外中心化仍为-0.396849，P=0.00373336。这个复核不计作新增独立队列。
新S486-S283对比：60对，均值-0.283105，95%均值区间[-0.511364,-0.058819]，P=0.0247147，FDR=0.0247147；逐例剔除方向不变且P范围0.0113433至0.0383385。此相对位点差值在两种中心化口径下数值一致至浮点误差。
跨疾病：H01同一对比129对，均值-0.235845，P=1.08458e-6，FDR=3.25374e-6；PDAC26对P=0.0593549，食管25对P=0.219986。H02/H03只有4/3对，按原规格仅描述。不能把C02写成已证明iCCA专属机制。

四、C03：HCC S283/S297相对模式
H01：34对，均值-0.894437，28降6升，P=1.94355e-6，FDR=2.91533e-6。
H02：17对，均值-0.596413，10降7升，P=0.108887，FDR=0.108887，未达到原始P阈值。逐例剔除P可从0.0335388到0.192810，保留全部17例作为主结果，不挑有利删样。
H03：27对，均值-1.366250，25降2升，P=6.40750e-7，FDR=1.92225e-6；但源GRR插补、定位与队列独立性边界仍在。
三队列相对差值均对已检查的共同中心化项不敏感，说明原单点方向冲突的一部分可通过区分估计量来理解。不能因此认定已排除全部技术/细胞组成因素。未做HCC队列间合并或分期交互检验。
跨疾病：PDAC35对S283-S297均值-0.229355，P=0.184250；胃癌G01/G02分别15/27对，均值+0.340853/+0.092353，P=0.0553589/0.238572，均无原始P<0.05支持。不能把方向或P值差别直接当癌种交互。

五、C04：非肿瘤分支仍在材料资格阶段
找到四个真正相关的蛋白组/磷酸化研究：T细胞转移结肠炎、TCT与TNFΔARE两模型、Pcdh20CKO肠炎以及人UC直肠组织。前两篇的处理表公开列出，但当前自动下载返回HTML验证页面而不是实际表；独立无界面浏览器也未取得附件。MassIVE两个完整目录分别38/39条，只有原始谱图和参数/空统计文件，没有可直接分析的定量矩阵。
Pcdh20研究Methods确认n3/group的蛋白/磷酸化测量，公开列出的补充表尚未提供所需位点量值。人UC论文已确认，处理矩阵登录号及访问仍待核对。
当前未确认上述材料中的CGNL1目标位点覆盖，因此没有对C04补造检验，也不把访问失败表述为没有生物学信号。总蛋白或RNA-only资料不替代磷酸化证据。

六、持续推进规则
下一步优先解决C04材料，并对C01/C02/C03的定位、处理来源和细胞组成竞争解释做有针对性的核查。保留原始阴性和所有处理敏感性，不以最小P选择主线。新的独立验证或功能实验材料若缺失，要记录具体缺口；不能靠重复计算提升证据等级。
bootstrap区间估计均值，精确符号秩P检验排序位置假设；二者不是同一个估计/检验。个别bootstrap区间不跨零不能替代预设P值判定。所有P与各自家族FDR完整保存在wave1/RESULTS.json和wave2/RESULTS.json。
'''
rf=O/'CGNL1_CONTINUOUS_PROGRESS_20261023_v1.txt';rf.write_text(report,encoding='utf-8');assert '\ufffd' not in rf.read_text(encoding='utf-8')
save('PROGRESS_ACCEPTANCE.json',dict(utc=plan['utc'],goal_status='active',wave1_QA=True,wave2_QA=True,
 valid_primary_tests_or_rechecks=14,unique_mainline_selected=False,C04_inferential_tests=0,
 report_sha256=sha(rf),goal_plan_sha256=sha(O/'GOAL_PLAN.json'),whole_project_complete=False))
checkpoint=dict(utc=plan['utc'],stage='CGNL1_four_candidate_continuous_exploration_active',workspace=str(O),report=str(rf),goal_status='active',whole_project_complete=False)
cp=P/'metadata/PROJECT_CHECKPOINT_CONTINUOUS_20261023.json';cp.write_text(json.dumps(checkpoint,indent=2)+'\n',encoding='utf-8')
statep=P.parent/'CGNL1_CURRENT_STATE_20261017.json';backup=O/'CURRENT_STATE_BEFORE_CONTINUOUS.json'
if not backup.exists():
 oldhash=sha(statep);s=read(statep);shutil.copy2(statep,backup);assert sha(backup)==oldhash
 s.update(utc=plan['utc'],stage=checkpoint['stage'],latest_checkpoint=str(cp),current_goal_workspace=str(O),whole_project_complete=False)
 assert sha(statep)==oldhash
 temp=statep.with_suffix('.json.continuoustmp');temp.write_text(json.dumps(s,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');temp.replace(statep)
print(json.dumps(dict(report=str(rf),sha256=sha(rf),goal_status='active'),ensure_ascii=False))
