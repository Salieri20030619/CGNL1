from pathlib import Path
import json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'HCC_source_independence_20261024';D.mkdir(exist_ok=True)
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==z
assert 'PMC10783603 TableS2 mmc3 and TableS5 mmc6'==read(P/'metadata/HCC_EXPLORATORY_SPEC571.json')['source']
source=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576/H01/S049_HCC_Clinical_Information_and_TMT11_Sample_Mapping_Gao2019.r1.xlsx');assert source.exists()
z=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),id='C03_SOURCE_INDEPENDENCE_01',status='partially_resolved',
 source_trace=[dict(cohort='H01',label='Gao2019',original_doi='10.1016/j.cell.2019.08.052',original_pmid='31585088',mapping_file=str(source),mapping_sha256=sha(source)),dict(cohort='H02',label='Jiang2019',doi='10.1038/s41586-019-0987-8',archive='PXD008373'),dict(cohort='H03',pmcid='PMC10783603',analysis_source='OwnTableS2 andTableS5',own_recruitment_center='Mengchao Hepatobiliary Hospital of Fujian Medical University',proteome_recruited=160,phosphoproteome_reported=132)],
 primary_evidence_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC10783603/',primary_evidence_sections=['Results:Cross-validation with Gao andJiangcohorts','Methods:Patients and tissue samples','SupplementaryTableS5 caption'],
 established='H03is presented as newlycollected owncohort,with Gao/Jiang as externalcomparators;ourH03inputsare its ownabundance tables,not their externalvalidation tables',
 does_not_establish=['Patient-by-patient identifierdeduplication across all3','H01vsH02 donor independence','Noimputation/no-localization-bias inH03','Three reliable independent confirmations ofCGNL1siteeffect'],
 bibliographic_note='H03reference25 points toGao correctionDOI10.1016/j.cell.2019.10.038 PMID31730861;originalpaper DOI2019.08.052 PMID31585088;do not count correction asanothercohort',new_tests=0)
save(D/'SOURCE_INDEPENDENCE_AUDIT.json',z)
text='''CGNL1 HCC来源独立性补充核查（2026-10-24）
结论：H03作为自身采集队列的来源依据已补强；此前笼统的“所有来源重叠未知”应细化，而不是直接升级为三次可靠独立验证。

本项目源文件追溯
H01临床/TMT映射文件明确标注Gao2019，对应原论文DOI 10.1016/j.cell.2019.08.052（PMID31585088）。H02来源为Jiang2019、DOI 10.1038/s41586-019-0987-8及PXD008373。H03原冻结方案明确使用PMC10783603的TableS2/mmc3与TableS5/mmc6。

新核实的来源证据
H03原文将自身队列与Gao、Jiang的外部队列分列比较，Methods说明自身患者来自福建医科大学孟超肝胆医院；TableS5是其132例患者的磷酸化数据。因此，本项目H03不是仅将Gao/Jiang的外部比较表重新命名为新来源。
来源：https://pmc.ncbi.nlm.nih.gov/articles/PMC10783603/
原文关于三队列来源差异的描述提供研究设计层面的支持，但不能替代逐患者ID去重；H01与H02之间的个体重叠仍未在本项目直接核实。作者单位、发表年份或人数不同本身均不能作为无重叠的充分证明。

对研究结论的影响
可以更明确地把H03视为作者报告的另一采集队列，保留其对跨来源方向比较的价值；不能据此解除H03的GRR插补/原始观测掩码问题，也不能忽略H02原17例阴性、其S297测量形式混合及仅6/34组织测量有当前标准下直接定位支持。
原文队列人数不是每个CGNL1对比的有效样本量；本项目S283−S297仍按各来源实际完整患者数分别报告，不合并或重新计算P值。
另核实H03参考文献中的Gao DOI 10.1016/j.cell.2019.10.038对应更正条目（PMID31730861）；原论文是上述2019.08.052。更正条目不是另一队列。

作者代码范围核查
已只读取得并校验v1.0归档的4个文本文件，未执行作者代码。代码和README将本队列、Jiang及Gao输入分开，用于蛋白亚型模型；明确的impute.RegImpute调用属于Jiang蛋白数据处理。未在这些文本中找到H03磷酸化观测掩码。大型RData内部未检查，因此不据文件名宣称其中完全没有相关材料。

状态
该问题标记为“部分解决”：H03自身来源与外部比较来源已区分；H01/H02个体独立性及测量可比性仍需具体证据。本轮没有新模型或新的疾病效应检验。C04全分馏下载继续在原生会话中运行，尚未完成整批验收；Goal保持active。
'''
f=O/'CGNL1_HCC_SOURCE_INDEPENDENCE_20261024_v2.txt';f.write_text(text,encoding='utf-8');assert f.read_text()==text
save(D/'ACCEPTANCE.json',dict(utc=z['utc'],source_paths_verified=True,new_tests=0,not_patientlevel_dedup_claim=True,report_sha256=sha(f)))
s=read(O/'STATUS.json');s.update(HCC_source_independence='H03owncohort vsGao/Jiangexternalcomparators clarified;H01-H02individualoverlapunverified;measurementlimitsretained',latest_HCC_source_report=str(f));save(O/'STATUS.json',s)
print(json.dumps(dict(report=str(f),sha256=sha(f))))
