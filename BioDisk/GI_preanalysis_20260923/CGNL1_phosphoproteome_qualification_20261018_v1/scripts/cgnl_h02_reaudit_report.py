from pathlib import Path
import json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'H02_localization_reaudit_20261024'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==z
r=read(D/'RESULTS.json');a=read(D/'ADDITIONAL_SOURCE_CHECK.json');obs=read(D/'OBSERVATION_AUDIT.json')
assert len(obs)==68 and r['all_four_mono_positive_n']==4 and all(x['localization_prob'] is None for x in obs)
b=next(x for x in a['B03_target_annotations'] if x['Phosphosite'].endswith(':s486'));assert b['Peptide']=='APsLGAQSK'
text='''CGNL1候选测量形式复核：H02定位字段与多磷酸化信号（2026-10-24）
本轮没有增加疾病效应检验，没有更改H02原17例结果。

H02定位证据表述修正
原始Phospho (STY)Sites.txt确实包含按实验命名的Localization prob列，因此此前笼统写“没有逐样本定位字段”不准确。更准确的结论是：S283、S297两行各有232个逐实验定位列，但这些列均为空。原C03_H02所用17对患者、两位点四测量共68个值，均有正强度，但对应逐实验定位概率均为空。
全局Localization prob分别为S283=1、S297=0.99995。全局值不能替代每个样本的定位证据；空的实验字段也不能单独证明所有强度都由match-between-runs产生。仍需要evidence.txt/msms.txt等更底层记录；当前H02归档目录未找到这两个文件。

H02不是单一磷酸化形式的直接重复
MaxQuant磷酸化位点表的___1/___2/___3强度分量分别区分单、双、三磷酸化形式；依据Cox实验室官方教程：https://cox-labs.github.io/coxdocs/modifications.html 。这里的含义应与SILAC标记数所称multiplicity区分。
在原17对患者中，S283的34个组织测量均有___1正信号，未见___2/___3正信号。S297的34个组织测量中，29个含___2/___3正信号，只有13个含___1正信号；两类可重叠，不能相加当作样本数。
若只要求两位点在肿瘤和正常组织均有___1正值，则只剩4对患者（L006、L012、L014、L032）。这仅是覆盖资格计数，没有对这4例重新计算P值或替换原17例模型。
原Intensity与三个分量之和存在最大相对差约0.00004076；保留源数值精度差异，不将它们写为逐位完全相等。
所以H02原S283−S297分析是位点层面总强度的比较，其中S297混合不同磷酸化形式；不能把它写成对PDAC或H01中S297单磷酸化形式的严格重复。原P=0.108887的阴性结果继续保留，也不能借此挑选4例寻求阳性。

C02及其他来源的比较
本轮重新从B03原表读取S486与S283标注：S486为APsLGAQSK；S283为RQDsAGPVLDGAR/QDsAGPVLDGAR，均为单标记源肽段。H01的既有序列核查亦将S297标注为RSsSSSTTPTSANSLYR。没有发现与C01 S198/S202完全相同的双标记解释问题，但这些标注仍不能替代B03/H01的逐谱图定位概率。
C02因测量对象较明确仍适合优先补证；C03应保留分来源、分磷酸化形式解释，不做当前不具备资格的合并推断。H03插补和患者重叠限制未消除。

下一步与状态
优先寻找H02底层鉴定记录或原始归档目录，以解释为何实验定位列为空并定位多磷酸化肽段。对于仅4例的单磷酸化完整子集，即便后续计算也只能预先限定为探索性描述，不能取代原主分析。继续C04可访问材料搜索。所有原始统计保留；Goal保持active。
证据目录：continuous_exploration_20261023_v1/H02_localization_reaudit_20261024。
'''
f=O/'CGNL1_H02_MEASUREMENT_REAUDIT_20261024.txt';f.write_text(text,encoding='utf-8');assert f.read_text()==text
save(D/'ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_hash_verified=r['source_sha256'],observations=68,no_perexperiment_probability_available=True,global_probability_not_substituted=True,new_tests=0,original_H02_null_retained=True,report_sha256=sha(f),artifacts=[dict(path=str(D/n),sha256=sha(D/n)) for n in ['SPEC_BEFORE_COUNTS.json','RESULTS.json','OBSERVATION_AUDIT.json','ADDITIONAL_SOURCE_CHECK.json','SOURCE_TARGET_ROWS.json']]))
s=read(O/'STATUS.json');s.update(H02_localization='Columns exist butS283/S297all232experimentfields eachblank;68primarymeasurementprobs unavailable',H02_phosphoforms='S29729of34observationscontainmulti-phosphosignal;mono-onlycomplete4patients;original17patientnullretained',latest_measurement_report=str(f),next_task='H02 evidence/msms source inventory and C04 accessible data sources');save(O/'STATUS.json',s)
print(json.dumps(dict(report=str(f),sha256=sha(f))))
