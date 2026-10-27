from pathlib import Path
import json,csv,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';W=O/'wave4_phosphoform'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==z
r=read(W/'RESULTS.json')
with (W/'independent_R_QA.csv').open() as h:q=next(csv.DictReader(h))
errors={k:abs(float(q[k])-r[k]) for k in ['n','mean','median','P','alternative_P']};assert max(errors.values())<1e-12
with (W/'patient_contrasts.csv').open() as h:rows=list(csv.DictReader(h))
assert len(rows)==len(set(x['case'] for x in rows))==29
for x in rows:
 assert abs(float(x['double_T'])-float(x['double_N'])-float(x['mono_T'])+float(x['mono_N'])-float(x['value']))<1e-12
spec=read(W/'SPEC_BEFORE_RESULTS.json')
for p,h in spec['sources'].items():assert sha(Path(p))==h
text=f'''CGNL1同序列磷酸化形式探索性比较（2026-10-24）
主要结论：当前数据不支持PDAC中S198/S202双磷酸化形式相对于S202单磷酸化形式发生配对差异。不能把此前跨肽段C01结果归因于孤立S198的独立变化。

问题与方法
基于新源表发现，比较同一肽段PSNSQPTSPSLEDPAK的双磷酸化S198/S202与单磷酸化S202形式。每个患者计算：(双磷酸化肿瘤−双磷酸化正常)−(单磷酸化肿瘤−单磷酸化正常)。UMich多位点导出为主要输入；BCM对同一批患者作处理敏感性核查。
本方案是在已知既有位点结果后提出的探索性细化，已在计算新对比前冻结。只纳入四个测量均完整且样本身份唯一对应的患者；无新增插补、无事后删除异常值。44对可映射患者中29对纳入，15对因至少一个测量缺失未进入此对比。未将技术通道当作生物学样本。
本轮仅一个预先限定的新对比，报告原始P；不重新计算原有检验家族FDR，不把本轮单项P称为校正后证据。检验为双侧精确符号秩；均值区间采用10000次患者bootstrap，seed20261024。R独立复算均值、中位数、P与敏感性P，差异小于1e-12。

结果
n={r['n']}；均值={r['mean']:.6f}；中位数={r['median']:.6f}；17例向上、12例向下。
原始P={r['P']:.12g}。95% bootstrap均值区间=[{r['bootstrap_mean_CI95'][0]:.6f},{r['bootstrap_mean_CI95'][1]:.6f}]。
BCM相同29例均值={r['alternative_BCM_mean']:.6f}，P={r['alternative_P']:.12g}，结论未改变。同队列处理比较不构成独立验证。
留一均值范围=[{r['LOO_mean_range'][0]:.6f},{r['LOO_mean_range'][1]:.6f}]；没有据此删除患者或挑选子集。

对候选方向的影响
1. C01原S198−S297的28对患者结果保留；其测量实质为双磷酸化S198/S202形式与单磷酸化S297形式的比较。本轮结果不否定该跨肽段相对模式，但不能将它缩窄为“S198特异调控”。
2. 本轮“不支持差异”不是两种形式等效或S198没有功能的证明；未进行等效检验，覆盖也不同（29对与原28对）。不能把两个不同样本集的均值直接相减来推导新的因果分量。
3. 后续资源排序不应再以孤立S198为优先机制靶点。C02 S486仍值得优先补齐定位与独立重复；C01/C03保留为磷酸化形式和S297相关相对模式问题，是否能解释为S297调控仍需要新证据。此为基于测量明确程度的判断，不是以本轮P值决定某条路线有效或无效。
4. 当前结果均为相对磷酸化形式丰度变化，不能换算为绝对位点占有率。离子响应、肽段汇总、逐谱图定位和组织组成等限制仍在。C04未产生新数据。

交付和状态
冻结方案、逐患者输入、15例排除原因、Python结果、R复算与验收保存在wave4_phosphoform。原分析和全部关键阴性保留。Goal保持active；下一步优先核对C02及其他候选位点是否也存在类似多位点汇总问题，并继续C04材料资格工作。
'''
f=O/'CGNL1_PHOSPHOFORM_RESULT_20261024.txt';f.write_text(text,encoding='utf-8');assert f.read_text()==text
save(W/'ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),independent_R_QA_pass=True,max_QA_error=max(errors.values()),source_hashes_unchanged=True,unique_patients=29,new_tests=1,primary_null_preserved=True,report_sha256=sha(f),artifacts=[dict(path=str(W/n),sha256=sha(W/n)) for n in ['SPEC_BEFORE_RESULTS.json','RESULTS.json','EXCLUSIONS.json','patient_contrasts.csv','independent_QA.R','independent_R_QA.csv','R_session.txt']]))
s=read(O/'STATUS.json');s.update(wave4='phosphoformcontrast n29 mean0.174445 P0.213245 independentRQApass;isolatedS198claimunsupported',latest_phosphoform_result=str(f),next_task='Candidate phosphoform interpretation audits beyondC01;C04 sourcequalification');save(O/'STATUS.json',s)
c=read(O/'CANDIDATE_COMPARISON_v2.json');c['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
for x in c['candidates']:
 if x['id']=='C01':x.update(refinement='Samepeptidedouble-vsmono n29P0.213245;retainoriginalcrosspeptidepattern,do notprioritizeisolatedS198mechanism',refinement_evidence=str(W/'RESULTS.json'))
save(O/'CANDIDATE_COMPARISON_v3.json',c)
print(json.dumps(dict(report=str(f),sha256=sha(f),QAerrors=errors)))
