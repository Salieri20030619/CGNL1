from pathlib import Path
import json,hashlib,datetime,copy
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1')
old=json.loads((O/'CANDIDATE_COMPARISON_v4.json').read_text())
for item in old['evidence_manifest']:assert hashlib.sha256((O/item['path']).read_bytes()).hexdigest()==item['sha256']
F=O/'C04_full_fraction_qualification';D=F/'disease_descriptive_20261026';R=F/'relative_peptide_pattern_20261026'
assert json.loads((D/'ACCEPTANCE.json').read_text())['pass_all'] and json.loads((R/'INDEPENDENT_R_QA.json').read_text())['pass']
ratio=json.loads((R/'RESULTS.json').read_text());effects=json.loads((D/'RESULTS.json').read_text())['results'];by={(x['endpoint'],x['model']):x for x in effects}
v=copy.deepcopy(old);v['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();v['supersedes']='CANDIDATE_COMPARISON_v4.json'
c=next(x for x in v['candidates'] if x['id']=='C04');c.update(priority='retain relative phosphopeptide pattern candidate; do not prioritize a large robust S283 decline',strength='13-file identification scope complete; mouse S284 measurement supports human S283 correspondence; complete eight-animal measurement; independent numerical QA',gap='S284 effect attenuates after normalization and primary LOO crosses zero; stronger peptide signal has unresolved sites; total protein and composition unresolved; n8; no functional evidence',next='finish protein pilot and qualification; preserve peptide-level pattern and develop localization-resolved targeted validation',new_result=dict(S284_N1_log2_difference=by[('D01','N1')]['mean_log2_difference'],S284_N1_LOO=[by[('D01','N1')]['LOO_min'],by[('D01','N1')]['LOO_max']],unresolved_peptide_N1_log2_difference=by[('D02','N1')]['mean_log2_difference'],posthoc_relative_peptide_log2_difference=ratio['mean_log2_ratio_difference'],posthoc_not_independent_validation=True,no_P_values=True))
new_evidence=[D/'ACCEPTANCE.json',D/'RESULTS.json',R/'SPEC_POSTHOC.json',R/'RESULTS.json',R/'INDEPENDENT_R_QA.json',F/'joint_spectral_qualification_20261026/HUMAN_MOUSE_PEPTIDE_MAPPING.json']
for p in new_evidence:v['evidence_manifest'].append(dict(path=str(p.relative_to(O)),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
v['working_question']='Does the relative composition of CGNL1 phosphorylation forms change with digestive disease context, and can a reliably localized pattern be linked to junction or microtubule organization?'
v['working_question_is_not_validated_mainline']=True;v['final_unique_mainline_selected']=False;v['whole_project_complete']=False
report=O/'CGNL1_CANDIDATE_REVIEW_20261026_v5.txt';v['report']=str(report)
lines=['CGNL1 候选方向比较 v5（2026-10-26）','',
'当前判断：保留多磷酸化形式和疾病环境差异的问题，不把主线锁定为 S283 在所有消化道疾病中同向下降。C04 新结果没有支持稳健的大幅 S283 对应位点下降，却提供了磷酸肽相对模式变化的探索性线索。尚未确定唯一最终主线。','',
'候选问题：CGNL1 不同磷酸化形式的相对构成是否随消化道疾病状态改变；在可靠定位后，这种变化是否关联上皮连接或微管组织功能？这是一项待验证研究问题，不是已证实机制。','',
'C01：PDAC 磷酸化形式重排','保留 S198/S202 双磷酸肽相对 S297 单磷酸肽的配对模式。原“单独 S198”解释已降级；同一肽段双磷酸与单磷酸的比较不支持原来的简单解释。优先缺口是同疾病独立队列和逐观测定位证据，不能把同队列不同处理当独立验证。','',
'C02：S486 与疾病环境相关的变化','测量资格较好，B03 与 Mayo 均有可用位点证据；但 Mayo 的小样本方向未重复 B03 总体方向，其中一例对参照和校正较敏感。保留环境依赖问题，不称普适 iCCA 标志物；需要独立、可比队列及功能验证。','',
'C03：HCC 位点模式异质性','保留跨队列模式比较和动态调控候选。H02 阴性及 S297 直接定位支持不足、H03 填补处理、H01/H02 受试者独立性未核清，都限制合并推断。已发表扰动线索不能代替所选位点的功能证明。','',
'C04：非肿瘤炎症中的磷酸肽相对模式','13 个原始文件及联合鉴定已完成，5 对照/3 炎症小鼠。小鼠 S284 在序列上对应人类 S283，但人鼠肽段有一处 I/V 差异，不是完全相同肽段。',
f"S284 主归一化模型平均 log2 差 {by[('D01','N1')]['mean_log2_difference']:.6f}；逐动物剔除范围 [{by[('D01','N1')]['LOO_min']:.6f}, {by[('D01','N1')]['LOO_max']:.6f}]，方向不稳。定位未确定磷酸肽的主模型 log2 差 {by[('D02','N1')]['mean_log2_difference']:.6f}，3/3 炎症动物低于对照中位数。",
f"看到上述结果后追加的事后比值诊断：未定位肽段/S284 肽段的组间 log2 比值差 {ratio['mean_log2_ratio_difference']:.6f}，相对比值为 {ratio['relative_ratio_fold_change']:.6f} 倍；逐动物剔除范围 [{ratio['LOO_min']:.6f}, {ratio['LOO_max']:.6f}]。同通道归一化因子在比值中代数抵消，不能把 N0/N1/N2 一致当作三次独立验证。",
'该比值是不同磷酸肽的相对模式，不是占有率，不能单凭比值证明总蛋白无关或排除细胞/异构体组成和肽段测量差异。总蛋白双批次读取试点仍待完成及验收；未扩展全范围，也未完成蛋白校正。','',
'共同决策依据与停止条件','按测量可靠性、独立重复、疾病解释、创新与功能可行性判断，不按最小 P 值排序。C04 本轮均为描述性探索，不计算 P/FDR/置信区间。阴性、方向冲突、敏感性和不可评估结果保留在对应完整报告。','当前较值得保留的是可测量的相对磷酸化模式及其环境差异；不能因为某些位点机制未报道就宣称具有创新性。须在定位可靠、同疾病独立重复及功能实验能够区分假说之后确定机制主线。','功能计划已提出，但尚未执行；基因级 WT 救援、位点状态和定位/表达控制、连接或微管功能读出及独立疾病验证均不能由本轮数值 QA 代替。','',
'验收说明','C01-C03 既有证据文件哈希与 v4 清单一致，本次未重跑或改写统计。C04 新分组比较及事后比值均经独立 R 从原始输入复算，数值误差为 0。v4 和原始分析保留。','本比较为阶段性证据更新，不是全项目完成或机制确认。']
assert not report.exists();report.write_text('\n'.join(lines)+'\n',encoding='utf-8');(O/'CANDIDATE_COMPARISON_v5.json').write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
receipt=dict(utc=v['utc'],report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),comparison_sha256=hashlib.sha256((O/'CANDIDATE_COMPARISON_v5.json').read_bytes()).hexdigest(),unchanged_C01_C03_evidence_verified=True,posthoc_label_preserved=True,whole_project_complete=False)
(O/'CANDIDATE_REVIEW_v5_ACCEPTANCE.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
sfile=O/'STATUS.json';s=json.loads(sfile.read_text());s['latest_candidate_comparison']=str(O/'CANDIDATE_COMPARISON_v5.json');s['C04_relative_pattern']=str(R/'RESULTS.json');s['latest_candidate_review']=str(report);tmp=sfile.with_suffix('.tmp');tmp.write_text(json.dumps(s,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');tmp.replace(sfile);print(json.dumps(receipt))
