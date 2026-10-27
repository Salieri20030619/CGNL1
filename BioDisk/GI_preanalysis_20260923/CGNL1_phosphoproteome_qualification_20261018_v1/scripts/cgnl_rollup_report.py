from pathlib import Path
import json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'material_search_restart_20261024';A=D/'rollup_audit'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==z
r=read(A/'RESULTS.json');assert all(x['max_value_difference']==0 and not x['missing_mismatch'] and x['reference_intensity_difference']==0 for x in r['comparisons'])
report='''CGNL1 C01测量对象澄清：S198行代表S198/S202双磷酸化形式
日期：2026-10-24。此为新增源表证据审计，不改写旧模型或其统计值。

已确认的结果
本轮从官方UMich包取得PDAC肽段表（133238028字节）和多位点表（214613421字节），ZIP成员CRC及本地/服务器SHA256一致。三个层级的CGNL1记录分别为单位点22行、多位点33行、肽段18行。
S198单位点行与多位点表中S198S202行：275个导出通道的缺失模式完全一致，209个共同数值逐一相等，参考强度也相等。源肽段PSNsQPTsPSLEDPAK包含S198/S202两个标记。多位点表中只有这一条明确列出S198，没有找到孤立S198形式的导出行。
S202单位点行与S202单磷酸化行：275通道缺失一致，220个共同数值逐一相等，源肽段PSNSQPTsPSLEDPAK。
S297单位点行与S297单磷酸化行：275通道缺失一致，209个共同数值逐一相等，源肽段RSsSSSTTPTSANSLYR;SsSSSTTPTSANSLYR。
275是含QC等列的导出通道数，不是275例独立患者；原C01主对比仍是28对独立患者。

对C01的修正
旧标签“S198−S297”应在当前源材料范围内解释为“S198/S202双磷酸化形式相对于S297单磷酸化形式的变化”。原28对患者效应、P值和FDR保留，不借此更换阳性结果。可以继续研究CGNL1磷酸化状态组合，但不能把该结果作为孤立S198位点的功能证据，不能据此计算S198占有率或只设计S198单突变就宣称解释了完整现象。
没有检出孤立S198导出行，不等于其在生物学上不存在。当前证据是汇总表对应，不是逐谱图定位概率，也不是新的独立队列。

算法依据与边界
TMT-Integrator方法论文描述：单位点报告优先采用单磷酸化形式；只在多磷酸化形式中出现的位点也可以进入单位点汇总。该说明与本次表间对应相容。引用：https://www.nature.com/articles/s41467-026-70118-7 （2026年方法说明，不能自动证明历史软件版本全部行为）。
官方2.4.0标签指向的公开树只含README等文件，没有可核查的Java实现；未假称已复现历史源码。我们依靠当前队列的实际源表数值建立上述对应。

下一步问题
一个更贴近测量对象的后续候选是：在同一肽段序列PSNSQPTSPSLEDPAK上，S198/S202双磷酸化形式相对于S202单磷酸化形式，是否随PDAC肿瘤状态改变？两种形式的响应不能直接当作绝对占有率，但该对比比跨不同肽段直接归因为S198更具体。
这属于由新源证据提出的探索性方案变更。必须先冻结新对比、样本纳入与检验家族，再查看其结果；不替换原主对比，也不把选中的新对比包装成预先独立验证。本轮未计算此新对比的P值或选择其结果。
C02与C03既有分析保留；C04材料问题仍未解决。Goal保持active，下一轮先冻结上述同序列磷酸化形式对比并完成适用分析及独立QA，同时继续材料来源资格检查。
'''
f=D/'CGNL1_PHOSPHOFORM_CLARIFICATION_20261024.txt';f.write_text(report,encoding='utf-8');assert f.read_text()==report
save(A/'ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),three_correspondences_exact=True,source_counts=r['source_row_counts'],new_tests=0,report=str(f),report_sha256=sha(f),artifacts=[dict(path=str(A/n),sha256=sha(A/n)) for n in ['SPEC_BEFORE_COMPARISON.json','RESULTS.json','TARGET_SOURCE_ROWS.json']]))
s=read(O/'STATUS.json');s.update(C01='S198/S202doublephospho versusS297monophospho sourceidentity confirmed;isolatedS198 interpretation unsupported',latest_phosphoform_report=str(f),next_task='Freeze exploratory same-sequence S198S202-versus-S202 phosphoform contrast;retainoriginalC01andnewtestfamily separately');save(O/'STATUS.json',s)
c=read(O/'CANDIDATE_COMPARISON_v1.json');c['utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
for x in c['candidates']:
 if x['id']=='C01':x.update(question='PDAC CGNL1 phosphoform remodeling: S198/S202 double versus S297 single;isolatedS198 not established',new_source_evidence=str(A/'RESULTS.json'))
save(O/'CANDIDATE_COMPARISON_v2.json',c)
print(json.dumps(dict(report=str(f),sha256=sha(f))))
