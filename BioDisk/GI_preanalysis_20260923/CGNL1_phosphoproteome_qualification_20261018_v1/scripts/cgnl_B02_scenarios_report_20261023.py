from pathlib import Path
import json,hashlib,datetime,shutil
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
O=P/'results/B02_exploratory_scenarios_20261023_v1'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
r=read(O/'RESULTS.json');qa=read(O/'INDEPENDENT_NUMERICAL_QA.json');a=read(O/'ACCEPTANCE.json');paper=read(O/'PAPER_DATA_COMPARISON.json')
assert qa['pass_QA'] and all(x['fit_valid'] for x in r)
assert all(sha(Path(x['path']))==x['sha256'] for x in a['artifacts'])
def write_json(name,d):
 f=O/name;f.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==d
labels={'S0':'214例，无批次项','S1':'209例，排除5例重复患者，无批次项','S2':'214例，重复患者归原批次，加入批次固定效应','S3':'214例，重复患者归第22批，加入批次固定效应'}
lines=['情景\t独立患者数\tS283 OR\t95%CI\t原始P\t四情景BH FDR']
for x in r:lines.append(f"{x['scenario']} {labels[x['scenario']]}\t{x['n']}\t{x['OR']:.6f}\t{x['CI95_low']:.6f}–{x['CI95_high']:.6f}\t{x['P']:.9f}\t{x['FDR_4_scenarios']:.9f}")
discrepancies=paper['other_batch_label_discrepancies'];cases=sorted(set(x['case'] for x in discrepancies))
write_json('EXPLORATORY_DECISION.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
 selected_source_label_fidelity_scenario='S2',selection_basis='214例源值不变，5例重复者最终列的原论文分组标签与原批次一致。结果计算前已记录此依据。',
 selection_not_based_on_P_or_effect_direction=True,
 same_endpoint_published_effect_available=False,matching_published_effect_not_claimed=True,
 conclusion='当前数据不支持S283在控制总CGNL1蛋白、年龄、性别、肿瘤大小后的血管侵犯关联。四情景方向相同但置信区间宽，不证明效应为零或四种模型等价。',
 nominal_significant_scenarios=0,all_four_scenarios_retained=True,
 label_issue=dict(id='CGNL-A10',severity='medium',status='source_header_discrepancy_documented',cases=cases,
  evidence='原Table S1相邻批次块重复出现Group06标签；对应NODE通道设计为iCCA_07。模型按NODE原始通道表编码，这一差异不静默抹除。',
  impact='S2代表原始通道设计加作者重复患者显示标签的一致性情景，不声称逐字复现所有原文批次表头。'),
 remaining_source_gaps=['重复值最终聚合或选择规则','精确重复制备层级','未取得原始通道级量值'],
 interpretation='新的探索性关联/敏感性分析；不是作者已发表CGNL1特异结论的成功复现，也不是机制验证',whole_project_complete=False))
table='\n'.join(lines)
report=f'''B02四情景探索性建模与原论文数据对照（2026-10-23）

结果：四种重复处理情景均未出现S283与血管侵犯的名义显著独立关联。OR为0.728–0.759，全部95%CI跨1。在本次预设范围内，更换5例重复患者的处理假设没有改变结论。置信区间较宽，不能据此证明不存在效应。

一、与原论文的对应关系
原文为Dong等，Cancer Cell 2022，DOI 10.1016/j.ccell.2021.12.006。原始Table S1含214例最终患者。逐例核对S283与总CGNL1两行，共428个数值与NODE作者最终表一致（容差1e-12）。
已核查主文、补充图说明及Table S1，没有同口径的CGNL1/S283—血管侵犯效应值作为复现参照。因此，本次复现的是作者数据输入；关联模型属于新开展的探索性分析，不宣称复现作者已有CGNL1结论。本轮没有把其他表型或癌种的方向用作匹配目标。

二、模型与四种情景
独立单位为患者。血管侵犯Yes=1、No=0。214例含97例Yes、117例No；排除5例后为209例，95例Yes、114例No。
所有模型使用同一S283连续变量，并控制总CGNL1蛋白、年龄（每10岁）、性别（男性相对女性）、肿瘤最大径（线性，每cm）。作者log2相对丰度不再归一化、不二分、不插补。S283 OR表示在其他变量相同时，位点相对丰度每增加1个作者log2单位的条件比值比，不表示绝对磷酸化占有率。
S2/S3加入22批次固定效应。不存在可用的分时点或分位置重复量值，因此没有虚构纵向或空间重复模型。

{table}

三、选择与原数据最一致的情景
保留S2作为“作者标签一致的探索性情景”：保留214例，并将5例重复患者归入原批次，与原始Table S1中这些最终患者列的分组标签一致。这项选择依据在拟合前已固定，不根据P值、效应大小或有利方向择优。S0保留为简约内部基准，S1、S3完整报告。
原Table S1另有批次表头差异：{len(cases)}例所在块标为Group06，但原始NODE设计为iCCA_07。此问题记为CGNL-A10；模型依照已验证的NODE通道设计编码。S2不等于已还原作者未知的重复数值合并过程。

四、模型诊断和验收
四种模型均收敛、设计矩阵满秩，无完全分离警告。Python statsmodels与独立R glm的目标系数、标准误及P值差异均小于1e-6。原始数据、此前结果及计划文件哈希未变。
无批次模型5个非截距参数；有批次模型26个非截距参数，每参数仅约3.73个较少结局类别观测。批次模型的参数容量有限，CI更宽，仅作探索性敏感性结果。数值收敛与独立计算一致不等于模型已获充分生物学验证。
四情景针对同一假设、使用高度重叠患者，并非四次独立验证。原始P全部报告；四次检验的BH FDR作为补充，不用于掩盖名义信号。本轮不存在原始P<0.05的情景。

五、结论边界及后续
当前数据不支持所设定的S283—血管侵犯条件关联。没有通过替换重复处理方案获得阳性结果。不能将同向OR小于1写成保护作用，也不能认为四种处理方法已被证明等价。重复值合并规则和制备层级仍未解决，原始通道量值仍不可用。
本轮四情景任务已完成，不重复更换结局、切点或协变量寻找阳性；整体CGNL机制主线仍未完成。

服务器完整结果目录：{O}
关键文件：SPECIFICATION_BEFORE_FITS.json、PAPER_DATA_COMPARISON.json、SCENARIO_RESULTS.csv、ALL_COEFFICIENTS.csv、INFLUENCE.csv、R_INDEPENDENT_QA.csv、EXPLORATORY_DECISION.json及ACCEPTANCE.json。
'''
f=O/'B02_EXPLORATORY_SCENARIOS_REPORT_20261023.txt';f.write_text(report,encoding='utf-8');assert '\ufffd' not in f.read_text(encoding='utf-8')
checkpoint=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage='B02_four_exploratory_scenarios_complete_no_nominal_S283_signal',
 result_root=str(O),report=str(f),numerical_QA_pass=True,selected_source_fidelity_scenario='S2',published_CGNL1_conclusion_reproduced=False,
 all_four_results_retained=True,whole_project_complete=False)
cf=P/'metadata/PROJECT_CHECKPOINT_B02_SCENARIOS_20261023.json';cf.write_text(json.dumps(checkpoint,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
statepath=P.parent/'CGNL1_CURRENT_STATE_20261017.json';backup=O/'CURRENT_STATE_BEFORE_B02_SCENARIOS.json'
if not backup.exists():
 oldhash=sha(statepath);state=read(statepath);shutil.copy2(statepath,backup);assert sha(backup)==oldhash
 state.update(utc=checkpoint['utc'],stage=checkpoint['stage'],latest_checkpoint=str(cf),latest_exploratory_B02_scenarios=str(O),whole_project_complete=False)
 assert sha(statepath)==oldhash
 temp=statepath.with_suffix('.json.B02tmp');temp.write_text(json.dumps(state,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');temp.replace(statepath)
write_json('DELIVERY_ACCEPTANCE.json',dict(utc=checkpoint['utc'],report_sha256=sha(f),models=4,QA_pass=True,
 result_acceptance_sha256=sha(O/'ACCEPTANCE.json'),decision_sha256=sha(O/'EXPLORATORY_DECISION.json'),
 exploratory_only=True,whole_project_complete=False))
print(json.dumps(dict(report=str(f),sha256=sha(f),selected='S2',header_discrepancy_patients=len(cases)),ensure_ascii=False))
