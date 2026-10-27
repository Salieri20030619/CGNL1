import pathlib,json,pandas as pd,numpy as np,csv,hashlib,gzip,zipfile
from scipy import stats
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');R=B/'GC_validation_v2';O=R/'04_results';D=R/'05_reports'
D.mkdir(exist_ok=True)
# Plain markdown tables without optional renderer dependencies.
def md(df):
 d=df.copy()
 for col in d:
  if pd.api.types.is_float_dtype(d[col]):d[col]=d[col].map(lambda x:f'{x:.6g}' if pd.notna(x) else 'NA')
 return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(str(x).replace('|','/') for x in row)+' |' for row in d.values)
main=pd.read_csv(O/'bulk_robustness'/'primary_effect_CI.tsv',sep='\t');lg=pd.read_csv(O/'bulk_robustness'/'leave_one_gene_out.tsv',sep='\t');lp=pd.read_csv(O/'bulk_robustness'/'leave_one_patient_out.tsv',sep='\t');co=pd.read_csv(O/'bulk_robustness'/'composition_proxy_sensitivity.tsv',sep='\t');cr=pd.read_csv(O/'bulk_robustness'/'paired_change_correlations.tsv',sep='\t');ps=pd.read_csv(O/'paired_scRNA'/'paired_gene_summary.tsv',sep='\t');ex=pd.read_csv(O/'independent_scRNA_refined'/'prespecified_endpoint_validation.tsv',sep='\t');corr=pd.read_csv(O/'independent_scRNA_refined'/'exploratory_state_correlations.tsv',sep='\t');summary=json.loads((O/'independent_scRNA_refined'/'summary.json').read_text());gate=pd.read_csv(O/'independent_scRNA'/'gating_reference_refined_patient_QC.tsv',sep='\t');initial=pd.read_csv(O/'independent_scRNA'/'gating_reference_patient_QC.tsv',sep='\t')
text='''# 胃癌第二轮验证内部证据附件

日期：2026-09-23。目录：sgs10:/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2。

## 决策

方向1暂缓，本轮未新增CRC分析。方向2支持启动有验收门槛的正式发现阶段，暂题为胃癌腹膜转移脂质相关髓系微环境的多组学解析与功能验证。当前不是新机制确认，也没有最终靶点或可用临床预测工具。

## 分析单位和预设范围

组织分析沿用17位独立患者、34份组织和固定8基因脂质摄取面板，患者与FFPE/冻存方式纳入模型。单患者排除只检查稳定性；单基因排除、组成代理敏感性和相关性分别在对应家族内BH校正。脂肪标记为ADIPOQ、PLIN1、CIDEC、CIDEA、LEP，与脂质面板不重叠，但仍为表达代理而非实测细胞比例。

GSE234129按作者标签、原QC和每样本每细胞类至少20细胞进行患者级聚合；4对原发灶/腹水巨噬细胞只报告方向、效应和一致性。同一大类中的变化可涉及亚群比例，尚未证明单细胞内诱导。

独立GSE228598共28位患者，所有84个归档成员均成功读取。139314个输入细胞，经UMI≥500、基因数200至7000、线粒体比例≤25%过滤后保留111681个细胞。样本混合恶性腹水与腹腔冲洗液；未把样本类型当作有无转移分组。

最终巨噬细胞富集规则要求至少2个髓系核心基因及至少2个C1QA/C1QB/C1QC表达，并排除T/B或上皮混合谱系；CD36、TREM2、FABP4不参与分群。每位患者的细胞富集群与其淋巴细胞群进行聚合比较。3个定位终点在外部表达读取前固定。其他标记检验另列为探索性家族，参与分群的基因不进入验证性检验。

## 组织证据完整结果

主模型与上一轮一致：差0.981819，95%CI 0.467144至1.496494，P=0.00101371；原8个面板家族FDR=0.00810971。不是新增独立队列。

### 单基因排除

'''+md(lg)+'''\n\n去除FABP4后的显著性减弱是关键结果，客户版予以保留。所有逐患者排除效应仍为正，但不能弥补基因组成依赖。

### 单患者排除

'''+md(lp)+'''\n\n### 组成代理敏感性

'''+md(co)+'''\n\n### 配对变化相关性

'''+md(cr)+'''\n\n脂肪表达代理可解释较多组织差异。过度校正与共线性也可能影响模型，因此这些分析定位问题来源，不能替代实测组成或因果判别。

## 配对单细胞完整结果

'''+md(ps[ps.celltype=='Macrophage'])+'''\n\n其中CD36、TREM2均4/4同向；FABP4为3对正向、1对零差。APOE、GPNMB、LIPA等不具有同样的一致性。未以稀疏细胞表达或细胞数获得推断性显著。

## 细胞富集规则的修订与限制

最初较宽的成熟髓系标记含CD163、MS4A7等，参考标签整体阳性预测值仅61.96%、召回率90.83%，提示混入单核细胞。该规则未作为最终客户证据。依据谱系特异性改用C1Q限制后，参考标签整体阳性预测值84.32%、召回率69.43%。修订依据参考标签的一致性和谱系特异性，不以外部终点显著性作为阈值选择标准；初始外部结果已查看，初始与修订结果分目录保留。

这是参考数据中的一致性，不是外部队列的金标准准确率。最终群体应称C1Q巨噬细胞富集群，不能宣称穷尽所有TAM或达到完整图谱注释精度。

'''+md(gate)+'''\n\n## 独立28例定位结果

最终规则获得43113个C1Q巨噬细胞富集细胞，28位患者均满足至少20个该群细胞。统计单位是28位患者。

'''+md(ex)+'''\n\nCD36/TREM2在28/28例该群中检出，患者中位细胞检出比例分别31.70%和66.61%；FABP4在27/28例中检出，但患者中位细胞检出比例仅1.17%。这些结果支持髓系定位，不支持把FABP4视为广泛巨噬细胞标志，也没有证明转移特异性或脂质摄取活性。

## 最终群体内的探索性关联

固定CD36/TREM2两个参照与APOE、LIPA、PLIN2、GPNMB、MRC1、SPP1六个伙伴，共12个相关性假设，独立成家族BH校正。以下结果保留全部阴性与门控敏感性，不作为客户版主要立项依据。

'''+md(corr)+'''\n\n初始宽规则下CD36–PLIN2的FDR为0.0558，最终C1Q规则下为0.0184；TREM2相关伙伴的显著性也随规则变化。效应是同队列中的探索性关联，不是新增外部重复，更不能当作调控方向。两版文件均保留。

## 查新与多组学现状

- 2023年腹水研究已描述APOE/TREM2相关TAM，不能将其存在当作新发现。
- Gut 2025已报道THBS2 CAF、C3–C3AR1与SPP1 TAM联系。
- JITC 2026已报道LysoPS–GPR34–ERK/c-Jun–AXL/CD36及胞葬相关免疫逃逸。
- 2026年脂质相关CAF研究已有FABP4/CD36/ABCA8功能证据；2022年研究提示癌细胞FABP4可抑制胃癌转移。
- 正式筛选应避开直接重复这些轴。CD36/TREM2只作定位与状态参照，最终候选尚未锁定。
- PXD064711表已于上一轮取得。本轮补查公开页面、PROXI与只读接口，A–E组和重复关系仍未核实；没有做未经授权的账户操作，也没有猜测分组检验。
- GSE162213为44株胃癌细胞系的H3K27ac ChIP相关数据，不能直接验证巨噬细胞调控，保留为资源审计，不强行拼接。
- 本轮完成的是转录组多分辨率证据。空间或组织定位、分组明确的蛋白数据和靶向脂质组属于下一阶段，未被写成已完成的分子多组学闭环。

## 下一阶段验收

应先获得可复核细胞状态与跨患者候选，再完成候选级查新和独立数据核查，最后收敛1条主轴及1条备选。满足来源明确、重复性、创新增量、可干预性后进入实验。若效应仅由脂肪组织取材比例解释而缺少可重复的髓系状态或功能依赖，则收缩对应机制题目。

客户版计划中20至30对组织和6至8个独立供者是资源规划，未作为已完成样本，也未冒称充分功效；具体规模须由先导变异和目标效应确定。

## 复核入口

00_docs含固定面板、门控和假设计划。04_results分别为bulk_robustness、paired_scRNA、independent_scRNA（初始规则）及independent_scRNA_refined（最终规则）。scripts和logs保留实现与完成记录。05_reports含客户建议书及本附件。
'''
(D/'internal_evidence_zh.md').write_text(text,encoding='utf-8')
checks=[]
def ck(n,ok):checks.append({'check':n,'pass':bool(ok)})
ck('17 independent tissue patients',len(lp)==17 and lp.omitted_patient.nunique()==17);ck('all patient leave-out directions positive',(lp.effect>0).all());ck('4 author-annotated macrophage pairs',(ps[(ps.celltype=='Macrophage')&ps.gene.isin(['CD36','TREM2'])].n_pairs==4).all());ck('external patient count 28',summary['patients']==28);ck('gating excludes endpoint genes',not set(['CD36','TREM2','FABP4'])&set(sum(json.loads((R/'00_docs'/'external_gating_plan_refined.json').read_text())['lineage'].values(),[])));ck('three endpoint BH matches',np.allclose(stats.false_discovery_control(ex.P.values,method='bh'),ex.FDR_BH_prespecified_3_endpoints.values));ck('external effects defined patient-wise',(ex.n_pairs==28).all());ck('no cell n substituted in primary tests',summary['macrophage_eligible_patients']==28);ck('initial failed gate retained',(O/'independent_scRNA'/'gating_reference_patient_QC.tsv').exists());ck('negative FABP4 exclusion retained',float(lg.loc[lg.omitted_gene=='FABP4','P.Pr(>|t|)'].iloc[0])>.05);ck('composition major attenuation retained',float(co.loc[co.model=='adipocyte','P.Pr(>|t|)'].iloc[0])>.05);ck('all main stages complete',all(p.exists() for p in [O/'bulk_robustness'/'COMPLETE.txt',O/'paired_scRNA'/'COMPLETE.json',O/'independent_scRNA_refined'/'COMPLETE.json']))
pd.DataFrame(checks).to_csv(R/'logs'/'validation_checks.tsv',sep='\t',index=False);assert all(c['pass'] for c in checks)
files=[]
for p in sorted(R.rglob('*')):
 if not p.is_file() or any(s in p.parts for s in ['__pycache__']) or p.name.endswith('.part') or p.suffix=='.zip':continue
 if p.parent==R/'02_raw' or p.parent==R/'01_metadata' or p.parent==R/'scripts':
  h=hashlib.sha256()
  with p.open('rb') as f:
   while True:
    b=f.read(8*1024*1024)
    if not b:break
    h.update(b)
  files.append({'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':h.hexdigest()})
pd.DataFrame(files).to_csv(R/'01_metadata'/'integrity_manifest.tsv',sep='\t',index=False);print(json.dumps({'checks_passed':len(checks),'report':str(D/'internal_evidence_zh.md'),'source_files_hashed':len(files)}))
