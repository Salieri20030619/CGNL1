
from pathlib import Path
import json,hashlib,zipfile
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(r'C:\Users\31937\CGNL1_review_consensus_retest_20261002');D=R/'results'
(R/'figures').mkdir(exist_ok=True)
m=pd.read_csv(D/'continuous_models.tsv',sep='\t');p=pd.read_csv(D/'patient_continuous_profiles.tsv',sep='\t');a=p[p.stage=='cluster_primary_n10'];mch=pd.read_csv(D/'matched_models.tsv',sep='\t');mp=pd.read_csv(D/'matched_patient_effects.tsv',sep='\t');de=pd.read_csv(D/'CGNL1_detection_and_depth_diagnostic.tsv',sep='\t')
coverage=de.groupby('dataset')[['cells','CGNL1_counts','CGNL1_detected_cells']].sum()
colors={'GSE245552':'#d97706','GSE285991':'#2563eb'}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
fig,ax=plt.subplots(2,2,figsize=(13,9),layout='constrained')
for ds,v in a.groupby('dataset'):
 ax[0,0].scatter(v.CGNL1,v.junction,label=ds+' (n='+str(len(v))+')',color=colors[ds],s=42)
 b=np.polyfit(v.CGNL1,v.junction,1);xx=np.linspace(v.CGNL1.min(),v.CGNL1.max(),80);ax[0,0].plot(xx,np.polyval(b,xx),color=colors[ds],alpha=.7)
ax[0,0].set(xlabel='CGNL1 log1p pseudobulk CP10K',ylabel='Four-gene junction score',title='A  Reviewed endothelial clusters: 18 patients')
ax[0,0].legend(fontsize=9);ax[0,0].set_title('A  Reviewed endothelial clusters: 18 patients\nStudy-adjusted beta = -0.113; P = 0.366',fontsize=11)
ss=['cluster_primary_n10','cluster_primary_n20','cluster_primary_n50'];labels=['At least 10 cells','At least 20 cells','At least 50 cells']
for j,model in enumerate(['study_adjusted','study_depth_adjusted']):
 v=m[(m.scope=='pooled')&(m.model==model)].set_index('stage').loc[ss];y=np.arange(3)+(j-.5)*.19;ax[0,1].errorbar(v.slope,y,xerr=np.vstack([v.slope-v.CI_low,v.CI_high-v.slope]),fmt='o',capsize=3,label=['Study adjusted','Study + depth adjusted'][j])
ax[0,1].axvline(0,color='grey',lw=1);ax[0,1].set(yticks=range(3),yticklabels=labels,xlabel='Unstandardized slope, 95% CI',title='B  Prespecified cell-count ladder (n = 18 / 16 / 8)');ax[0,1].invert_yaxis();ax[0,1].legend(fontsize=8,loc='center right',bbox_to_anchor=(1,.31))
wide=mp[mp.stage.isin(['unmatched_n10','matched_primary_n10_pairs3'])].pivot(index='donor',columns='stage',values='effect')
for donor,row in wide.iterrows():ax[1,0].plot([0,1],[row.unmatched_n10,row.matched_primary_n10_pairs3],'-o',color=colors[donor.split(':')[0]],alpha=.8)
ax[1,0].axhline(0,color='grey',lw=1);ax[1,0].set(xticks=[0,1],xticklabels=['Unmatched\nP = 0.000663','Depth matched\nP = 0.286883'],ylabel='CGNL1+ minus CGNL1- junction score',title='C  Same 7 patients and same 7 eligible lesions')
st=['unmatched_n10','matched_caliper1_n10_pairs3','matched_caliper05_n10_pairs3','matched_primary_n10_pairs3'];v=mch[mch.scope=='pooled'].set_index('stage').loc[st]
ax[1,1].errorbar(v.effect,range(4),xerr=np.vstack([v.effect-v.CI_low,v.CI_high-v.effect]),fmt='o',capsize=4,color='#166534');ax[1,1].axvline(0,color='grey',lw=1);ax[1,1].set(yticks=range(4),yticklabels=['Unmatched','Caliper 1.0','Caliper 0.5','Caliper 0.25'],xlabel='Mean patient difference, 95% CI',title='D  Signal attenuates at first depth-matching step');ax[1,1].invert_yaxis()
fig.suptitle('CGNL1: reviewer-consensus identity retest | frozen primary and sensitivity analyses',fontsize=14)
fig.savefig(R/'figures/reviewed_identity_retest.png');fig.savefig(R/'figures/reviewed_identity_retest.svg');plt.close(fig)
summary="""CGNL1｜采用审阅共识后的复测说明
2026-10-02；内部完整分析与审计材料

结论
已采信主要谱系审阅共识，并完成原计划 A05 连续关联、A06 检测深度匹配的复测。当前两个已审阅队列的主分析不支持稳定的 CGNL1–连接评分正向关联。预设小样本分支存在名义显著的正向探索性信号，保留为候选线索；不能替代主结果。该结论不等于证明 CGNL1 无生物学价值，也不等于认定数据集质量不合格。

审阅来源与使用范围
第一份匿名审阅披露接触系统历史摘要，因此登记为有上下文暴露的匿名审阅。DeepSeek 声明未接触原标签、预测或此前结论；独立性按该声明记录，其交付未附输入包哈希和独立冻结清单。两者主要谱系协调后 43/44 簇一致，不代表所有亚型、纯度和置信度一致。
本次采用：GSE245552 簇21（匿名 D02/C028），482细胞；GSE285991 簇14（D01/C018），1,618细胞。44簇共识已经回填到新版本，不改写原注释。C003保留T细胞谱系及增殖状态；C006谱系未定。聚类及参考证据均排除CGNL1和连接面板，接受身份规则冻结在新目标复测前。

Methods
使用作者提供的非负整数计数，以原提供基因全集的文库总量计算CP10K；不将decontX调整计数当作原始检测计数。主分支纳入审阅认可内皮簇全部细胞，基础QC为UMI≥200、基因≥100、可评估线粒体比例≤35%。不再要求两个身份标记。病灶先达到≥10个内皮细胞，每个病灶计算CGNL1 log1p伪bulk CP10K；连接评分为CLDN5、ESAM、TJP1、OCLN四基因同尺度均值。患者有多个合格病灶时等权平均，患者作为统计样本；不合并低细胞数病灶以凑过门槛。
主连续模型为两研究合并的研究校正模型，以及各研究独立模型；三检验单独BH校正。预设敏感性包括≥20/50细胞、参考注释分歧、逐步QC和测序深度协变量；每个分析范围内单独保留对应敏感性家族。主模型另用HC3和4,999次研究内置换核查。
匹配按病灶进行，CGNL1检测阳性/阴性细胞的log1p UMI和log1p检测基因数分别限定caliper≤0.25，不放回最近邻；每病灶≥3对，然后以患者为单位比较差异。预设阶梯包含未匹配、caliper1/0.5/0.25、至少5/10对、≥20/50细胞。少于4名患者或模型自由度不足不作推断检验。主匹配合并及队列检验在可估计范围内独立BH校正。
技术比较复用此前在发现队列冻结的20个丰度/检测率匹配基因，不根据本次P值选基因；不同基因仅比较共同患者，结果为描述性。
所有敏感性、主要阴性、原始P、FDR和家族均完整保留；事后共同患者与逐一去除患者诊断独立标记，不替代冻结主分析。

Results：审阅后主连续模型未达到名义显著
主分析20个合格病灶，18位独立患者（GSE245552=10、GSE285991=8）。
两队列合并：标准化β=-0.113404，P=0.366035，FDR=0.412958；非标准化斜率=-0.264698，95%CI [-0.869973,0.340577]。HC3 P=0.619501，研究内置换P=0.369400。
GSE245552：β=-0.464099，P=0.176639，FDR=0.412958。
GSE285991：β=+0.337936，P=0.412958，FDR=0.412958；n=8，探索性。
两队列方向不同，但斜率交互P=0.132500，当前数据不能确认队列间斜率差异。
≥10/20/50细胞的研究校正P依次为0.366035/0.640814/0.061474，患者数18/16/8。严格细胞级参考支持分支及逐步QC均未使宽阈值主结论变为正向名义显著。全部模型见continuous_models.tsv。

Results：同一组患者的未匹配信号在深度匹配后减弱
未匹配：7位患者7/7同向，平均差异+0.304250，95%CI [0.188655,0.419846]，P=0.000663，敏感性家族FDR=0.005305；精确符号翻转P=0.015625。属于小样本探索性差异。
首次宽松匹配caliper1：仍为同7位患者、同7个病灶，差异+0.081494，P=0.297165；caliper0.5：+0.097252，P=0.287128。
主匹配caliper0.25：7位患者4/7同向，差异+0.101121，95%CI [-0.110611,0.312852]，P=0.286883，FDR=0.476173，精确符号翻转P=0.343750。
因此，该变化发生于对同一批病灶重新选择深度可比细胞时，不是减少患者数引起；支持检测深度或与其相关细胞组成影响未匹配差异，不能仅凭此认定所有差异均为技术伪影。
单队列主匹配：GSE245552仅2患者（不作推断），GSE285991为5患者，P=0.476173。具体患者效应及配对数完整保留。

Results：预设≥50细胞深度校正分支保留名义显著线索
全部认可内皮簇，≥50细胞且研究/深度校正：8患者，β=+0.186717，非标准化斜率+0.596087，95%CI [0.081941,1.110233]，P=0.032313，敏感性家族FDR=0.166697。
排除双参考免疫细胞、仅保留双参考内皮细胞的对应分支P=0.038717/0.046082，FDR均为0.166697。它们高度重叠，不是三次独立验证。
按用户标准，这些达到P<0.05的结果可以进入候选列表，明确登记为小样本、名义显著、FDR未通过的探索性信号。事后HC3 P=0.090189，逐一去除患者P范围0.008009–0.140618，8次中的6次未达0.05；方向均保持正向。主结果与该局部信号须同时呈现。
未作研究校正的部分合并模型出现负向名义显著，研究校正后消失，相关完整结果仍保留，不解释为稳定的反向生物学效应。

诊断：身份变化、患者覆盖及低检测共同影响证据
仅比较这两个队列，旧宽筛选23患者β=+0.504153，P=0.026391；固定为新分析覆盖的18患者后，旧选择β=+0.464743，P=0.090130；同18患者换成审阅认可内皮簇后β=-0.113404，P=0.366035。可见既有患者覆盖变化，也有细胞构成变化；不能把全部变化归因于样本量或某个单一阈值。
≥50细胞子集的8患者，用≥10细胞病灶汇总时研究校正P=0.025093；进一步限制病灶至≥50细胞时P=0.061474。反映病例子集与多病灶汇总共同影响。该诊断是事后分解，不能用其中最小P取代预设结果。
CGNL1在认可内皮中检测稀疏，见下方全簇汇总和逐病灶表；当前证据不足以区分真实弱效应、亚型/状态异质性与测量功效不足。不能据此断言基因本身无效，或任一数据集/方法有根本错误。

Discussion与后续接口
两个队列已用于身份诊断及探索，不能称为新独立验证。小样本分支只作探索性解释；病例数不能以细胞数替代。
GSE245552所提供17套feature列表缺少可识别MT基因，线粒体比例不可评估；新输出为NA，旧0值保留为legacy，未据此过滤。尚未证明是作者预过滤或基因映射原因。面板中的真实观测0与特征缺失明确区分。
下一步按原计划进入A04：在固定认可内皮中，以排除CGNL1及连接面板的证据处理亚型/状态，先冻结注释，再检验候选关联；同时补齐其余队列身份复核及独立验证。组织层面、空间和bulk路线保留，不能因本次阴性直接作废，也不能引用其阳性掩盖本次验证失败。
湿实验接口保留“细胞身份与表达定位”“连接表型与CGNL1关联”“后续干预/救援验证”三个待验证问题；当前不把计算关联写成CGNL1调控连接的已确立机制，不据此宣称临床可用。

验收与交付
27套原始矩阵缓存；独立重算20个主分析病灶的目标和面板表达，与正式表一致；108行连续模型、27行匹配模型；核验病灶/患者对应、MT缺失、两组匹配前后病灶一致和原始整数计数。
任务退出码0，有COMPLETE.json和QA.json。首轮空样本数值类型错误已修复，原失败脚本与日志保留；修复未变更冻结规则。下载后44文件SHA256核验通过。
原始全基因稀疏缓存留在sgs10，交付包包含统计表、全细胞新注释、冻结方案、审阅证据、脚本、日志与结果图。patient_continuous_profiles.tsv中的cells、counts等数值字段是病灶均值；病灶实际整数计数见lesion_continuous_profiles.tsv，不把均值当总数。
"""
summary+='\n全内皮簇检测计数（不是独立样本量）\n'+coverage.to_string()+'\n'
(R/'复测说明.txt').write_text(summary,encoding='utf-8')
(R/'NEXT_STEP.json').write_text(json.dumps({'current_module':'review consensus applied and A05/A06 retest complete','next_original_plan_module':'A04 target-independent endothelial subtype/state analysis','positive_candidate':'prespecified n50 study-depth model, n8, P0.032313, FDR0.166697, exploratory only','primary_negative_retained':True,'wetlab_bridge_reserved':True},indent=2))
print(coverage.to_string());print('REPORT_AND_FIGURES_WRITTEN')
