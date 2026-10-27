from pathlib import Path
import os
os.environ.update(MPLBACKEND='Agg',OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4')
import json,datetime,hashlib,html,re,zipfile
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/increment_context_20261012_v1');O=R/'results';F=R/'figures'
S=json.loads((R/'SPECIFICATION.json').read_text());qa=json.loads((R/'NUMERICAL_QA.json').read_text());assert qa['all_pass']
e1=pd.read_csv(O/'E1_measurement_slopes.tsv',sep='\t');e1=e1[e1.definition.eq('author_epithelial_QC')]
het=pd.read_csv(O/'E1_heterogeneity.tsv',sep='\t');het=het[het.definition.eq('author_epithelial_QC')]
inter=pd.read_csv(O/'E2_state_interactions.tsv',sep='\t')
perf=pd.read_csv(O/'E3_prediction_performance.tsv',sep='\t').fillna({'modifier':'none'})
state=pd.read_csv(O/'E5_CGNL1_state_models.tsv',sep='\t')
sp=pd.read_csv(O/'E5_state_predictive_increment.tsv',sep='\t');sp=sp[sp.predictor.eq('CGNL1')]
cond=pd.read_csv(O/'E6_marker_conditioned_models.tsv',sep='\t');cp=pd.read_csv(O/'E6_marker_conditioned_increment.tsv',sep='\t')
# Supplementary conditioning verification uses an independent raw formula model.
v=pd.read_csv(O/'patient_profiles_all_variants.tsv',sep='\t')
v=v[v.variant.eq('pooled_common_genes')&v.definition.isin(['discovery_frozen','author_epithelial_QC'])]
err=[]
for _,row in cond[cond.P.notna()].iterrows():
 u=v if row.scope=='all38' else v[v.dataset.ne('Ogden2025')] if row.scope=='discovery23' else v[v.dataset.eq(row.scope)]
 markers=S['ordinary_markers'] if row.marker=='all_six_joint' else [row.marker]
 formula='cycling ~ CGNL1 + depth + epithelial_identity + mature_colonocyte + C(dataset) + '+' + '.join(markers)
 fit=smf.ols(formula,data=u).fit(cov_type='HC3',use_t=True)
 err.extend([abs(fit.params.CGNL1*u.CGNL1.std()/u.cycling.std()-row.beta),abs(fit.pvalues.CGNL1-row.P)])
for family,g in cond.groupby('family'):
 g=g[g.P.notna()]
 if g.empty:continue
 order=np.argsort(g.P.to_numpy());q=np.empty(len(g));q[order]=np.minimum(1,np.minimum.accumulate((g.P.to_numpy()[order]*7/np.arange(1,len(g)+1))[::-1])[::-1])
 err.extend(abs(q-g.FDR).tolist())
ep=pd.read_csv(O/'E6_marker_conditioned_predictions.tsv',sep='\t')
assert not ep.duplicated(['scope','marker','patient']).any()
for (scope,marker),g in ep.groupby(['scope','marker']):
 row=cp[cp.scope.eq(scope)&cp.marker.eq(marker)].iloc[0]
 err.append(abs(1-((g.observed-g.augmented)**2).mean()/((g.observed-g.baseline)**2).mean()-row.improvement))
assert max(err)<1e-8
qa['marker_conditioning_independent_max_error']=max(err)
qa['marker_amendment_unchanged']=hashlib.sha256((R/'POST_PRIMARY_MARKER_CONDITIONING.json').read_bytes()).hexdigest()==json.loads((R/'POST_PRIMARY_MARKER_FREEZE.json').read_text())['sha256']
assert qa['marker_amendment_unchanged']
(R/'FINAL_NUMERICAL_QA.json').write_text(json.dumps(qa,indent=2))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
colors={'discovery23':'#087f8c','Ogden2025':'#cc7040','GSE245552':'#227d9a','GSE285991':'#718c3a'}
def save(name):
 plt.savefig(F/(name+'.png'),dpi=180,bbox_inches='tight');plt.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close()
variants=S['measurement_variants']
fig,ax=plt.subplots(figsize=(10,5.1))
for i,variant in enumerate(variants):
 for j,scope in enumerate(['discovery23','Ogden2025']):
  row=e1[e1.variant.eq(variant)&e1.scope.eq(scope)].iloc[0];ratio=row.beta/row.slope
  lo=row.CI_low*ratio;hi=row.CI_high*ratio;y=i+(j-.5)*.22
  ax.errorbar(row.beta,y,xerr=[[row.beta-lo],[hi-row.beta]],fmt='o',capsize=3,color=colors[scope],label=('Discovery:23 patients' if j==0 else 'Multiome:15 patients') if i==0 else None)
ax.axvline(0,color='#8c9a9f',ls='--');ax.set_yticks(range(4),['Original measurement','Pool libraries / all genes','Pool libraries / common genes','Within-cohort ranks']);ax.invert_yaxis();ax.set_xlabel('Standardized CGNL1 coefficient (95% HC3 CI)');ax.set_title('The actomyosin direction conflict survives measurement checks');ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.20),ncol=2,frameon=False);fig.tight_layout();save('01_measurement')
fig,axes=plt.subplots(1,2,figsize=(12,4.8))
for ax,variant,title in zip(axes,['legacy','pooled_common_genes'],['Original measurement','Common-gene measurement']):
 rows=[]
 for validation,scope,label in [('leave_one_patient','discovery23','Discovery LOPO (n=23)'),('leave_one_patient','Ogden2025','Multiome LOPO (n=15)'),('leave_one_cohort','Ogden2025','Transfer to multiome (n=15)'),('leave_one_patient','all38','Pooled LOPO (n=38)')]:
  row=perf[perf.variant.eq(variant)&perf.validation.eq(validation)&perf.scope.eq(scope)&perf.modifier.eq('none')].iloc[0];rows.append((label,100*row.additive_improvement))
 bars=ax.barh([x[0] for x in rows],[x[1] for x in rows],color=['#087f8c' if x[1]>0 else '#cc7040' for x in rows]);ax.axvline(0,color='#333',lw=.8)
 for b,(_,val) in zip(bars,rows):ax.text(val+(1.5 if val>=0 else -1.5),b.get_y()+b.get_height()/2,f'{val:+.1f}%',va='center',ha='left' if val>=0 else 'right')
 ax.set_xlim((-145,50) if variant=='legacy' else (-30,60));ax.invert_yaxis();ax.set_title(title);ax.set_xlabel('Held-out MSE improvement (%)')
fig.suptitle('Internal increment persists; a shared transportable model is not established',fontsize=13);fig.tight_layout();save('02_prediction')
names={'mature_colonocyte':'Maturation','regenerative_proxy':'Regeneration','cycling':'Cycling','hypoxia':'Hypoxia'}
fig,axes=plt.subplots(1,4,figsize=(14,4.1),sharey=True)
for ax,outcome in zip(axes,S['modifiers']):
 for i,scope in enumerate(['GSE245552','GSE285991','Ogden2025']):
  row=state[state.scope.eq(scope)&state.outcome.eq(outcome)].iloc[0];ratio=row.beta/row.slope;lo=row.CI_low*ratio;hi=row.CI_high*ratio
  ax.errorbar(row.beta,i,xerr=[[row.beta-lo],[hi-row.beta]],fmt='o',capsize=3,color=colors[scope])
 ax.axvline(0,color='#9ba7ae',ls='--');ax.set_yticks(range(3),['Discovery A (13)','Discovery B (10)','Multiome (15)']);ax.invert_yaxis();ax.set_xlabel('Standardized coefficient');ax.set_title(names[outcome])
fig.suptitle('Post-primary four-state follow-up: cycling is an exploratory candidate',fontsize=13);fig.tight_layout();save('03_states')
fig,ax=plt.subplots(figsize=(9,4.8))
p=cp[cp.scope.eq('Ogden2025')].set_index('marker').reindex(S['ordinary_markers'])
ax.barh(p.index,p.improvement*100,color=['#087f8c' if x>0 else '#cc7040' for x in p.improvement])
for i,(marker,row) in enumerate(p.iterrows()):
 q=cond[cond.scope.eq('Ogden2025')&cond.marker.eq(marker)].iloc[0]
 ax.text(row.improvement*100+(1 if row.improvement>=0 else -1),i,f'{row.improvement*100:+.1f}% | P={q.P:.3f}',ha='left' if row.improvement>=0 else 'right',va='center',fontsize=9)
ax.axvline(0,color='#777');ax.set_xlim(-50,57);ax.invert_yaxis();ax.set_xlabel('Additional LOPO MSE improvement from CGNL1 (%)');ax.set_title('Cycling increment depends on the ordinary marker in the baseline\n15 multiome patients; seven-model family FDR retained in table');fig.tight_layout();save('04_marker_conditioning')
def table(t):return t.to_html(index=False,border=0,float_format=lambda x:f'{x:.5g}',na_rep='不可估计',classes='data')
body='''<header><div class="tag">CGNL1 · 上皮增量信息追查 · 2026-10-12</div><h1>原队列骨架状态增量得到加强；跨队列通用性未建立，增殖信息形成新的候选线索</h1><p class="lead">统一患者聚合与共同基因分母后，原23患者中的actomyosin负向关联和患者外增量保留。该关系在15患者多组学队列仍偏正，治疗和预设状态交互未解释分歧。补充检验发现CGNL1与较低增殖转录状态相关，但其额外信息依赖普通上皮标志物的控制方式。</p><p>固定靶点 CGNL1｜原发现队列13+10患者｜多组学15患者，身份交集敏感性13患者｜所有推断以患者为单位。</p></header>
<section><h2>统一测量口径后，原23患者的负相关与增量仍然保留</h2>
<p>原始结果被数值重现：标准化β=−0.453，P=0.0112；原三程序家族FDR=0.0336，留一患者预测误差降低29.6%。本轮统一同一患者的合格文库计数聚合、以20,819个共同测量基因为分母后，β=−0.373，P=0.00233，本轮四分析范围家族FDR=0.00932；留一患者误差降低42.9%。这些变化支持原数据中的增量并非简单由文库平均方式或基因覆盖差异产生。</p>
<p>单独统一患者计数聚合、仍用原全转录组分母，β=−0.438，P=0.0145、FDR=0.0580，为名义显著。秩转换后仍为负，但P=0.0656、FDR=0.262；不能将所有测量敏感性写成均统计显著。</p>
<figure><img src="figures/01_measurement.png"><figcaption>图1：两类患者来源在四种测量口径中的条件关联。每个测量口径保留所有预设分析范围。原始FDR与本轮不同检验家族的FDR分别记录。</figcaption></figure>'''+table(e1[['variant','scope','n','beta','P','FDR']])+'''</section>
<section><h2>方向差异超过简单归一化问题，尚无跨队列通用的actomyosin预测模型</h2>
<p>共同基因口径下，多组学15患者的系数仍为正（β=+0.193，P=0.298、FDR=0.298）。允许三个队列具有各自的协变量效应后，多组学斜率相对两个发现队列等权平均斜率的差异为+0.700，P=0.0168、两检验家族FDR=0.0336；三队列总体异质性检验P=0.0510。原始尺度的预设斜率对比为名义显著，秩敏感性未达到名义阈值。</p>
<p>多组学队列内部重新拟合时，CGNL1仍带来小幅患者外改善：原口径5.3%，共同基因口径4.2%，但对应关联方向为正。将原两队列模型直接迁移到多组学队列，原口径误差增加104.6%，共同基因口径仍增加10.2%。相对误差随评分尺度变化，不能把104.6%到10.2%的变化当作机制恢复。</p>
<p>共同基因口径下，转移预测的基线MSE为0.171857，加入CGNL1后为0.189465；去除平均偏差后的误差分别为0.047136与0.091352。该分解仅作诊断，没有用验证集结果重新校准预测。合并38患者使用共同CGNL1效应时，留一患者误差增加8.1%。</p>
<figure><img src="figures/02_prediction.png"><figcaption>图2：正数为误差下降，负数为误差上升。LOPO为留一患者；跨队列预测只使用训练队列拟合和标准化。</figcaption></figure>'''+table(het[['variant','test','contrast','P','FDR']])+'''</section>
<section><h2>治疗和四类状态交互没有解释已观察到的方向分歧</h2>
<p>多组学队列包含8位未治疗和7位已治疗患者。合并15患者、控制深度、上皮身份与成熟分化后的CGNL1×治疗交互P=0.595。两个亚组只报告描述性效应，不以细胞数扩大样本量。原发现队列的现有患者级治疗字段不足以完成同口径联合治疗模型，未据论文标题补造标签。</p>
<p>四个预设状态的交互检验均未达到原始P&lt;0.05。增加缺氧交互在合并患者留一预测中有改善，但对多组学留出队列仍使误差增加18.1%；不能据内部改善认定它解释了跨队列差异。</p>'''+table(inter[['modifier','n','beta','P','FDR']])+'''</section>
<section><h2>增殖状态产生名义显著的新线索，但方向一致尚未转化为一致增量</h2>
<p>预设状态关联检查中发现增殖线索后，先记录事后补充方案，再对成熟、再生、增殖、缺氧四类状态完整检验。共同基因口径、多组学15患者中，CGNL1与增殖分数负相关：β=−0.583，P=0.0164、四状态家族FDR=0.0657，定位为名义显著。两个发现队列也为负方向，单队列未达到名义阈值；合并38患者β=−0.264，P=0.0403、FDR=0.161，同样为名义显著。</p>
<p>多组学队列内，加入CGNL1后留一患者误差降低20.8%，8/15患者改善，在冻结的27个预测者中排名第二。原两队列训练、预测多组学患者时误差降低9.9%。然而原23患者内部误差增加3.0%，留出GSE245552时增加9.2%，留出GSE285991时降低4.6%；合并38患者内部改善仅2.8%。这支持保留候选线索，未建立一致可迁移的增殖增量。</p>
<figure><img src="figures/03_states.png"><figcaption>图3：完整四类状态的患者层面关联；CGNL1不进入任何结果评分。成熟分化为结局时从协变量中移除自身。</figcaption></figure>'''+table(state[['scope','outcome','n','beta','P','FDR_target']])+'''</section>
<section><h2>增殖增量受到普通上皮标志物控制方式的影响</h2>
<p>普通标志物VIL1对多组学增殖分数的留一患者改善为54.6%，高于CGNL1的20.8%。因此追加六个原定普通标志物逐一控制及六者联合控制，完整保留七种条件。</p>
<p>加入VIL1后，CGNL1效应仍为负（β=−0.378），P=0.147、七模型家族FDR=0.205；额外预测误差改善为11.6%。加入KRT20后，P=0.475、FDR=0.554，额外预测误差增加23.2%。加入CDH1、CLDN3、CLDN4或CDX2后仍有原始P&lt;0.05的名义信号，各FDR=0.0832。因而不能把增殖线索表述为已全面超越普通成熟上皮标志。</p>
<p>15患者的六标志物联合模型仅4个残差自由度，按方案不作推断。合并38患者的联合模型可估计，β=−0.209，P=0.194、FDR=0.194，额外留一患者改善仅0.9%；原23患者联合模型P=0.167、FDR=0.691，预测改善12.0%。预测改善和目标效应的推断证据分别呈现。</p>
<figure><img src="figures/04_marker_conditioning.png"><figcaption>图4：在相同基础协变量加一个普通上皮标志物之后，CGNL1的额外患者外增量；P对应条件模型的CGNL1目标效应。</figcaption></figure>'''+table(cond[['scope','marker','n','df','beta','P','FDR']])+'''</section>
<section><h2>方法</h2>
<p>本轮使用GSE245552和GSE285991既有冻结上皮标签、合格文库伪批量缓存以及Ogden多组学作者去污染矩阵。发现队列分别13和10患者，合格文库至少20个上皮细胞；同患者文库在新敏感性中合计计数，深度以细胞数加权。多组学主定义为作者上皮且满足既有QC，身份标记交集另作敏感性。表达按ln(1+CP10k)计算；分数为冻结基因的均值。共享基因分母由两类矩阵共同测得的20,819个基因组成，面板覆盖及历史值重建均核验。</p>
<p>M1控制深度、EPCAM/KRT8/KRT18/KRT19上皮身份及CA1/CA2/GUCA2A/SLC26A3/AQP8/CEACAM7/CLCA4成熟分化代理，合并分析增加研究效应。系数使用HC3稳健标准误及t推断，满秩且残差自由度至少5；异质性模型允许队列特异的协变量与CGNL1斜率，检验总体斜率差异及多组学相对两个发现队列等权均值的差异。状态交互使用队列内标准化，保留各队列独立主效应；治疗模型包含治疗主效应及CGNL1交互。</p>
<p>患者留一与队列留一预测的拟合、均值、标准差均仅来自训练数据。预测模型至少3个残差自由度；在成对模型均合格时比较MSE，保留全部患者预测。四类状态及普通标志物条件检验的新增范围、已知结果和冻结时间分别记录为事后方案；其结局、对照及不可估计模型均完整输出。</p>
<p>筛选采用目标效应原始P&lt;0.05，FDR在预先记录的同类假设家族内计算；P显著而FDR未通过者称为名义显著。测量斜率每定义/口径四分析范围；异质性两对比；状态交互四状态；四状态目标关联按范围分家族；标志物条件模型七个；表达匹配基因沿用原冻结20基因，不按新结果重选。单个预设治疗交互不作FDR。预测误差不生成推断性P。</p></section>
<section><h2>讨论与下一步：保留上皮状态主线，尚不能确定分歧原因</h2>
<p>本轮把已有事实与待解释问题进一步分开：原23患者的骨架状态增量可重现，统一测量口径后保留；多组学队列也有小幅内部信息增量，但actomyosin方向不同。当前数据不支持一种已建立的通用负向机制，且未确定治疗或某类上皮状态是分歧的原因。</p>
<p>平台与队列完全混杂：前者为整细胞原始计数，后者为单核去污染后的非整数矩阵，后者没有可用的raw或layer原始计数。统一计数聚合和基因分母不能消除核/胞质捕获差别或去污染处理影响，不能据此宣称平台因素已排除。跨队列患者身份未通过基因型独立核验，沿用研究内患者去重；本轮所有分析均复用旧患者，并非新的独立验证。条件关联不等于细胞自主机制，状态代理、RNA增殖分数、蛋白定位、真实增殖率与肝定植是不同终点。</p>
<p>增殖线索可以按原始P&lt;0.05进入候选清单，但证据仍为事后探索性、普通标志物依赖，不能替代未成功复现的actomyosin终点。下一步真正有判别力的是固定CGNL1、冻结actomyosin和增殖两个读出及普通标志物基线，在未参与当前筛选的新患者或匹配平台材料中检验；若进入机制实验，应将连接处肌球蛋白定位、增殖和RNA读出分别测量。当前记录不把这些尚未实施的工作标为完成。</p></section>
<section><h2>内部验收与完整复现材料</h2>
<p>旧来源文件哈希未变；原β及29.6%增量被精确重现。独立公式模型重算系数与P、独立BH实现、所有预测误差及CGNL1四状态留出预测通过核验。E3有6个整分支因训练样本自由度不足不可估计，均位于10患者队列的附加状态交互模型；没有仅保留部分患者预测的分支。</p>'''
links=[('冻结主方案','SPECIFICATION.json'),('增殖等四状态补充方案','POST_PRIMARY_STATE_INCREMENT.json'),('普通标志物条件补充方案','POST_PRIMARY_MARKER_CONDITIONING.json'),('最终数值验收','FINAL_NUMERICAL_QA.json'),('来源哈希','SOURCE_FINGERPRINTS.json'),('测量全部定义/模型','results/E1_measurement_slopes.tsv'),('异质性完整检验','results/E1_heterogeneity.tsv'),('全部状态交互','results/E2_state_interactions.tsv'),('全部留出预测','results/E3_heldout_predictions.tsv'),('预测模型资格','results/E3_prediction_coverage.tsv'),('表达对照模型','results/E4_benchmark_models.tsv'),('四状态全部模型','results/E5_all_state_models.tsv'),('四状态预测与27对照','results/E5_state_predictive_increment.tsv'),('普通标志物条件模型','results/E6_marker_conditioned_models.tsv'),('普通标志物条件预测','results/E6_marker_conditioned_increment.tsv')]
body+='<ul>'+''.join('<li><a href="'+p+'">'+label+'</a></li>' for label,p in links)+'</ul><p>数据来源：<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245552">GSE245552</a>；<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285991">GSE285991</a>；<a href="https://github.com/EfremovaLab/CRC_metastasis">Ogden研究作者代码与数据说明</a>。</p></section>'
css='''body{margin:0;background:#f2f5f7;color:#18313c;font:16px/1.85 "Microsoft YaHei","Noto Sans CJK SC",sans-serif}main{max-width:1140px;margin:auto;padding:28px 22px 55px}header{background:#163f4d;color:#fff;border-radius:12px;padding:30px}.tag{font-size:14px;color:#a4dcda}h1{font-size:29px;line-height:1.5}.lead{font-size:18px}section{margin-top:22px;background:#fff;border:1px solid #dce5e9;border-radius:9px;padding:25px}h2{font-size:22px;line-height:1.5;border-left:5px solid #087f8c;padding-left:12px}table{display:block;overflow:auto;width:100%;font-size:12px;border-collapse:collapse}td,th{border:1px solid #dce5e9;padding:7px;text-align:left}th{background:#e9f3f4}figure{margin:24px 0}img{max-width:100%;height:auto}figcaption{font-size:14px;color:#56707e}a{color:#087f8c}@media(max-width:720px){main{padding:12px}header,section{padding:16px}h1{font-size:23px}}'''
doc='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CGNL1 上皮增量信息继续调查</title><style>'+css+'</style></head><body><main>'+body+'</main></body></html>'
(R/'CGNL1_increment_context_report.html').write_text(doc,encoding='utf-8')
for p in re.findall(r'(?:src|href)="([^"]+)"',doc):
 if not p.startswith('http'):assert (R/p).is_file(),p
decision=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),target='CGNL1',
 task='epithelial_information_increment_context_investigation',
 original_increment='reproduced; retained after patient count pooling and common gene normalization',
 general_actomyosin_negative_relation='not_supported_across_cohorts',
 cohort_slope_contrast=dict(measurement='pooled_common_genes',P=.0167865,FDR=.033573),
 cause_of_disagreement='unresolved; platform and cohort confounded; therapy and4 state interactions not nominal',
 cycling_candidate='nominal external association; post-primary; conditional marker-dependent; no consistent cross-cohort increment',
 whole_project_complete=False,this_investigation_complete=True,
 next_step='Independent matched-platform epithelial validation of fixed CGNL1,actomyosin,cycling readouts and ordinary-marker conditional baselines.',
 results_directory=str(R))
# Use exact table statistics in machine-readable decision.
q=het[het.variant.eq('pooled_common_genes')&het.test.eq('external_vs_equal_discovery_mean')].iloc[0]
decision['cohort_slope_contrast'].update(P=float(q.P),FDR=float(q.FDR),difference=float(q.contrast))
(R/'INVESTIGATION_DECISION.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2))
(R/'REPORT_VALIDATION.json').write_text(json.dumps(dict(figures=4,relative_links_valid=True,numerical_QA_pass=True,report_bytes=len(doc.encode('utf-8'))),indent=2))
print('REPORT_COMPLETE',flush=True)
