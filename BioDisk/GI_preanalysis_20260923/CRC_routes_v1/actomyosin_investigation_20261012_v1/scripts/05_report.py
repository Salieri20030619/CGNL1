from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt,html,re
F=R/'figures';assert json.loads((R/'NUMERICAL_QA.json').read_text())['all_pass']
t=pd.read_csv(O/'D1_D2_component_models.tsv',sep='\t');v=pd.read_csv(O/'D1_D3_model_input.tsv',sep='\t');ext=pd.read_csv(O/'D4_external_models.tsv',sep='\t');ep=pd.read_csv(O/'D4_patient_profiles.tsv',sep='\t')
lo=pd.read_csv(O/'D3_leave_patient_out.tsv',sep='\t');ko=pd.read_csv(O/'D5_KO_gene_effects.tsv',sep='\t');q=json.loads((O/'D5_KO_QUALIFICATION.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
gs=S['genes'];p=t[t.scope.eq('pooled')&t.model.eq('M1')].set_index('outcome')
fig,ax=plt.subplots(figsize=(8.6,5.5))
for i,g in enumerate(gs):
 a=p.loc[g];fac=v.CGNL1.std()/v[g].std();low=a.CI_low*fac;high=a.CI_high*fac
 ax.errorbar(a.beta,i,xerr=[[a.beta-low],[high-a.beta]],fmt='o',capsize=4,color='#177e89')
ax.axvline(0,color='#888',ls='--');ax.set_yticks(range(7),gs);ax.invert_yaxis();ax.set_xlabel('Standardized CGNL1 coefficient (95% HC3 CI)');ax.set_title('Constituent-gene associations after differentiation control\n23 discovery patients; seven-gene FDR family')
fig.tight_layout();fig.savefig(F/'components.png',dpi=180,bbox_inches='tight');fig.savefig(F/'components.svg',bbox_inches='tight');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11.4,4.8))
for i,g in enumerate(gs):axes[0].plot(p.loc['minus_'+g,'beta'],i,'o',color='#177e89')
axes[0].axvline(p.loc['full','beta'],ls='--',color='#bf6c3b',label='Original full score');axes[0].axvline(0,color='#aaa',lw=.8);axes[0].set_yticks(range(7),['Drop '+g for g in gs]);axes[0].invert_yaxis();axes[0].set_xlabel('Standardized coefficient');axes[0].set_title('Leave one gene out');axes[0].legend()
axes[1].scatter(lo.beta,lo.improvement*100,color='#177e89',s=34);axes[1].axvline(0,color='#aaa',ls='--');axes[1].axhline(0,color='#aaa',ls='--');axes[1].set_xlabel('Coefficient after excluding one patient');axes[1].set_ylabel('Nested patient-held-out error improvement (%)');axes[1].set_title('All23 patient exclusions retained')
fig.tight_layout();fig.savefig(F/'robustness.png',dpi=180,bbox_inches='tight');fig.savefig(F/'robustness.svg',bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4.1));plotrows=[]
a=p.loc['full'];fac=v.CGNL1.std()/v.score_full.std();plotrows.append(('Original cohorts (n=23)',a.beta,a.CI_low*fac,a.CI_high*fac))
for definition,label in [('author_epithelial_QC','Multiome author epithelial (n=15)'),('frozen_identity_intersection','Multiome identity intersection (n=13)')]:
 a=ext[ext.definition.eq(definition)&ext.model.eq('M1')].iloc[0];u=ep[ep.definition.eq(definition)];fac=u.CGNL1.std()/u.actomyosin.std();plotrows.append((label,a.beta,a.CI_low*fac,a.CI_high*fac))
for i,(label,beta,low,high) in enumerate(plotrows):ax.errorbar(beta,i,xerr=[[beta-low],[high-beta]],fmt='o',capsize=4,color='#177e89' if i==0 else '#bd713b')
ax.axvline(0,color='#888',ls='--');ax.set_yticks(range(3),[x[0] for x in plotrows]);ax.invert_yaxis();ax.set_xlabel('Standardized coefficient (95% HC3 CI)');ax.set_title('The negative association does not reproduce in the multiome cohort')
fig.tight_layout();fig.savefig(F/'external_replication.png',dpi=180,bbox_inches='tight');fig.savefig(F/'external_replication.svg',bbox_inches='tight');plt.close(fig)
gene_order=['CGNL1','CGN']+q['core_symbols_used'];groups=['CGNL1_KO_22','CGN_KO_1G8','double_KO_21D3']
mat=ko.pivot(index='gene',columns='group',values='log2FC_CPMplus1').reindex(index=gene_order,columns=groups)
fig,ax=plt.subplots(figsize=(8.6,5.6));limit=max(1.5,float(abs(mat).max().max()));im=ax.imshow(mat,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto')
for i in range(len(mat)):
 for j in range(3):ax.text(j,i,f'{mat.iloc[i,j]:+.2f}',ha='center',va='center',color='white' if abs(mat.iloc[i,j])>.65*limit else '#182f3e')
ax.set_yticks(range(len(mat)),gene_order);ax.set_xticks(range(3),['CGNL1-KO clone22','CGN-KO clone1G8','Double-KO clone21D3']);ax.set_title('Actual MDCKII knockout: mixed transcript responses\nOne clone/genotype; descriptive, no genetic-effect P values')
fig.colorbar(im,ax=ax,label='log2[(mean KO CPM+1)/(mean WT CPM+1)]',shrink=.85)
fig.tight_layout();fig.savefig(F/'knockout.png',dpi=180,bbox_inches='tight');fig.savefig(F/'knockout.svg',bbox_inches='tight');plt.close(fig)
tr=pd.read_csv(O/'D4_transfer_predictions.tsv',sep='\t');trs=[]
for direction,g in tr.groupby('direction'):trs.append(dict(direction=direction,**performance(g)))
pd.DataFrame(trs).to_csv(O/'D4_transfer_performance.tsv',sep='\t',index=False)
def table(d):return d.to_html(index=False,border=0,float_format=lambda x:f'{x:.5g}',na_rep='未匹配/不可估计',classes='data')
comp=t[t.scope.eq('pooled')&t.model.eq('M1')&t.outcome.isin(gs+list(S['submodules']))][['outcome','n','beta','P','FDR','family']]
drops=t[t.scope.eq('pooled')&t.model.eq('M1')&t.outcome.str.startswith('minus_')][['outcome','beta','P','FDR']]
ex=ext[['definition','model','n','beta','P','FDR']]
kg=ko[ko.group.eq('CGNL1_KO_22')&ko.gene.isin(S['genes'])][['gene','WT_mean_CPM','KO_mean_CPM','log2FC_CPMplus1','delta_ln1pCP10k']]
refs=json.loads((R/'LITERATURE_EVIDENCE.json').read_text())
body='''<header><div class="eyebrow">CGNL1 · actomyosin 关联追查 · 2026-10-12</div><h1>原队列负相关稳定保留，但跨队列及扰动证据未建立普遍转录抑制规律</h1><p class="lead">完成了文献核查、七基因拆解、患者影响诊断、另一批肝转移多组学患者检验，以及真实 CGNL1 敲除表达数据复核。当前可保留原队列的轻链相关探索性线索；不能将其升级为“CGNL1 普遍压低收缩装置”的机制结论。文献更直接支持 CGNL1 对连接处肌球蛋白定位、张力传递与局部 GTPase 调节的作用。</p><p>固定靶点：CGNL1｜原队列 23 患者｜多组学上皮 15 患者，身份交集敏感性 13 患者｜敲除数据 12 个培养样本、每种 KO 基因型仅一个克隆。</p></header>
<section><h2>文献支持连接组织作用，尚未补齐“协调转录抑制”这一步</h2>
<p>2008 年原始研究显示，CGNL1 参与 GEF-H1/Tiam1 的连接募集，成熟单层与连接形成期的 RhoA/Rac1 反应不同。2017 年内皮研究同样将其联系到 GTPase 活性和黏附组织。2023 年直接结合、敲除与回补研究进一步表明，CGNL1 可募集 NM2A/NM2B 到连接处；缺失后的黏附连接碎片化与 actomyosin 作用有关。该研究的构建体表达实验中，内源 NM2 总丰度没有同步改变。<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC2555940/">2008</a>、<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5852532/">2017</a>、<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC10202830/">2023</a></p>
<p>据此，下一步机制必须区分三个层次：MYH/MYL/ROCK 转录量、连接/基底区的肌球蛋白分布、RhoA-GTP/pMLC 与实际张力。我们的七基因 RNA 分数只测到第一层。本文未将另外两层当作已经测定的结果。</p>
<p>本轮还核实了一个容易遗漏的数据入口：GSE228955 的标题仍为 WT MDCKII，但最新样本和完整 txt 矩阵包含 WT、CGN 单 KO、CGNL1 单 KO 和双 KO。Excel 附件仅有 WT，不能与完整矩阵混合使用。<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228955">GEO 原始记录</a></p></section>
<section><h2>方法：先固定组成和验证策略，再查看新结果</h2>
<p>本轮以已有负相关为起点，冻结七个组成基因 MYH9、MYH10、MYL9、MYL12A、MYL12B、ROCK1、ROCK2；逐基因、七种逐一删除配方和三个预定义子模块分别分析。M1 控制研究、深度、上皮身份和成熟分化；M2 另加再生/增殖/缺氧；组成敏感性另加成纤维、内皮和肝细胞代理。模型满秩且至少五个残差自由度时使用 HC3 t 推断。同类组成、删除配方和子模块分别构成 7、7、3 个假设的 FDR 家族。原总体程序仅用于重现参照，其既有三程序家族 FDR=0.0336 保持。</p>
<p>逐一排除所有 23 位患者，重算系数与嵌套留一患者预测，不选择性删除患者。多组学队列使用作者去污染后的非整数表达矩阵：患者内合计后转为 ln(1+CP10k)，不当作原始 UMI。主定义为作者上皮且通过原 QC，每患者至少 20 个细胞；冻结身份标记交集为敏感性。另控制作者治疗状态；患者×亚型分析使用患者及亚型固定效应，以患者为聚类单位计算标准误。</p>
<p>GSE228955 以完整 txt 文件的 12 列为准，按精确符号匹配，不用旁系同源基因替换。CGNL1 KO、CGN KO 与双 KO 分别对 WT，保留每个培养样本和效应量。每个 KO 基因型只有一个克隆，因此不据三个培养样本构造独立遗传 KO 效应的 P 值。</p></section>
<section><h2>负相关以轻链部分为主，七基因程序不是统一方向的收缩开关</h2>
<p>原 23 位患者的分化控制模型被精确重现。MYL12B 的负向关联达到名义显著（β=-0.455，P=0.0180），七基因家族 FDR=0.126；轻链子模块 β=-0.412，P=0.0319，三个子模块家族 FDR=0.0957。两者保留为探索性线索。MYH10 的关联接近零，ROCK 子模块整体也接近零。</p>
<p>线性系数分解中，MYL12A/MYL12B 合计对应约 69.7% 的净负向系数；该比例是统计分解，不能解释为因果贡献率。逐一删除任一基因后，七种配方的效应方向均仍为负，但删除 MYL9 或 MYL12B 后名义统计证据减弱。未从这些配方中择优定义新的“最佳模块”。</p>
<figure><img src="figures/components.png"><figcaption>图 1．同一分化控制模型中的单基因效应。区间为 HC3 95% 区间。</figcaption></figure>'''+table(comp)+'''<details><summary>完整逐一删除配方结果</summary>'''+table(drops)+'''</details></section>
<section><h2>原队列信号不由单个患者单独造成</h2>
<p>逐一排除 23 位患者后，全部目标系数仍为负，标准化 β 范围约 -0.552 至 -0.338。重新进行嵌套患者外预测后，误差改善范围为 16.5%–38.0%。这支持该关联在原数据中的内部稳定性；不增加独立患者数量，也不替代外部复现。</p>
<figure><img src="figures/robustness.png"><figcaption>图 2．左：七种删除基因配方；右：全部患者排除敏感性。所有排除结果均保留。</figcaption></figure></section>
<section><h2>另一批肝转移多组学患者未复现负相关</h2>
<p>作者上皮定义下 15 位患者，M1 标准化 β=+0.199，P=0.279，FDR=0.419。采用冻结身份标记交集后，13 位患者仍为正向（β=+0.320，P=0.198，FDR=0.297）。加入治疗状态也没有恢复负相关。</p>
<p>在 14 位患者、100 个患者×亚型伪批量中，控制患者与作者亚型后，系数仍偏正，P=0.315；统计单位为 14 位患者，不是 100 个亚型条目。该结果没有支持一种被亚型混合掩盖的共同负相关。</p>
<p>在总体复现失败后，补充检查了全部七个组成基因（仅 M1，两种冻结细胞定义，分别按七基因家族校正）。多组学作者上皮中的 MYL12A/MYL12B 均偏正；身份交集中的 MYL12A 为正向名义信号（β=+0.446，P=0.00903，FDR=0.0632）。因此，轻链部分也未形成跨队列一致的负向关系。该检查为事后诊断，未替换原总体终点。<a href="results/D4_postprimary_all7_components.tsv">完整七基因补充结果</a></p><p>将原 23 位患者模型直接迁移到多组学队列，加入 CGNL1 后误差增加约 104.6%。反向迁移的误差有所改善，但系数方向不一致，不能据此认定相同机制获得复现。多组学队列已用于项目其他问题，本轮是该终点的新检验，并非前瞻性独立验证。</p>
<figure><img src="figures/external_replication.png"><figcaption>图 3．相同 M1 调整下的队列比较。新的队列检验未复制原负向关联。</figcaption></figure>'''+table(ex)+'''</section>
<section><h2>真实 CGNL1 敲除出现轻链变化，但没有整体同向转录反应</h2>
<p>完整矩阵包含四组各三个培养样本，CGNL1 KO 为 clone22。其 CGNL1 RNA 均值低于 WT；该体系的蛋白敲除资格来自作者原始研究。六个可精确匹配的核心基因中，MYH9、MYH10、MYL9、MYL12B 升高，ROCK1/ROCK2 降低。MYL9 约为 WT 的 1.94 倍，MYL12B 约为 1.08 倍；这些是该克隆比较的描述性效应。</p>
<p>MYL12A 没有精确符号匹配，因此不生成虚构的完整七基因分数。六基因可测部分的均值变化仅为 +0.021（ln1pCP10k 分数），支持进一步查看具体读出，未建立普遍收缩转录装置上调。CGN KO 与双 KO 的完整结果也保留。</p>
<figure><img src="figures/knockout.png"><figcaption>图 4．三种 KO 对 WT 的描述性变化。每种 KO 仅一个克隆，色差不是独立基因型推断的显著性。</figcaption></figure>'''+table(kg)+'''</section>
<section><h2>当前应保留连接张力组织假说，停止将 RNA 负相关写成普遍转录抑制</h2>
<p>本轮给出的实际推进是：排除了单个患者独占原关联的解释，定位了轻链相关组成，完成了另一队列的未复现检查，并补入真实 CGNL1 KO 数据。原队列关联继续保留；其跨队列普遍性和直接转录抑制机制没有成立。</p>
<p>在继续固定 CGNL1 的前提下，最有判别价值的问题是：<strong>CGNL1 是否改变连接处与基底区的肌球蛋白分布/张力传递，而其 RNA 轻链读出仅随细胞状态变化？</strong>这一问题的分子基础来自既有文献，在 CRC 中仍需新的机制证据。</p>
<p>后续验证应把匹配细胞密度与上皮状态作为条件，在 CRC 模型中比较两个独立靶点扰动及全长回补，并同时查看连接 NM2A/NM2B、基底应力纤维、RhoA-GTP/pMLC 和 MYL9/MYL12B RNA。只有定位合格的结合区功能分离构建体，才能进一步区别肌球蛋白募集与其他 CGNL1 功能。MYL12B/MYL9 在这里是读出，不是更换靶点。</p></section>
<section><h2>讨论：不同证据层次和无法跨越的推断边界</h2>
<p>原数据与多组学数据的整细胞/单核平台、去污染尺度、治疗与细胞状态均有差别。它们可能影响关联，但当前尚未证明哪一种因素导致方向差异，不能用“情境依赖”自动解释验证失败。当前结果只支持保留限定队列内的探索性关联。</p>
<p>组成拆解与删除诊断是基于已知总体信号的后续分析，所有分支完整保留。轻链名义信号未通过对应 FDR 家族；这一点限制证据强度，不抹去观察到的方向。敲除体系为犬肾上皮，每基因型一个克隆，克隆与基因型混杂且非 CRC；三个培养样本不构成三个独立 KO 克隆。RNA、RhoA-GTP、pMLC、肌球蛋白定位、机械力和肝定植分别是不同终点。</p></section>
<section><h2>原始证据与复现文件</h2><ul>'''+''.join('<li><a href="'+x['url']+'">'+html.escape(x['title'])+'</a></li>' for x in refs)+'''</ul>
<p><a href="SPECIFICATION.json">冻结方案</a> · <a href="NUMERICAL_QA.json">数值核验</a> · <a href="results/D1_D2_component_models.tsv">组成与配方全模型</a> · <a href="results/D3_leave_patient_out.tsv">患者影响明细</a> · <a href="results/D4_external_models.tsv">另一队列全模型</a> · <a href="results/D5_KO_individual_gene_values.tsv">KO 个体表达值</a> · <a href="results/D5_KO_QUALIFICATION.json">KO 资格与缺失符号</a></p></section>'''
css='''body{margin:0;background:#f2f5f7;color:#19313d;font:16px/1.85 "Microsoft YaHei","Noto Sans CJK SC",sans-serif}main{max-width:1120px;margin:auto;padding:28px 22px 55px}header{background:#163f4d;color:white;border-radius:12px;padding:30px}.eyebrow{color:#a0dcda;font-size:14px}h1{font-size:28px;line-height:1.5}h2{font-size:22px;border-left:5px solid #177e89;padding-left:13px;line-height:1.5}.lead{font-size:18px}section{margin-top:22px;background:white;border:1px solid #dce5e9;border-radius:9px;padding:25px}table{font-size:12px;border-collapse:collapse;width:100%;display:block;overflow:auto}th,td{border:1px solid #dce5e9;padding:7px;text-align:left}th{background:#e9f3f4}figure{margin:24px 0}img{max-width:100%;height:auto}figcaption{font-size:14px;color:#55707e}a{color:#087c88}footer{font-size:13px;color:#647a86;padding:20px 0}summary{cursor:pointer}@media(max-width:720px){main{padding:12px}header,section{padding:17px}h1{font-size:24px}}'''
doc='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CGNL1 actomyosin 深入调查</title><style>'+css+'</style></head><body><main>'+body+'<footer>生成时间 UTC：'+datetime.datetime.now(datetime.timezone.utc).isoformat()+'</footer></main></body></html>'
(R/'CGNL1_actomyosin_investigation_report.html').write_text(doc,encoding='utf-8')
for link in re.findall(r'(?:href|src)="([^"]+)"',doc):
 if not link.startswith('http'):assert (R/link).exists(),link
decision={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'fixed_target':'CGNL1','original23_patient_signal':'internally_stable_but_component_and_model_dependent','light_chain_candidate':'MYL12B nominal association; no constituent FDR significance','external15_patient_negative_replication':'failed_direction_not_reproduced','within_subtype_replication':'not_supported','actual_KO':'mixed6_gene_descriptive_response;MYL12A_unmatched;one_clone_per_genotype','general_transcriptional_suppression_mechanism':'not_established','next_question':'CGNL1 junctional myosin localization and force redistribution in CRC, with transcript levels as separate readouts','prior_outputs_modified':False,'all_requested_investigation_modules_complete':True,'whole_CGNL1_project_complete':False}
(R/'INVESTIGATION_DECISION.json').write_text(json.dumps(decision,indent=2))
(R/'REPORT_VALIDATION.json').write_text(json.dumps({'figures':4,'relative_links_pass':True,'numerical_QA_pass':True},indent=2))
print('REPORT_COMPLETE',len(doc))
