
from pathlib import Path
import json
import numpy as np,pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(r'C:\Users\31937\CGNL1_A04_subtype_check_20261002');D=R/'results';(R/'figures').mkdir(exist_ok=True)
# Fixed-baseline lesion diagnostic, clearly separate from broad subtype thresholds.
baseline=pd.read_csv(D/'lesion_profiles_composition.tsv',sep='\t');sub=pd.read_csv(D/'lesion_subtype_profiles.tsv',sep='\t');sp=sub[(sub.branch=='all')&(sub.mincells==5)&sub['sample'].isin(baseline[baseline.branch=='all']['sample'])]
pt=sp.groupby(['dataset','donor','subtype'],as_index=False)[['CGNL1','junction','mean_logUMI']].mean();rows=[]
for label,g in pt.groupby('subtype'):
 row={'subtype':label,'n':len(g),'beta':np.nan,'slope':np.nan,'CI_low':np.nan,'CI_high':np.nan,'P':np.nan,'status':'not_estimable'}
 if len(g)>=4 and g.CGNL1.std()>1e-12:
  X=sm.add_constant(pd.concat([g[['CGNL1']],pd.get_dummies(g.dataset,drop_first=True,dtype=float)],axis=1),has_constant='add')
  if len(g)-X.shape[1]>=3 and np.linalg.matrix_rank(X.to_numpy())==X.shape[1]:
   m=sm.OLS(g.junction,X).fit();ci=m.conf_int().loc['CGNL1'];row.update(beta=m.params.CGNL1*g.CGNL1.std()/g.junction.std(),slope=m.params.CGNL1,CI_low=ci.iloc[0],CI_high=ci.iloc[1],P=m.pvalues.CGNL1,status='estimated')
 rows.append(row)
f=pd.DataFrame(rows);f['family']='posthoc_fixed_baseline_subtype_n5';f['FDR']=np.nan;k=f.P.notna();f.loc[k,'FDR']=multipletests(f.loc[k,'P'],method='fdr_bh')[1];f.to_csv(D/'fixed_baseline_subtype_diagnostic.tsv',sep='\t',index=False)
p=pd.read_csv(D/'patient_profiles_composition.tsv',sep='\t');p=p[p.branch=='all'].copy();m=pd.read_csv(D/'composition_adjusted_models.tsv',sep='\t');mp=pd.read_csv(D/'matched_patient_effects.tsv',sep='\t');p['short']=p.donor.str.split(':').str[-1]
colors={'GSE245552':'#d97706','GSE285991':'#2563eb'}
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(2,2,figsize=(14,10),layout='constrained')
xx=np.arange(len(p));sf=p.sinusoidal_fraction.to_numpy();lf=p.lymphatic_fraction.to_numpy()
ax[0,0].bar(xx,sf,label='Sinusoidal-like',color='#2a9d8f');ax[0,0].bar(xx,1-sf-lf,bottom=sf,label='Other blood endothelial-like',color='#8da0cb');ax[0,0].bar(xx,lf,bottom=1-lf,label='Lymphatic-like',color='#e9c46a')
ax[0,0].set(xticks=xx,xticklabels=p.short,ylabel='Fraction (equal lesion weights)',title='A  Coarse reference-like composition in 18 patients',ylim=(0,1.25));ax[0,0].tick_params(axis='x',rotation=65);ax[0,0].axvline(9.5,color='grey',ls='--',lw=1);ax[0,0].legend(fontsize=8,ncol=1,loc='upper right',frameon=False)
ax[0,1].bar(xx,p.counts,color=[colors[x] for x in p.dataset]);ax[0,1].set(xticks=xx,xticklabels=p.short,ylabel='Total provided CGNL1 counts',title='B  11 / 18 patients have at most 5 target counts',yscale='symlog',ylim=(0,230));ax[0,1].tick_params(axis='x',rotation=65);ax[0,1].axhline(5,color='grey',ls='--',lw=1)
for i,v in enumerate(p.counts):ax[0,1].annotate(str(v),(i,v),xytext=(0,3),textcoords='offset points',ha='center',fontsize=8)
keys=['study','study_sinusoidal','study_lymphatic','study_both_fractions','study_both_depth'];v=m[m.branch=='all'].set_index('model').loc[keys];yy=np.arange(5)
ax[1,0].errorbar(v.slope,yy,xerr=np.vstack([v.slope-v.CI_low,v.CI_high-v.slope]),fmt='o',capsize=3,color='#264653');ax[1,0].axvline(0,color='grey',lw=1);ax[1,0].set(yticks=yy,yticklabels=['Study','+ sinusoidal fraction','+ lymphatic fraction','+ both fractions','+ both fractions + depth'],xlabel='CGNL1 slope, 95% CI',title='C  Composition adjustment does not recover a positive association');ax[1,0].invert_yaxis()
wide=mp[(mp.branch=='all')&(mp.minpairs==3)].pivot(index='donor',columns='constraint',values='effect')
for donor,row in wide.iterrows():ax[1,1].plot([0,1],[row.depth_only,row.same_coarse_subtype],'-o',color=colors[donor.split(':')[0]],alpha=.8)
ax[1,1].axhline(0,color='grey',lw=1);ax[1,1].set(xticks=[0,1],xticklabels=['Depth only\nMean +0.101; P=0.287','Depth + same subtype\nMean +0.041; P=0.648'],ylabel='Patient junction-score difference',title='D  Same 7 patients, at least 3 matched pairs')
fig.suptitle('CGNL1 | A04 subtype, composition and information-content check',fontsize=15);fig.savefig(R/'figures/A04_checks.png',dpi=170);fig.savefig(R/'figures/A04_checks.svg');plt.close(fig)
refs=pd.read_csv(R/'identity/reference_LODO_confusion.tsv',sep='\t');refs['correct']=refs.truth.eq(refs.predicted);validation={k:float(g.loc[g.correct,'cells'].sum()/g.cells.sum()) for k,g in refs.groupby('reference')}
note="""CGNL1｜原计划A04阶段性检查：亚型、组成与有效信息量
日期：2026-10-02；范围：GSE245552和GSE285991两个已经接受主要谱系审阅的队列

主要结论
本轮不能把此前阴性归因于“只是患者人数不够”。目标转录本低检出造成的信息不足确实存在，但当前粗亚型组成校正没有恢复稳定的正向关联；同一批患者的同亚型匹配效应进一步减弱。少数原阳性患者仍保持正向，因此也不能断言CGNL1无效或所有旧信号均为伪影。

与原v1.1的对应
本轮补充第5页A04“区分内皮亚型与内皮数量”，并将其接入第6页A05患者级主检验与A06同亚型匹配。完成的是两个队列的粗亚型/组成诊断及其目标复测，不是全队列A04收口；动静脉及毛细血管细分仍以描述性证据保留，完整内皮数量变化、治疗背景解释及其余队列复核尚未完成。之前18例主检验阴性以及n8、P=0.032313/FDR=0.166697的小样本敏感性正向线索均保留，本轮不替换原检验家族。

Methods
沿用认可内皮簇的全部482和1,618细胞，未新增身份门槛。亚型聚类、参考特征选择、相似性计算均排除CGNL1、CLDN5、ESAM、TJP1、OCLN，以及PECAM1、VWF、KDR、EMCN、CDH5。每队列内HVG1500、按sample选择、PCA20、15邻居、cosine距离、Leiden0.3/0.6、随机种子20261002。
疾病参考使用既有CRCLM multiome作者内皮四亚型，正常肝参考使用GSE115469作者内皮三类标签。逐供者等权构建参考中心，在参考中选择每类前100对比特征，以Pearson相似度转移标签；每次留一供者验证均重新选特征。疾病参考SEC映射为肝窦样，Tip-like和Stalk-like合为其他血管内皮样，Lymphatic映射为淋巴内皮样。不以相似度差值淘汰细胞；全部相似度/分差保留。
亚型身份表及哈希在分层目标检验前冻结。主要谱系的外部审阅不等于本轮亚型经过独立盲审；本轮属于目标排除、参考支持的计算注释。参考LODO属于参考内部泛化检查，不是目标队列的独立身份真值验证。
主要组成检验固定此前18位患者及20病灶（每病灶≥10细胞），先病灶后患者等权汇总，逐步加入肝窦样比例、淋巴内皮样比例、两比例及深度；其余血管内皮样为组成参照。各分支5个目标效应检验组成自己的BH家族。另对双比例模型使用HC3。
亚型内连续分析以每病灶每亚型≥5细胞为宽检查，并保留3/10细胞阶梯。该门槛低于整体≥10细胞，故其他血管内皮样分支可以有21位患者，不代表新增了21位独立验证患者；固定原20病灶的事后对照单独输出并单独BH校正。
患者内匹配沿用log1p UMI与log1p基因数caliper0.25及不放回规则，比较仅深度匹配与同时要求同一粗亚型匹配；每病灶至少3对作为宽检查、至少10对作为原v1.1可比门槛。每患者仅一个差值，少于4患者不作推断检验。
另使用此前两个主要谱系参考均判定为内皮的细胞重复组成和同亚型匹配。这一敏感性在目标分层检验前记录，不按CGNL1或P值选择细胞。FDR按组成、分层阈值、匹配和影响诊断分别记录。

Results：粗亚型标记有支持，精细标签不能强行确定
GSE245552：肝窦样126、其他血管内皮样348、淋巴样8细胞。
GSE285991：肝窦样345、其他血管内皮样1,178、淋巴样95细胞。
肝窦样群具有FCN3、OIT3、FCGR2B等程序，GSE285991淋巴样簇具有CCL21、PROX1、PDPN、FLT4。动静脉、毛细血管及状态程序有梯度，但不足以给所有细胞作稳定精细命名。部分小簇出现免疫或基质程序，宽分支仍保留，双参考支持分支单独核对。正常肝参考缺少淋巴类别，并将大量疾病内皮映射到同一近门管肝窦标签，因此不把两套参考当作细亚型完全一致的证据。

Results：校正亚型比例未恢复正向总体关联
固定18位患者：
仅研究校正：β=-0.113404，P=0.366035。
再加肝窦样比例：β=-0.112772，P=0.399364。
再加淋巴样比例：β=-0.114388，P=0.378472。
同时加入两比例：β=-0.116684，P=0.404555；非标准化斜率=-0.272352，95%CI[-0.955299,0.410594]。
两比例加深度：β=-0.134481，P=0.375270。
上述5个检验的本次A04家族FDR均为0.404555。原A05同一基线P的FDR仍保留0.412958，二者校正家族不同，不能互相覆盖。
双参考支持细胞分支仍为18患者，两比例校正β=-0.174734，P=0.229388，FDR=0.245723。
因此当前数据不支持“只是不同粗亚型比例掩盖了一个普遍正向关联”。不代表已经排除所有精细状态或治疗因素。

Results：同亚型匹配的改变不来自减少患者
宽分支、至少3对：
仅深度匹配：同7位患者，4/7正向，平均差异+0.101121，P=0.286883（原A06主家族FDR=0.476173）。
加入同粗亚型：仍为完全相同7位患者，4/7正向，平均差异+0.041391，95%CI[-0.169321,0.252103]，P=0.647785，当前匹配家族FDR=0.777343。
双参考支持分支：同7位患者，同亚型差异+0.018369，P=0.843211，FDR=0.903079。
原v1.1至少10对：同亚型匹配只剩P01/P03两位，差异分别+0.151673/+0.262954，均为正向，只作描述，不作推断P值。双参考支持分支也仅2位、均正向，平均差异+0.151432。
原有两位阳性病例在亚型控制后保留，但无法据此宣称更广患者具有一致效应。

Results：亚型内分析尚未出现名义显著信号
宽分支、每病灶每亚型至少5细胞：
肝窦样13患者，β=-0.011846，P=0.933725，FDR=0.992262。
其他血管内皮样21患者，β=-0.107165，P=0.372315，FDR=0.896446。
淋巴样仅3患者，不作推断。
3/5/10细胞阶梯及双参考支持分支全部可估计模型均未达到P<0.05，不登记新的正式阳性候选。双参考支持的GSE245552肝窦样6患者CGNL1无变异，该模型不可估计，不将其写成已检验阴性。上述分析为粗分类，不能替代后续稳定精细分类或独立队列。

Results：有效信息量不足可以被直接观察
整体主分析18患者中，11位在全部合格内皮细胞里CGNL1总计数≤5，2位为0。
GSE245552的10位患者计数中位数1.5、范围0–9；其中8位≤3。
GSE285991的8位患者中位数13.5、范围2–155。
这些是作者提供整数计数的观测值，不代表真实完全不表达，也不把单细胞检出0等同生物学缺失。增加细胞或改善测量只能改善每位患者的表达估计，不能替代新增独立患者。
整体主模型逐一去除18位患者后，P=0.198431–0.578059，均未达名义显著。GSE285991单独的逐一去除仍为正向但均不显著；GSE245552多数为负向。当前整体阴性不能归为某一个患者单独造成。

Discussion：可以确定与尚不能确定的事项
已观察到患者内目标测量稀疏、严格同亚型匹配有效人数很少、宽分支中患者方向不一致；未发现粗亚型组成校正后恢复稳定正向的证据。因此合理表述是“病例有效数与测量信息均不足，且效应不稳定”，不能承诺增加病例就一定显著，也不能由阴性证明基因无效。
当前参考为疾病/正常肝资料，存在跨平台和病态迁移问题。粗类混合了部分状态；免疫/基质样细胞、亚型不确定性、治疗差异以及GSE245552线粒体特征不可评估仍需分别记录。没有依据把任一整个数据集称为有根本质量问题。
后续优先补齐原预审GSE231559:L6及发现队列的同口径复核，再评价尚未用于目标的新队列；新增资料应同时考察独立患者数、内皮覆盖及CGNL1可测变异。组织、空间、临床路线按原v1.1继续，不能把本轮两队列阴性当作全部项目失败。
湿实验接口继续保留细胞来源/定位和连接状态问题；本轮不提供促转移或连接调控机制已成立的结论。

验收
身份阶段、关联阶段退出码均0；身份哈希在目标检验前后保持一致；18患者基线表达和7患者原深度匹配效应重现上一版本；原始计数用于目标，归一化/调整参考仅用于标签。下载后44个文件逐一SHA256核验通过。完整分支、P/FDR、不可估计原因、逐患者数据及脚本随包交付。

参考来源
正常肝内皮标签与生物学背景：
https://www.nature.com/articles/s41467-018-06318-7
肝窦/大血管及淋巴相关标记：
https://www.nature.com/articles/s41586-019-1373-2
跨器官血管细胞状态及标记参考：
https://www.nature.com/articles/s41591-024-03376-x
疾病参考直接复用项目已下载CRCLM_multiome_GEX_decontaminated.h5ad的作者Cell_subtype，未将其调整计数作为本轮CGNL1检测终点。
"""
note+='\n参考内部留一供者总体准确率（不等于目标数据准确率）\n'+json.dumps(validation,ensure_ascii=False,indent=2)+'\n\n固定原20病灶、每亚型≥5细胞的事后对照\n'+f.to_string(index=False)+'\n'
(R/'A04检查说明.txt').write_text(note,encoding='utf-8')
(R/'MODULE_STATUS.json').write_text(json.dumps({'module':'A04 subtype-composition diagnosis linked to A05/A06','status':'two_reviewed_cohort_coarse_check_complete','whole_original_A04_complete':False,'new_nominal_positive_subtype_models':0,'next_original_plan':['same-framework review of remaining original cohorts including GSE231559 L6','qualify independent data by patients, EC coverage, target measurement','continue original spatial and bulk branches'],'wetlab_space_preserved':True},indent=2))
print('FIXED_BASELINE_SUBTYPE',f.to_string(index=False));print('REFERENCE_ACCURACY',validation);print('REPORT_READY')
