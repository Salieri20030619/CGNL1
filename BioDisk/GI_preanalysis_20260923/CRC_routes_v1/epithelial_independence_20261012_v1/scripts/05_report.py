from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4',MPLBACKEND='Agg')
import json,datetime,html,re,hashlib
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/epithelial_independence_20261012_v1');O=R/'results';F=R/'figures'
qa=json.loads((R/'NUMERICAL_QA.json').read_text());assert qa['all_pass']
a=pd.read_csv(O/'A_CGNL1_models.tsv',sep='\t');cv=pd.read_csv(O/'A_predictive_increment.tsv',sep='\t');inp=pd.read_csv(O/'A_model_input.tsv',sep='\t')
b=pd.read_csv(O/'B_paired_effects.tsv',sep='\t');bs=pd.read_csv(O/'B_paired_summary.tsv',sep='\t')
c=pd.read_csv(O/'platform_adapted/C_patient_coefficients.tsv',sep='\t');cs=pd.read_csv(O/'platform_adapted/C_gene_summary.tsv',sep='\t')
mix=pd.read_csv(O/'lineage_mixture_diagnostic/A_CGNL1_models.tsv',sep='\t');mixcv=pd.read_csv(O/'lineage_mixture_diagnostic/A_predictive_increment.tsv',sep='\t')
blocks=pd.read_csv(O/'block_diagnostic/C_gene_summary.tsv',sep='\t')
s=json.loads((R/'SPECIFICATION.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','savefig.facecolor':'white'})
colors={'M0':'#8495a4','M1':'#087f8c','M2':'#d47737'}
outs=list(s['outcomes']);names={'focal_adhesion':'Focal adhesion','microtubule_junction':'Microtubule / junction','actomyosin':'Actomyosin'}
fig,axes=plt.subplots(1,3,figsize=(12.5,4.7),sharey=True)
for ax,out in zip(axes,outs):
 for j,model in enumerate(['M0','M1','M2']):
  row=a[a.scope.eq('pooled')&a.variance.eq('HC3')&a.model.eq(model)&a.outcome.eq(out)].iloc[0]
  ratio=inp.CGNL1.std()/inp[out].std();beta=row.beta_standardized;lo=row.CI_low*ratio;hi=row.CI_high*ratio
  ax.errorbar(beta,j,xerr=[[beta-lo],[hi-beta]],fmt='o',color=colors[model],capsize=4,lw=2)
 ax.axvline(0,color='#9ba7ae',ls='--',lw=1);ax.set_yticks(range(3),['Depth + study','+ identity + maturation','+ regeneration/cycle/hypoxia']);ax.set_title(names[out]);ax.set_xlabel('Standardized coefficient (95% HC3 CI)');ax.invert_yaxis()
fig.suptitle('CGNL1 conditional associations in metastatic epithelial profiles (23 patients)',y=1.01)
fig.tight_layout();fig.savefig(F/'A_conditional_associations.png',dpi=180,bbox_inches='tight');fig.savefig(F/'A_conditional_associations.svg',bbox_inches='tight');plt.close(fig)
predict=['CGNL1']+s['ordinary_marker_comparators']
heat=cv[cv.scope.eq('pooled')&cv.model.eq('M1')].pivot(index='predictor',columns='outcome',values='relative_MSE_improvement').reindex(index=predict,columns=outs)*100
fig,ax=plt.subplots(figsize=(8.2,5.1));im=ax.imshow(heat.to_numpy(),cmap='RdBu',norm=TwoSlopeNorm(vmin=-40,vcenter=0,vmax=40),aspect='auto')
for i in range(len(heat)):
 for j in range(3):ax.text(j,i,f'{heat.iloc[i,j]:+.1f}%',ha='center',va='center',color='white' if abs(heat.iloc[i,j])>23 else '#152c38')
ax.set_yticks(range(len(heat)),heat.index);ax.set_xticks(range(3),[names[x] for x in outs]);ax.set_title('Patient-held-out MSE improvement after differentiation controls')
fig.colorbar(im,ax=ax,label='Positive = lower held-out error (%)',shrink=.8)
fig.tight_layout();fig.savefig(F/'A_increment_marker_comparison.png',dpi=180,bbox_inches='tight');fig.savefig(F/'A_increment_marker_comparison.svg',bbox_inches='tight');plt.close(fig)
v=b[b.group.eq('Tumor_no_Tu10')];metrics=['CGNL1','CGNL1_minus_epithelial','CGNL1_minus_mature']
fig,ax=plt.subplots(figsize=(8.1,4.5))
for patient,g in v.groupby('patient'):
 vals=[g[g.metric.eq(x)].delta.iloc[0] for x in metrics];ax.plot(range(3),vals,'o-',alpha=.8,label=patient)
ax.axhline(0,color='#758b97',ls='--');ax.set_xticks(range(3),['CGNL1','CGNL1 - epithelial identity','CGNL1 - mature colonocyte']);ax.set_ylabel('Liver - primary difference\n[ln(1+CP10k), or difference of fixed scores]');ax.set_title('Paired primary / liver epithelial contrasts (5 patients)')
ax.legend(title='Patient',bbox_to_anchor=(1.02,1),loc='upper left');fig.tight_layout();fig.savefig(F/'B_paired_specificity.png',dpi=180,bbox_inches='tight');fig.savefig(F/'B_paired_specificity.svg',bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(9,4.9))
z=c[c.min_micro.eq(10)&c.gene.eq('CGNL1')]
for patient,g in z.groupby('patient'):
 vals=[g[g.model.eq(x)].beta_logCP10k.iloc[0] for x in ['C0','C1','C2']];ax.plot(range(3),vals,'o-',alpha=.65,lw=1,label=str(patient))
means=z.groupby('model').beta_logCP10k.mean().reindex(['C0','C1','C2'])
ax.plot(range(3),means,'D-',color='#172b3a',lw=3,ms=7,label='Equal-patient mean')
ax.axhline(0,color='#758b97',ls='--');ax.set_xticks(range(3),['Section','+ depth, identity, maturation','+ state and tissue mixture']);ax.set_ylabel('Micro - macro coefficient [ln(1+CP10k)]');ax.set_title('Spatial association attenuates with state and composition controls\nCommon-measured-panel exploratory sensitivity; 7 patients')
ax.legend(title='Patient',bbox_to_anchor=(1.01,1),loc='upper left');fig.tight_layout();fig.savefig(F/'C_spatial_attenuation.png',dpi=180,bbox_inches='tight');fig.savefig(F/'C_spatial_attenuation.svg',bbox_inches='tight');plt.close(fig)
# Summaries retain patient-level outcomes and complete model-family statistics.
tr=pd.read_csv(O/'A_cross_cohort_predictions.tsv',sep='\t');trans=[]
for key,g in tr.groupby(['train_cohort','test_cohort','outcome','predictor']):
 mb=float(((g.observed-g.baseline)**2).mean());ma=float(((g.observed-g.augmented)**2).mean())
 trans.append(dict(zip(['train_cohort','test_cohort','outcome','predictor'],key),n=len(g),MSE_baseline=mb,MSE_augmented=ma,relative_MSE_improvement=1-ma/mb))
trans=pd.DataFrame(trans);trans.to_csv(O/'A_cross_cohort_performance.tsv',sep='\t',index=False)
rank=cv[cv.scope.eq('pooled')&cv.model.eq('M1')&cv.outcome.eq('actomyosin')].sort_values('relative_MSE_improvement',ascending=False).copy();rank['rank']=np.arange(1,len(rank)+1);rank.to_csv(O/'A_actomyosin_benchmark_rank.tsv',sep='\t',index=False)
def table(df):
 return df.to_html(index=False,border=0,float_format=lambda x:f'{x:.5g}',na_rep='不可估计',classes='data')
at=a[a.scope.eq('pooled')&a.variance.eq('HC3')][['model','outcome','n','beta_standardized','P','FDR_target','FDR_benchmark']]
act=a[a.outcome.eq('actomyosin')&a.variance.eq('HC3')&a.model.eq('M1')][['scope','n','beta_standardized','P','FDR_target']]
mxt=mix[mix.scope.eq('pooled')&mix.variance.eq('HC3')][['model','outcome','n','beta_standardized','P','FDR_target']]
ct=cs[cs.gene.eq('CGNL1')&cs.min_micro.eq(10)][['model','n','mean_effect','positive','negative','P_signflip','FDR_target_signflip','P_ttest']]
blockt=blocks[blocks.gene.eq('CGNL1')&blocks.min_micro.eq(10)][['model','n','mean_effect','positive','negative','P_signflip','FDR_target_signflip']]
bt=bs[bs.group.eq('Tumor_no_Tu10')][['metric','n','mean_delta','median_delta','positive','negative','P','FDR']]
links=[('冻结方案','SPECIFICATION.json'),('空间面板修订','MEASUREMENT_AMENDMENT01.json'),('事后分块诊断','POST_PRIMARY_BLOCK_AUDIT.json'),('完整模型统计','results/A_all_models.tsv'),('患者外预测','results/A_heldout_predictions.tsv'),('原发/转移全部配对','results/B_paired_effects.tsv'),('原空间模型资格','results/C_model_eligibility.tsv'),('补充空间模型统计','results/platform_adapted/C_gene_summary.tsv'),('数值核验','NUMERICAL_QA.json')]
body='''<header><div class="eyebrow">CGNL1 · 独立信息验证 · 2026-10-12</div><h1>CGNL1 保留部分独立上皮状态信息，肝定植阶段特异性尚未建立</h1><p class="lead">分化控制后，CGNL1 与较低肌动蛋白—肌球蛋白转录状态相关，并改善患者外预测。原发/肝转移的小样本方向线索仍在；微/宏转移关联对细胞状态和组织组成敏感。本轮结果支持继续研究 CGNL1 的上皮骨架状态联系，尚不足以将其定义为肝定植阶段驱动因子。</p><p>固定靶点：CGNL1｜27 个肝转移文库形成 23 位合格患者的上皮测量；配对原发/肝转移 5 位患者；空间主比较 7 位患者。各分析分别使用患者作为统计单位，未将其简单相加为总样本量。</p></header>
<section><h2>方法：以分化控制、标志物对照和患者外预测检验独立信息</h2>
<p>本轮先冻结方案，再计算新目标结果。A 模块复用 GSE245552/GSE285991 的冻结上皮标签，文库至少 20 个上皮细胞，保留零表达细胞；先合计原始计数，再计算 ln(1+CP10k)，同患者文库等权平均，分别获得 13 和 10 位患者。上皮身份面板为 EPCAM/KRT8/KRT18/KRT19；成熟结肠细胞面板为 CA1、CA2、GUCA2A、SLC26A3、AQP8、CEACAM7、CLCA4。各面板排除 CGNL1。</p>
<p>M0 控制研究和深度；M1 增加上皮身份与成熟分化；M2 再增加再生、增殖和缺氧代理。主终点为黏着斑程序，微管/连接和 actomyosin 为预设次终点。主要方差估计为 HC3，常规 OLS 单列。每个模型/队列/方差方法内对三个程序作 BH 校正；六个普通上皮标志和二十个表达/检测率匹配基因组成独立比较家族。匹配基因仅依据表达与检测率选定，未查看目标关联。系数模型要求满秩、至少五个残差自由度。患者留一法中的拟合和标准化均仅使用训练患者。</p>
<p>B 模块分析 GSE225857 配对样本的 CGNL1、CGNL1 减上皮身份分数、CGNL1 减成熟分化分数，保留固定 CNV95/99 及共享作者亚型敏感性。小样本使用患者层面的精确符号检验。C 模块使用 GSE294385 非复发基线切片的独立作者病理标签，要求 Layer2 与 Layer3 一致、UMI≥500、检测基因≥200；每患者至少 10 个微转移和 20 个宏转移 spot，并保留至少 20 个微转移 spot 的敏感性。先在患者内估计区域系数，再等权汇总患者，使用精确符号翻转检验；不使用 spot 层面的推断 P 值。</p>
<p>原定 C1/C2 因空间面板覆盖不足未通过资格。补充模型统一采用所有切片共同测得的原面板基因；上皮身份保留 EPCAM/KRT8，缺氧保留 CA9/SLC2A1/VEGFA，再生保留 TACSTD2/ANXA1/CLU/LY6E，其余面板保持。该修订在原始 C0 结果已知后、补充调整结果计算前固定。后续组织组成与状态分块检查单列为事后诊断。上述不同检验家族分别保留原始 P 和 FDR。</p></section>
<section><h2>预设黏着斑主终点未获支持，actomyosin 次终点显示独立状态增量</h2>
<p>合并 23 位患者，M1 中 CGNL1 与黏着斑程序关联接近零（标准化 β=0.029，P=0.936）；加入 CGNL1 后患者外均方误差增加 22.3%。当前数据不支持以该黏着斑 RNA 主终点承接主线。</p>
<p>预设 actomyosin 次终点呈负向关联：M1 标准化 β=-0.453，P=0.0112，三个程序家族 FDR=0.0336；患者外误差降低 29.6%，23 位中 15 位预测改善。两个队列的 M1 系数均为负，单队列统计证据较弱。跨队列迁移中，误差改善分别为 65.4% 和 5.1%，幅度明显不对称。</p>
<p>M2 加入再生、增殖、缺氧后，负向效应量保留（β=-0.505），HC3 P=0.0685、FDR=0.206；患者外误差仍降低 24.4%。这一结果表现为有预测增量但统计证据依赖调整集合的状态关联。</p>
<figure><img src="figures/A_conditional_associations.png"><figcaption>图 1．合并患者的目标效应和 HC3 区间。三个程序均为 RNA 代理，不是蛋白定位、机械张力或酶活性测量。</figcaption></figure>'''+table(at)+'''<p>两个队列的分化控制模型：</p>'''+table(act)+'''</section>
<section><h2>CGNL1 优于六个普通上皮标志，但并非表达匹配基因中的唯一强预测者</h2>
<p>在相同 M1 基线下，六个普通上皮标志对 actomyosin 的加入均未降低留一患者误差，而 CGNL1 降低 29.6%。在 CGNL1、六个上皮标志及二十个表达匹配基因的 27 个预测者中，CGNL1 的该项改善列第一，但最接近的表达匹配对照为 28.9%，与 CGNL1 很接近。排名是预测性能比较，不是基因特异性的生物学 P 值，也不是更换靶点依据。</p>
<p>微管/连接程序在 GSE285991 有名义信号（M1 P=0.0329，FDR=0.0986），另一队列方向不同；合并预测改善仅约 1.5%，低于 CDX2 的约 39.6%。当前数据未将微管 RNA 程序确立为 CGNL1 的稳定独立读出。</p>
<figure><img src="figures/A_increment_marker_comparison.png"><figcaption>图 2．同一分化基线中加入各预测者后的患者外误差改善。正值为改善，负值为误差增加。</figcaption></figure>
<p>事后加入非上皮来源代理的诊断中，M1+组成的 actomyosin 负向关联仍在（β=-0.622，P=0.00101，FDR=0.00304；误差降低 41.2%）。M2+组成的效应仍负向，但 P=0.0847、FDR=0.254。该诊断支持信号不易被简单非上皮转录混合解释，同时保留了状态调整后的不确定性。</p>'''+table(mxt)+'''</section>
<section><h2>配对部位线索超过一般上皮身份变化，但尚未超过成熟分化解释</h2>
<p>5 对患者中 CGNL1 本身 4/5 在肝转移升高；扣除上皮身份分数后 5/5 同向升高，平均差 0.386。扣除成熟结肠细胞分数后只有 3/5 升高，中位差约 0.0066。因而“区别于上皮身份量”与“区别于成熟分化”得到的证据强度不同。</p>
<p>共享作者肿瘤亚型、要求每个位点每亚型至少 20 个细胞并至少两个共有亚型时，保留 4 位患者，CGNL1 及相对上皮身份的变化均 4/4 上升；相对成熟分化为 3/4。CNV 敏感性方向见完整表。5 个非零配对的双侧符号检验最小可达 P 为 0.0625，4 个配对为 0.125；小样本方向不能因未达到名义阈值而写成“无变化”，也不能称为显著机制证据。</p>
<figure><img src="figures/B_paired_specificity.png"><figcaption>图 3．每条线代表一位患者。固定分数差是特异性敏感性比较，不等同于完整多因素调整。</figcaption></figure>'''+table(bt)+'''</section>
<section><h2>微转移关联在成熟分化控制后保留，加入状态与组成后接近零</h2>
<p>符合患者内区域可识别条件的 7 位患者中，原始区域模型 6/7 显示微转移 CGNL1 较高，平均微减宏系数 +0.115。原方案仅 C0 可估计，精确 P=0.03125；按三个预定义目标模型保留家族大小后 FDR=0.09375。</p>
<p>平台适配补充模型 C1 控制深度、上皮身份与成熟分化后，方向仍为 6/7，平均效应 +0.059。完整 C2 加入再生、增殖、缺氧及成纤维/内皮/肝细胞代理后，平均效应为 -0.0053，4/7 正向、3/7 负向。两种微转移 spot 下限保留相同七位可估计患者，结果一致。</p>
<p>相对六个普通上皮标志中位标准化效应的优势，在 C1/C2 中未获名义支持。现有空间数据支持区域共变线索，当前不支持超出这些状态与组成变量的独立微/宏转移差异。</p>
<figure><img src="figures/C_spatial_attenuation.png"><figcaption>图 4．每条细线是一位患者，粗线是等患者权重均值。微转移和宏转移是病理类别，并非同一病灶的纵向早晚时间点。</figcaption></figure>'''+table(ct)+'''
<p>事后分块诊断显示，两类因素均能削弱 C1 的关联。单独增加状态块后平均效应 +0.0175；单独增加组织组成块后 +0.0079，均变为 4/7 正向。组成块在本数据中对应更大幅度的衰减；这不是组成造成全部生物学效应的因果证明。</p>'''+table(blockt)+'''</section>
<section><h2>当前主线应落在上皮收缩/骨架转录状态，肝定植特异性继续作为待验证问题</h2>
<p><strong>已获得的增量：</strong>CGNL1 在现有转移灶上皮数据中提供了超出所选普通上皮身份和成熟分化代理的部分状态信息；actomyosin 负向关联具有患者外增量，六个普通上皮标志未复制该增量。</p>
<p><strong>尚未建立的增量：</strong>CGNL1 的独立肝定植阶段作用、微转移细胞内来源、蛋白定位机制，以及促定植或抑制定植的因果方向。原定黏着斑主终点未获支持，不能将次终点正结果改写为原主假设成功。</p>
<p>在保留 CGNL1 的前提下，后续最有判别价值的是验证“CGNL1 对应何种上皮收缩/骨架组织状态”，并在真正控制细胞来源和病灶阶段的数据或实验中判断其作用。当前可用工作表述为：<strong>“CGNL1 关联的上皮骨架状态及其在结直肠癌肝转移中的情境依赖性”。</strong>暂不使用“肝定植驱动因子”作为结论。</p></section>
<section><h2>讨论：证据范围与适用边界</h2>
<p>本轮为既有数据的探索性深化，不是新患者的独立外部验证。两队列研究内按患者去重，来自不同研究登记；未进行跨研究患者身份的基因型核验。A 模块仅有转移灶样本，不能单凭 A 的状态关联判断器官或定植阶段特异性。冻结上皮标签支持来源分层，不能代替每个细胞的恶性身份确证。</p>
<p>固定面板是可审计的生物学代理，并不覆盖全部分化/应激状态；空间补充身份面板只有两个共同测得的基因。组织组成分数不是实测细胞比例，空间 spot 的 CGNL1 不能直接归属肿瘤上皮。部分状态变量也可能位于真实生物路径中，因此调整后的衰减说明独立关联未建立，不能宣称该蛋白没有功能。</p>
<p>配对部位样本有治疗背景，空间微/宏比较没有纵向定植时间。患者少、微转移覆盖有限；无数 spot 不能补足独立患者。RNA 程序负相关不等于 CGNL1 抑制肌球蛋白定位、机械收缩或蛋白活性。完整模型中证据减弱、单队列差异、比较基因接近的性能及失败资格记录均保留。</p></section>
<section><h2>数值与完整性检查已通过</h2>
<p>27 个文库的 CGNL1 和既有黏着斑测量与原记录完全一致；768 项患者外预测用未标准化训练模型独立重算，最大误差低于 3×10⁻¹⁴；目标系数经残差化方法核验，误差低于 10⁻⁸。冻结方案哈希未改变，各分析有退出码和完成标记。</p>
<p>实现检查修正了共享作者亚型选择器，排除 Tumor_* 合计行；同时保留原空间三模型家族大小，避免因两个模型不可估计而缩小 FDR 家族。初期进程等待接口不兼容已恢复，分析以记录的退出状态和完整输出验收；相关日志与修正记录保留。</p></section>
<section><h2>可复核文件</h2><ul>'''+''.join('<li><a href="'+path+'">'+label+'</a></li>' for label,path in links)+'''</ul><p>全部原始 P、FDR、效应、预测、患者级明细及不可估计模型均在结果目录。大型原始输入和缓存留在 SGS10，交付包保存来源记录和分析脚本。</p></section>'''
css='''body{background:#f1f5f7;color:#19303d;font:16px/1.85 "Microsoft YaHei","Noto Sans CJK SC",sans-serif;margin:0}main{max-width:1120px;margin:auto;padding:28px 22px 60px}header{background:#133d4c;color:white;border-radius:12px;padding:32px}.eyebrow{color:#99d9d9;font-size:14px}h1{font-size:29px;line-height:1.5}h2{font-size:22px;line-height:1.5;border-left:5px solid #087f8c;padding-left:14px}section{margin-top:22px;background:white;padding:26px;border:1px solid #dbe5ea;border-radius:10px}.lead{font-size:18px}p{margin:12px 0}table{border-collapse:collapse;width:100%;font-size:12px;line-height:1.6;display:block;overflow:auto}th,td{border:1px solid #dce5e8;padding:8px;text-align:left;overflow-wrap:anywhere}th{background:#e8f2f4}figure{margin:24px 0}img{max-width:100%;height:auto}figcaption{font-size:14px;color:#546c79}a{color:#087f8c}code{font-size:13px}footer{padding:22px;color:#607680;font-size:13px}@media(max-width:720px){main{padding:12px}section,header{padding:17px}h1{font-size:24px}}'''
doc='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CGNL1 独立信息验证</title><style>'+css+'</style></head><body><main>'+body+'<footer>CGNL1 epithelial independence analysis · '+datetime.datetime.now(datetime.timezone.utc).isoformat()+'</footer></main></body></html>'
(R/'CGNL1_independent_information_report.html').write_text(doc,encoding='utf-8')
decision={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'target':'CGNL1','ordinary_epithelial_marker_only':'not_sufficient_to_explain_all_observed_state_information','retained_signal':'negative actomyosin transcriptional association and patient-held-out increment after differentiation control','primary_focal_adhesion_endpoint':'not_supported','liver_colonization_stage_specific_independence':'not_established','spatial_signal':'positive micro/macro difference after limited differentiation adjustment; attenuates with state or mixture blocks','protein_function_causality':'not_measured','A_patients':23,'B_paired_patients':5,'C_patients':7,'post_primary_adaptations_disclosed':True,'new_goal_completed':'existing_data_independent_information_test','whole_research_project_complete':False}
(R/'INDEPENDENCE_DECISION.json').write_text(json.dumps(decision,indent=2))
for path in re.findall(r'(?:href|src)="([^"]+)"',doc):
 if not path.startswith('http'):assert (R/path).exists(),path
(R/'REPORT_VALIDATION.json').write_text(json.dumps({'figures':4,'local_links_valid':True,'numerical_QA_pass':True,'decision_saved':True},indent=2))
print('REPORT_COMPLETE',len(doc),flush=True)
