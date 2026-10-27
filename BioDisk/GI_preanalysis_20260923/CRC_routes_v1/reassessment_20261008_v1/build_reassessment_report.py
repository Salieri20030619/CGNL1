from pathlib import Path
import json,datetime,hashlib,html,base64,zipfile,shutil
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'reassessment_20261008_v1';O=R/'results';F=R/'figures';F.mkdir(exist_ok=True)
def read(n):return pd.read_csv(O/(n+'.tsv'),sep='\t')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
freeze=json.loads((R/'SPEC_FREEZE.json').read_text());assert sha(R/'REASSESSMENT_SPEC.json')==freeze['sha256']
exits={n:int((R/'logs'/f'{n}.exit').read_text()) for n in ['spatial','lineage','dll4','ncor1']};assert all(v==0 for v in exits.values())
for f in ['SPATIAL_SPECIFICITY_COMPLETE.json','LINEAGE_PAIRS_COMPLETE.json','DLL4_COMPLETE.json','NCOR1_DIAGNOSTICS_COMPLETE.json']:assert (O/f).exists()
target=read('R1_target_summary');rank=read('R1_control_rank');specific=read('R1_specificity_contrast');r2=read('R2_state_models');perf=read('R2_incremental_performance');pair=read('R3_paired_disease_effects');ps=read('R3_paired_summary');dll=read('R4_DLL4_target_and_controls');nc=read('R4_NCOR1_target_and_controls');controls=read('R1_selected_controls')
qa=[]
old=pd.read_csv(B/'feasibility_audit_20261008_v1/results/GSE294385_v3_patient_summary.tsv',sep='\t');old=old[old.scope.eq('baseline')&old.radius_um.eq(200)&old.branch.eq('all_QC')&old.hypothesis.eq('H1')].iloc[0];new=target[target.cohort.eq('GSE294385')&target.radius_um.eq(200)&target.model.eq('A_prior')].iloc[0]
qa.append({'check':'prior_GSE294385_H1_reproduced','passed':bool(np.isclose(old.mean_beta,new.mean_beta) and np.isclose(old.P,new.P)),'old_FDR_9_scale_region_family':old.FDR,'new_target_FDR_3_model_family':new.target_FDR,'new_benchmark_FDR_51_gene_family':new.benchmark_FDR,'FDR_families_not_interchangeable':True})
for cohort,g in target[target.radius_um.eq(200)].groupby('cohort'):
 qa.append({'check':'R1_target_BH_'+cohort,'passed':bool(np.allclose(multipletests(g.P,method='fdr_bh')[1],g.target_FDR))})
allgenes=read('R1_all_gene_models')
for (cohort,model),g in allgenes[allgenes.radius_um.eq(200)].groupby(['cohort','model']):
 valid=g.P.notna();qa.append({'check':'R1_benchmark_BH_'+cohort+'_'+model,'passed':bool(np.allclose(multipletests(g.loc[valid,'P'],method='fdr_bh')[1],g.loc[valid,'benchmark_FDR']))})
assert all(x['passed'] for x in qa)
retry=json.loads((R/'perturbation_inputs/NCOR1_RETRY_RECORD.json').read_text());assert retry['integrity'];assert sha(R/'perturbation_inputs'/retry['url'].split('/')[-1])==retry['sha256']
decisions=[
 ('R1','CGNL1具有独立、突出的壁细胞邻域关联','未通过主线资格','11患者完整模型β=0.01715、P=0.14456、FDR=0.14456；4患者队列2正2负；绝对效应在50比较基因中52百分位。中间B模型相对比较基因中位差P=0.02571、FDR=0.05886保留为名义信号。','保留原空间事实；不以其充当核心机制'),
 ('R2-EC','普遍内皮连接/骨架状态故事','未得到新测试支持','16患者三程序均未达名义显著；留一患者误差增加约10.4%–16.8%。','不继续用RNA评分替代蛋白连接功能'),
 ('R2-stroma','转成基质/CAF主线','依据不足','18患者主模型三程序均未达名义显著；配对来源5患者仅3/5上升。高表达来源本身不是机制证据。','不因表达较高直接换成CAF故事'),
 ('R2/R3-epithelial','肿瘤上皮转移适应/黏着斑候选','值得优先独立验证；尚未通过跨队列资格','5配对患者4/5上升；CNV99_noTu10同4例4/4上升。GSE285991上皮黏着斑P=.01055、FDR=.02109，LOPO误差下降22.2%；另一队列方向相反，HC3 P=.24254。','列为新的候选问题，不作为已成立的新机制故事'),
 ('R4-DLL4','Notch激活可诱导CGNL1','本次最直接的响应证据','DLL4 6h后HUVEC约1.38倍、HREC约2.09倍；各体系3个处理值均高于全部对照值；固定Notch阳性对照同步响应。','确认培养体系中的刺激响应；不外推为CRCLM细胞机制'),
 ('R4-NCOR1','任一Notch相关上游改变都会升高CGNL1','当前观察不支持','siNCOR1后CGNL1约0.98倍；NCOR1下降而NOTCH4上升。归一化表仅描述，不新增供者P。','不把NCOR1加入已确立CGNL1调控链'),
 ('总体','CGNL1是否保留为项目主角','保留候选资格；新的完整机制故事尚未成立','有真实诱导反应、疾病差异和局部状态信号，也有来源/模型/队列依赖。','撤下既定内皮免疫叙事，重新评估立项；不靠缩减终点维持故事')
]
dm=pd.DataFrame(decisions,columns=['test','question','verdict','evidence','decision']);dm.to_csv(O/'STORY_REASSESSMENT_MATRIX.tsv',sep='\t',index=False)
decision={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'task':'CGNL1 story reassessment and competing-explanation tests','tests_completed':['R1 conditional spatial association and50 expression-matched comparators','R2 three fixed lineages and3 transcription programs with patient prediction','R3 paired primary/metastatic contrasts and frozen RNA-CNV sensitivities','R4 DLL4 raw-count/author-result validation and NCOR1 descriptive perturbation'],'CGNL1_status':'retain_as_conditional_candidate; no_complete_new_mechanistic_story_ready','previous_endothelial_immune_mainline':'not_supported_as_formal_story','priority_disease_hypothesis':'epithelial metastatic adaptation/focal-adhesion association merits independent verification, not yet qualified','most_direct_response_evidence':'DLL4-induced CGNL1 in cultured HUVEC/HREC','no_inference':'Notch induction in endothelium cannot be transplanted into tumor-epithelial mechanism; RNA score nulls do not refute protein function','original_results_preserved':True,'exit_codes':exits,'scope_reduction_not_used_as_rescue':True}
(R/'REASSESSMENT_DECISION.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2));(O/'NUMERICAL_QA.json').write_text(json.dumps(qa,indent=2))

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(11.8,4.3),sharey=True)
for ax,cohort in zip(axes,['GSE294385','GSE217414']):
 v=target[target.cohort.eq(cohort)&target.radius_um.eq(200)].set_index('model').loc[['A_prior','B_focal','C_stress']];ax.axhline(0,color='#888',lw=1);ax.errorbar(range(3),v.mean_beta,yerr=[v.mean_beta-v.CI_low,v.CI_high-v.mean_beta],fmt='o-',color='#1f7887',capsize=4);ax.set_xticks(range(3),['A: prior','B: focal + region' if cohort=='GSE294385' else 'B: focal','C: + hypoxia']);ax.set_title(cohort+' | n='+str(v.n.iloc[0]));ax.set_ylabel('Mean standardized mural-neighborhood coefficient')
fig.tight_layout();fig.savefig(F/'R1_conditional_spatial_effect.png',dpi=170);plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(12,4.3),sharey=True)
patient_colors={p:plt.get_cmap('tab10')(i) for i,p in enumerate(sorted(pair.patient.unique()))}
for ax,group,title in zip(axes,['Tumor_no_Tu10','Fibroblast','CNV99_noTu10'],['Author epithelial, no Tu10','Fibroblast','RNA-CNV99 epithelial, no Tu10']):
 v=pair[pair.group.eq(group)]
 for _,row in v.iterrows():ax.plot([0,1],[row.primary,row.liver],'-o',alpha=.8,lw=1.4,color=patient_colors[row.patient],label=row.patient)
 ax.set_xticks([0,1],['Primary','Liver metastasis']);ax.set_title(title+'\nn='+str(len(v)));ax.set_ylabel('CGNL1 pseudobulk log2(CPM+1)')
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=5,frameon=False);fig.tight_layout(rect=[0,.06,1,1]);fig.savefig(F/'R3_paired_site_change.png',dpi=170);plt.close(fig)
dv=read('R4_DLL4_individual_values');nv=read('R4_NCOR1_individual_values');fig,axes=plt.subplots(1,3,figsize=(11.8,4.0),sharey=True)
for ax,cell in zip(axes[:2],['HUVEC','HREC']):
 v=dv[dv.cell_type.eq(cell)&dv.gene.eq('CGNL1')];den=v[v.condition.eq('control')].author_normalized_count.mean()+1
 for k,condition in enumerate(['control','DLL4']):
  yy=(v[v.condition.eq(condition)].author_normalized_count.to_numpy()+1)/den;ax.scatter(k+np.linspace(-.07,.07,len(yy)),yy,s=45,color=['#777','#1f7887'][k]);ax.hlines(yy.mean(),k-.18,k+.18,color='black')
 ax.set_xticks([0,1],['Fc','DLL4, 6 h']);ax.axhline(1,color='#bbb',ls=':');ax.set_title(cell+' | assay triplicates');ax.set_ylabel('(Normalized CGNL1 + 1) / control mean')
v=nv[nv.gene.eq('CGNL1')];den=v[v.condition.eq('control')].normalized_count.mean()+1
for k,condition in enumerate(['control','siNCOR1']):
 yy=(v[v.condition.eq(condition)].normalized_count.to_numpy()+1)/den;axes[2].scatter(k+np.linspace(-.07,.07,len(yy)),yy,s=45,color=['#777','#c17d35'][k]);axes[2].hlines(yy.mean(),k-.18,k+.18,color='black')
axes[2].set_xticks([0,1],['Control siRNA','siNCOR1, 48 h']);axes[2].axhline(1,color='#bbb',ls=':');axes[2].set_title('HUVEC | normalized values');fig.tight_layout();fig.savefig(F/'R4_actual_perturbations.png',dpi=170);plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4.2));v=r2[r2.lineage.eq('epithelial')&r2.program.eq('focal_adhesion')&r2.variance.eq('OLS')].set_index('scope').loc[['GSE245552','GSE285991','pooled']]
sc=v.beta/v.slope;lo=v.CI_low*sc;hi=v.CI_high*sc;ax.axvline(0,color='#888');ax.errorbar(v.beta,np.arange(3),xerr=[v.beta-lo,hi-v.beta],fmt='o',color='#1f7887',capsize=4);ax.set_yticks(range(3),[f'{idx} | n={int(row.n)}' for idx,row in v.iterrows()]);ax.set_xlabel('Standardized CGNL1 coefficient, focal-adhesion RNA program');ax.set_title('Epithelial association depends on cohort');fig.tight_layout();fig.savefig(F/'R2_epithelial_focal_adhesion.png',dpi=170);plt.close(fig)

def tbl(t):return '<div class="table">'+t.to_html(index=False,border=0,escape=True,float_format=lambda v:f'{v:.5g}',na_rep='—')+'</div>'
def img(n):return '<img alt="'+n+'" src="data:image/png;base64,'+base64.b64encode((F/n).read_bytes()).decode()+'">'
def link(n,title):return '<a href="'+html.escape(n)+'">'+html.escape(title)+'</a>'
bulk=read('R3_historical_bulk_all_models');bulk['dataset']=bulk.cohort.map({'bulk_validation':'GSE50760','second_validation':'GSE14297'});bulk['n_pairs']=np.where(bulk.dataset.eq('GSE50760'),18,bulk.n_pairs)
parts=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CGNL1重新评估与竞争解释测试</title><style>body{font-family:system-ui,"Microsoft YaHei",sans-serif;max-width:1180px;margin:36px auto;padding:0 24px;color:#253640;line-height:1.7}h1,h2{color:#155e70}h1{font-size:28px}h2{font-size:21px;margin-top:34px}.lead{background:#edf5f6;padding:20px;border-left:5px solid #155e70}.note{background:#fff5e8;padding:15px;border-left:5px solid #bf853d}.table{overflow:auto;margin:16px 0}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;vertical-align:top;padding:9px;border-bottom:1px solid #d9e2e6}th{background:#edf5f6}img{max-width:100%;height:auto}a{color:#086d88}small{color:#566b76}</style><body>',
 '<h1>CGNL1重新评估与竞争解释测试</h1><p>2026-10-08 · 既有结果保留 · 新增比较方案先记录、后测试 · sgs10</p>',
 '<div class="lead"><b>重新评估结论：保留CGNL1候选资格，撤下已预设的“内皮连接—免疫—肝内定植”主故事。新的完整机制故事尚未成立。</b><p>真实DLL4干预支持CGNL1的Notch响应性；配对上皮及单队列黏着斑关联提供了值得独立验证的疾病线索。完整条件模型中的空间关联、跨队列转录状态一致性和特异增量不足以支撑原主线。当前疾病方向的优先候选是“肿瘤上皮的转移适应”，其地位仍是待验证问题，不是已确认的新结论。</p></div>',
 '<p>本轮没有按阳性结果挑选亚群、半径或细胞来源：固定测试三类细胞、三种RNA程序、三套空间调整模型、50个表达匹配比较基因、两套既有CNV阈值和三个公开培养干预条件。新测试均为探索性；已分析队列不重命名为独立验证。</p>',
 '<h2>竞争解释的验收结果决定研究优先级</h2>',tbl(dm),
 '<h2>控制中心自身状态后，独立邻域效应明显减弱</h2>',
 '<p>200 μm固定主尺度下，GSE294385的11患者原模型A复现β=0.03858、P=0.032997。加入中心自身壁细胞分数及作者区域截距后，β=0.01931、P=0.10475；进一步加入缺氧分数后，β=0.01715、95%CI −0.00699至0.04130、6/11正向、P=0.14456、当前三模型家族FDR=0.14456。GSE217414四患者始终2正2负，完整模型均值接近零。本轮所有P均按患者效应计算，没有spot伪重复。</p>',
 img('R1_conditional_spatial_effect.png'),tbl(target[target.radius_um.eq(200)]),
 '<div class="note">FDR口径必须保留：上一轮九个尺度/区域H1检验家族中，原模型结果FDR=0.03712；本轮三种嵌套调整模型的目标家族中，该原模型FDR=0.09899；在51基因比较家族中为0.08414。三个值回答不同的多重检验口径，不能覆盖彼此。原P=0.032997不变，本轮对完整模型的判断主要依据其目标效应P=0.14456及一致性，而非仅因FDR未通过。</div>',
 '<p>比较基因在查看邻域效应前，仅按表达均值与检出率选择。完整模型中CGNL1处于正向效应80百分位、绝对效应52百分位。目标相对比较基因中位效应之差为+0.02054、P=0.05542、FDR=0.05886。中间B模型这一差值为+0.02351、P=0.02571、FDR=0.05886，属于应保留的名义显著信号；它未在完整模型或四患者队列形成一致支持。比较基因不是已验证无功能的阴性基因，排序不是生物学P值。</p>',tbl(rank),tbl(specific),
 '<h2>没有形成普遍转录状态主线，但上皮黏着斑候选必须保留</h2>',
 '<p>使用两个既有审阅队列的原始计数，按每样本每谱系≥20细胞汇总，得到23位上皮、18位成纤维/基质、16位内皮合格患者。固定模型调整研究、深度与缺氧。合并主分析中，连接支架、肌动球蛋白、黏着斑三类读出均未达到目标原始P&lt;0.05。完整模型系数、区间和检验家族见下表。</p>',tbl(r2[r2.scope.eq('pooled')&r2.variance.eq('OLS')][['lineage','program','n','beta','slope','CI_low','CI_high','P','FDR']]),
 '<p>单队列有两项需保留的结果：GSE245552上皮肌动球蛋白关联β=−0.62915、P=0.02438、FDR=0.07315，为名义显著；GSE285991上皮黏着斑关联β=1.31544、P=0.01055、FDR=0.02109，通过其相应家族FDR。后者留一患者后10/10斜率仍为正，9/10重拟合P&lt;0.05；加入CGNL1使该队列留一患者MSE降低22.2%，但仅5/10患者个体误差改善。不能因为稳健性问题将这一阳性隐藏。</p>',
 '<p>其不确定性也明确：另一队列同一黏着斑读出方向相反（β=−0.29840、P=0.22759、FDR=0.63234）；合并23患者P=0.50624、FDR=0.61639。GSE285991的HC3 P=0.24254、FDR=0.47241，最大杠杆值0.840、最大Cook距离4.273；移除该高影响患者后效应仍为正且增大，不能称阳性仅由这一患者造成。这是一条单队列候选，尚未成为稳健、跨队列的状态证据。</p>',img('R2_epithelial_focal_adhesion.png'),tbl(r2[r2.P.lt(.05)]),
 '<p>合并预测中，上皮肌动球蛋白、基质肌动球蛋白和基质黏着斑的误差改善分别只有约1.84%、0.58%、0.14%；其他六项均恶化。三项内皮读出的误差增加约10.4%–16.8%。这些结果评价固定RNA读出的预测增量，不是临床预测。</p>',tbl(perf[perf.scope.eq('pooled')]),
 '<h2>配对数据支持保留转移相关上皮差异线索，尚不足升级机制</h2>',
 '<p>GSE225857五位已接受术前治疗的配对患者中，去除Tu10基质样簇后，CGNL1在肝转移上皮群相对原发灶4/5升高，中位Δlog2(CPM+1)=+0.97428。RNA-CNV95_noTu10为同4例3/4升高、中位+0.93880；CNV99_noTu10为同4例4/4升高、中位+0.78481。两种CNV阈值是同患者敏感性分析，不计两次验证。作者内皮4例3/4升高、中位+1.39896；成纤维5例3/5升高、中位+0.30721，因此疾病差异也不是上皮专属。</p>',
 '<p>使用预定双侧精确符号检验，n=4和n=5的最小可达P分别为0.125和0.0625。不能以这种数学上不可能达到P&lt;0.05的结果否定明确方向和效应；其统计地位仍为小样本探索性，也没有因此登记新的名义显著候选。</p>',img('R3_paired_site_change.png'),tbl(ps),
 '<p>两个旧组织验证队列的阳性也完整保留。GSE50760主模型18对样本CGNL1效应+1.03078、P=0.000148、固定126候选家族FDR=0.002216；加入肝组织代理后效应+0.01569、P=0.97177。GSE14297严格结肠13对主模型效应+1.25055、P=3.12×10⁻⁶、固定候选家族FDR=0.000173；加入肝组织代理后+0.79530、P=0.08366。这说明组织组成解释值得重视，但代理协变量和部位高度相关，不能据此断言所有组织差异均为混合伪影。GSE50760效应来自log2(FPKM+1)模型，GSE14297为作者VST/spline强度尺度，两者不强行当作相同logFC单位。</p>',tbl(bulk[['dataset','model','n_pairs','logFC','P.Value','adj.P.Val','FDR_fixed126']]),
 '<h2>真实DLL4刺激确认CGNL1响应，上游作用具有条件依赖</h2>',
 '<p>本轮新取得GSE163568两种人内皮体系的原始整数计数及作者DESeq2表，并核验六列样本、对比方向和逐样本归一化因子。DLL4包被刺激6小时后，HUVEC的CGNL1约升高1.38倍，HREC约2.09倍；各体系三个处理值均高于全部对照值。固定HES1、HEY1、HEY2、NRARP、RND1、UNC5B阳性对照均同向响应。<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC8804000/">原始研究</a>；<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163568">GEO数据</a>。</p>',
 '<p>供者与重复不能混写：论文说明HUVEC来自两个供者池，HREC来自单一供者，各做三次TLA测量。作者P/FDR保留为原分析的实验重复层面记录，本轮不重新生成供者级P，也不把三列当成三位独立患者。作者xlsx中部分其他基因被Excel转为日期，本轮隔离这些异常行、不猜测修复；CGNL1及固定对照映射唯一且不受影响。</p>',img('R4_actual_perturbations.png'),tbl(dll[dll.gene.eq('CGNL1')][['cell_type','donor_material','n_control_assay_replicates','n_treated_assay_replicates','control_mean_normalized','treated_mean_normalized','log2_ratio_means_plus1','author_log2FoldChange','author_P','author_FDR','all_treated_above_all_controls']]),
 '<p>GSE237653的siNCOR1处理48小时后，CGNL1约0.98倍，未观察到DLL4实验中的上升；NCOR1自身下降，NOTCH4上升。该表为作者归一化计数，虽然值看似整数，也不能用作原始计数拟合负二项模型；列到供者/批次的配对未明确，本轮只报告个体值及效应。该观测不支持把NCOR1降低直接加入必然诱导CGNL1的链条。<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237653">GEO原始记录</a>。</p>',tbl(nc[nc.gene.isin(['CGNL1','NCOR1','NOTCH4','HEY1'])]),
 '<h2>Methods</h2>',
 '<p>本轮在REASSESSMENT_SPEC.json记录竞争解释、细胞来源、面板、模型、比较基因选择和统计家族，随后生成带时间及SHA256的SPEC_FREEZE.json。历史结果已知，记录为探索性再评估。所有新计算写入独立目录，旧方案、图表及统计值保留。</p>',
 '<p>空间测试读取GSE294385非复发肝切片以及GSE217414四张切片的计数与坐标。前者按作者区域和条码对齐，采用既有UMI≥500、基因≥200 QC及100 μm阵列坐标；后者采用组织内UMI≥500及经核验的65/spot_diameter_fullres换算。主半径200 μm，100/300 μm为固定敏感性；排除自身及无邻居中心。三套模型依次增加中心壁细胞、作者区域和缺氧。系数在患者内标准化，先切片后患者等权；患者t检验及精确符号翻转分别报告。比较基因按共享可测基因中的表达和检出匹配选择，目标及冻结身份/响应面板全部排除；50个比较基因在关联计算之前冻结。</p>',
 '<p>细胞状态测试使用两个已冻结主要谱系注释的队列，读取27份原始矩阵，每样本每谱系≥20细胞后计算伪bulk log1p(CP10k)面板均值，重复基因符号先在原计数合并。研究、平均logUMI和缺氧为固定协变量，分别在合并与单队列拟合，记录OLS、HC3、影响诊断和留一患者预测。每个程序、研究范围和方差估计口径内，对可估计谱系目标项进行BH校正。预测每次仅用训练患者拟合。配对疾病比较复用原始计数聚合和目标排除的RNA-CNV分类，固定CNV95/99及Tu10排除；以患者为单位报告配对差异和精确符号检验。</p>',
 '<p>有有效推断时，原始目标P&lt;0.05用于名义候选筛选；完整FDR、检验家族、方向、效应和样本一致性并列保存。无有效生物重复或重复层级未明确时不生成新的生物学P。原作者培养实验P/FDR单列来源和推断层级，不与患者模型合并校正。</p>',
 '<h2>Discussion：新定位需要改变假说与证据要求</h2>',
 '<p><b>当前最强的新证据是刺激响应，而非疾病机制闭环。</b>DLL4可诱导内皮CGNL1这一事实得到原始数据支持；它与既有周细胞/Notch文献相容，但该诱导本身已有机制背景，不能直接充当CRCLM创新。更不能把内皮培养体系中的诱导与肿瘤上皮配对升高拼接为“Notch→肿瘤CGNL1→肝转移”的已验证链。</p>',
 '<p><b>疾病问题的优先候选转向肿瘤上皮的转移适应，但不直接立为新主故事。</b>同患者配对变化、RNA-CNV敏感性和GSE285991黏着斑信号使它值得独立验证；另外一个队列的反向关联、HC3不确定性和混合组织影响阻止了机制升级。该候选不能靠反复挑亚型、剔除患者或更换评分来获得“通过”。</p>',
 '<p><b>不建议继续以普遍内皮连接增强或免疫排斥作为默认主线，也不因基质表达较高就转为CAF故事。</b>连接RNA、骨架RNA及黏着斑RNA均不是蛋白定位、局部张力、Rac/Rho活性或通透性的直接读出；这些RNA关联不成立不能反向证明CGNL1蛋白没有功能。与之对应，真实RNA诱导也不能证明该蛋白对肝转移有必要性或充分性。</p>',
 '<p>证据边界集中在：既有患者数据已探索；单细胞身份为冻结粗谱系，RNA-CNV提供恶性支持而非病理真值；配对患者少且已治疗；空间仍为混合spot和转录代理；培养实验供者数少、非肝窦内皮且无CRC肝转移条件。本轮尚无直接CGNL1扰动在疾病相关系统中的结局数据。因此新的完整机制故事尚未通过立项资格；CGNL1保留候选，下一阶段的意义在于取得能够区分这些解释的独立证据，而非把原问题缩小后继续叙述。</p>',
 '<h2>复现与交付</h2>',
 '<p>空间、谱系/配对、DLL4及NCOR1/诊断四个任务均退出0。NCOR1文件第一次下载中断，随后以HTTP206从记录偏移继续，完整xlsx压缩结构与SHA256验证通过；失败与重试记录保留。原xlsx日期误写行已隔离并记录，未修改原输入。</p>',
 '<p>'+link('REASSESSMENT_DECISION.json','重新评估决定')+' · '+link('REASSESSMENT_SPEC.json','固定测试方案')+' · '+link('SPEC_FREEZE.json','方案冻结记录')+' · '+link('results/STORY_REASSESSMENT_MATRIX.tsv','竞争解释验收表')+'</p>',
 '<ul>'+''.join('<li>'+link('results/'+p.name,p.name)+'</li>' for p in sorted(O.glob('*.tsv')))+ '</ul>',
 '<p>处理版公开输入与元数据见perturbation_inputs、metadata；原项目大矩阵未重复打包。scripts为本目录中的全部.py文件，日志及退出码位于logs。各检验完整原始P、FDR和不可评估原因保留在统计表中。</p></body></html>']
(R/'CGNL1_reassessment_and_tests.html').write_text('\n'.join(parts),encoding='utf-8')
from html.parser import HTMLParser
class Links(HTMLParser):
 def __init__(self):super().__init__();self.href=[]
 def handle_starttag(self,tag,attrs):
  if tag=='a':self.href.extend(v for k,v in attrs if k=='href' and not v.startswith(('https:','http:','#')))
l=Links();l.feed((R/'CGNL1_reassessment_and_tests.html').read_text());assert all((R/x).is_file() for x in l.href)
validation={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'spec_hash_matches_freeze':True,'exit_codes':exits,'numeric_checks':len(qa),'numeric_checks_all_pass':True,'relative_links_checked':len(l.href),'NCOR1_resumed_file_integrity_pass':True,'original_results_modified':False,'figures':[p.name for p in F.glob('*.png')]}
(R/'DELIVERY_VALIDATION.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation),flush=True)
