from pathlib import Path
import json,datetime,hashlib,shutil,html,base64,zipfile
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'feasibility_audit_20261008_v1';O=R/'results';F=R/'figures';F.mkdir(exist_ok=True);H=R/'historical_evidence';H.mkdir(exist_ok=True)
def read(name):return pd.read_csv(O/(name+'.tsv'),sep='\t')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
copies={
 'historical_spatial_summary.tsv':'spatial_v3/results/screen/summary.tsv',
 'historical_A05_models.tsv':'review_consensus_retest_20261002_v1/results/continuous_models.tsv',
 'historical_A04_composition.tsv':'A04_subtype_check_20261002_v1/results/composition_adjusted_models.tsv',
 'historical_A04_same_subtype.tsv':'A04_subtype_check_20261002_v1/results/same_subtype_matched_models.tsv',
 'mixture_primary_reference.tsv':'spatial_v3/background_validation/queue_20260924_reference_v3/V03_MIXTURE/summary.tsv',
 'mixture_hepatic_reference_sensitivity.tsv':'spatial_v3/background_validation/queue_20260924_reference_v3/V03_HEPREF_SENS/summary.tsv',
 'historical_relative_distance.tsv':'spatial_v3/background_validation/final_exploration_20260924_v2/F02_SPATIAL/summary.tsv',
 'historical_source_state.tsv':'spatial_v3/background_validation/final_exploration_20260924_v2/F01_SOURCE_STATE/summary.tsv',
 'historical_marker_sensitivity.tsv':'spatial_v3/background_validation/final_exploration_20260924_v2/F04_MARKER_SENS/summary.tsv',
 'v3_plan.json':'v3_pilot_20261006/plan_v3.json',
 'v3_plan_sections.json':'v3_pilot_20261006/plan_sections.json',
 'M242T_acceptance.json':'fullslice_20261007_v1/results/ANALYSIS_ACCEPTANCE.json',
 'M242T_scale_stability.tsv':'fullslice_20261007_v1/results/scale_stability.tsv',
 'M242T_annotation_amendment.json':'annotation_panel_20261007_v1/standards/ANNOTATION_AMENDMENT_01.json',
 'NicheCompass_parameters.json':'v3_pilot_20261006/results/NICHECOMPASS_PARAMETERS.json'
}
history=[]
for dest,src in copies.items():
 p=B/src;shutil.copy2(p,H/dest);history.append({'source':str(p),'copy':'historical_evidence/'+dest,'sha256':sha(p)})
(O/'HISTORICAL_EVIDENCE_MANIFEST.json').write_text(json.dumps(history,indent=2))
# Independent recheck of original spatial BH families, including planned unestimable hypotheses as P=1 only for multiplicity bookkeeping.
hist=pd.read_csv(H/'historical_spatial_summary.tsv',sep='\t');fdr=[]
expected={'secondary_genes3':3,'baseline_distance_sensitivity':8,'all_timepoints_distance_sensitivity':12,'baseline_source_comparisons':12,'all_timepoints_source_comparisons':12,'baseline_score_sensitivity':24,'all_timepoints_score_sensitivity':24}
for fam,g in hist.groupby('family'):
 if fam=='single_primary':continue
 N=expected[fam];ps=np.r_[g.P.fillna(1).to_numpy(),np.ones(N-len(g))];qs=multipletests(ps,method='fdr_bh')[1][:len(g)];qs[g.P.isna()]=np.nan
 for j,(_,r) in enumerate(g.iterrows()):fdr.append({'family':fam,'scope':r.scope,'gene':r.gene,'contrast':r.contrast,'method':r.method,'FDR_recalculated':qs[j],'FDR_old':r.FDR,'passed':bool(np.isclose(qs[j],r.FDR,equal_nan=True))})
pd.DataFrame(fdr).to_csv(O/'GSE294385_historical_FDR_reproduction.tsv',sep='\t',index=False)
checks=json.loads((O/'CHECKS.json').read_text());follow=json.loads((O/'FOLLOWUP_CHECKS.json').read_text());prov=json.loads((O/'AUTHOR_ANNOTATION_PROVENANCE.json').read_text());liver=json.loads((O/'LIVER_AUDIT_COMPLETE.json').read_text());atlas=json.loads((O/'SOURCE_ATLAS_COMPLETE.json').read_text())
exit_codes={k:int((R/'logs'/f'{k}.exit').read_text().strip()) for k in ['numeric','followup','liver','atlas']}
assert all(v==0 for v in exit_codes.values())
assert all(c['passed'] for c in checks if c['check']!='M242T_raw_counts_metadata_alignment')
assert all(c['passed'] for c in follow)
assert liver['raw_all24_match'] and liver['geometry_all_match'] and liver['historical_CGNL1_section_effects_reproduced'] and liver['historical_patient_effects_and_P_reproduced']
assert atlas['cache_matches_all'] and prov['GitHub_blob_matches'] and prov['author_RDS_to_analysis_labels_exact_match'] and prov['patient_map_matches_published_table_serialization']
assert all(x['passed'] for x in fdr)

decision_rows=[
 ('01','v3阶段定位','记录解释不合格','S0/S1技术小试；全切片计算不是S2患者级主分析。','冻结方案及局部完成标记','修正进度叙述；保留旧状态文件'),
 ('02','GSE217414测量与条码','通过','四张矩阵、QC、CGNL1计数及12个原系数独立复现。','visium_original_reproduction.tsv','保留原数值并使用校正尺度'),
 ('03','Visium物理尺度','发现实现错误','55被用于理论65um字段；旧200um实际236.36um。四张校正后pitch99.35–100.14um。','visium_geometry.tsv','独立审计重算；原输出不覆盖'),
 ('04','无邻居测量','发现实现错误','旧脚本将无邻居赋零；校正200um排除8个中心。','visium_effects_all_fixed_variants.tsv','无邻居按不可评估；不当作真实零暴露'),
 ('05','GSE294385是否可用','通过；纠正前次遗漏','24张已下载，11患者映射；341328条作者原始RDS标签与分析表完全相同。','AUTHOR_ANNOTATION_PROVENANCE.json','不能继续称该队列无区域或仅候选登记'),
 ('06','历史空间结果','复现通过','24原始矩阵、距离、CGNL1切片效应、患者P及完整FDR家族核验。','LIVER_AUDIT_COMPLETE.json','保留组织空间阳性和组成敏感性'),
 ('07','H1壁细胞邻域','单队列探索性支持','GSE294385固定9种尺度/区域组合均P<.05且本家族FDR<.05；效应小；四张GSE217414方向不一致。','GSE294385_v3_patient_summary.tsv','保留为混合spot区域关联，不升级内皮机制'),
 ('08','CGNL1细胞来源','不支持内皮专属前提','27个原始单细胞矩阵及24供者标签复核；上皮、基质和部分内皮均贡献目标。','source_atlas_summary.tsv','重新评价目标定位，不将组织表达直接归属内皮'),
 ('09','旧A05/A04/A06','复现通过；稳定功能主张不支持','旧阴性/不稳定不是程序计算错误；n8名义信号仍保留，匹配与亚型控制后不能恢复普遍正向。','A05_full_model_reproduction.tsv','保留全部阳性、阴性、敏感性及各自检验家族'),
 ('10','H2连接/黏附/T细胞','未建立跨证据链支持','16患者功能轴没有名义显著目标项；11患者空间H2主尺度未见名义显著。','EC_exploratory_axes_models.tsv','不能写成免疫排斥、外渗或屏障机制'),
 ('11','H3增量','简化诊断不支持；完整H3未验收','三个读出的留一患者MSE均恶化；本轮基线仅研究+深度，非完整亚型/丰度/应激模型。','EC_incremental_prediction_diagnostic.tsv','不能宣称特异增量或预测价值，也不把此简化诊断当完整H3反证'),
 ('12','H4调控链','技术可测；机制未建立','LIANA校正尺度后仍可运行；CGNL1与Notch读出未建立稳定关系；FlowSig及真实扰动链未执行。','LIANA_corrected_scale_retest.tsv','现有共现分数不作因果或方向证据'),
 ('13','M242T注释和邻域','计算通过；生物用途受限','原始Gene Expression计数、坐标、QC、冻结标签及1024个抽检组成均复现；目标未测。','M242T_independent_neighborhood_reproduction.tsv','可用于技术背景；不能替代CGNL1复现'),
 ('14','离散生态位/NicheCompass','未具备生物学主方法资格','六群ARI0.528；NicheCompass原生输出可复核，但未证优于基线且其空间输入继承旧尺度。','M242T_scale_stability.tsv','不调K/种子寻找阳性；不以技术完成充当生物验收'),
 ('15','SMURF/OSDR/SpaCEy等','资格未满足','HD两个公开Zenodo记录当前仍restricted且无公开文件列表；动力学/蛋白临床输入未落实。','metadata/DIRECT_METADATA_CHECK.json','这些状态不是CGNL1生物阴性'),
 ('16','故事定位','应重新评估机制主线','真实证据是多来源组织情境与弱壁细胞邻域关联；内皮特异连接—Notch—免疫功能链未成立。','全部证据','不默认扩展v3整套机制模型；不以缩减终点保存故事')
]
dec=pd.DataFrame(decision_rows,columns=['item','question','verdict','evidence','reference','decision']);dec.to_csv(O/'ITEMWISE_VERIFICATION.tsv',sep='\t',index=False)
decision={'scope':'completed itemwise feasibility audit, not completion of research project','recommendation':'Reconsider CGNL1-centered endothelial mechanistic story; retain observed multi-lineage regional/weak mural-neighborhood association as evidence, not as a replacement paper story. Do not continue full v3 mechanism expansion by default.','technical_pass_does_not_establish_biology':True,'no_gene_wide_null_claim':True,'H1_GSE294385_exploratory_FDR_signal_preserved':True,'H2_H3_H4_mechanistic_chain_not_established':True,'prior_assessment_correction':'GSE294385 author pathology and existing24 archives were omitted from previous assessment; restored and validated.','original_project_files_modified':False,'exit_codes':exit_codes,'audit_validator_correction':'Initial M242T check wrongly compared controls+gene counts with gene-only transcript_counts; original implementation was correct; corrected validator passes.','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(R/'AUDIT_DECISION.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2))

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
vis=read('visium_effects_all_fixed_variants');vp=vis[(vis.radius_um==200)&vis.policy.eq('valid_neighbors_only')&vis.geometry.eq('corrected65')]
lp=read('GSE294385_v3_patient_effects');ls=read('GSE294385_v3_patient_summary')
fig,ax=plt.subplots(1,2,figsize=(12,4.2));ax[0].axhline(0,color='#777',lw=1);ax[0].scatter(range(len(vp)),vp.beta,c=['#287f8e' if x>0 else '#bb6248' for x in vp.beta],s=65);ax[0].set_xticks(range(4),['Patient 1','Patient 2','Patient 3','Patient 4']);ax[0].set_ylabel('Standardized CGNL1-mural coefficient');ax[0].set_title('GSE217414 | corrected 200 um | descriptive')
colors=['#166b80','#c47c34','#71528b'];branches=['all_QC','all_QC_plus_author_region_intercepts','macro_tumor_QC_only']
for j,branch in enumerate(branches):
 v=lp[lp.scope.eq('baseline')&lp.radius_um.eq(200)&lp.hypothesis.eq('H1')&lp.branch.eq(branch)];rr=ls[ls.scope.eq('baseline')&ls.radius_um.eq(200)&ls.hypothesis.eq('H1')&ls.branch.eq(branch)].iloc[0];ax[1].scatter(j+np.linspace(-.10,.10,len(v)),v.beta,s=22,alpha=.65,color=colors[j]);ax[1].errorbar(j+.23,rr.mean_beta,yerr=[[rr.mean_beta-rr.CI_low],[rr.CI_high-rr.mean_beta]],fmt='D',color='black',capsize=4)
ax[1].axhline(0,color='#777',lw=1);ax[1].set_xticks(range(3),['All QC','Region adjusted','Macro tumor']);ax[1].set_title('GSE294385 | 11 mapped patients | 200 um');ax[1].set_ylabel('Patient mean standardized coefficient');fig.tight_layout();fig.savefig(F/'H1_cross_cohort_audit.png',dpi=170);plt.close(fig)
ss=read('source_atlas_summary');fig,ax=plt.subplots(figsize=(8.3,4.1));labels=['epithelial','fibroblast_stromal','endothelial'];x=np.arange(3)
for j,ds in enumerate(['GSE245552','GSE285991']):
 v=ss[ss.dataset.eq(ds)].set_index('lineage').loc[labels];ax.bar(x+(j-.5)*.32,v.median_patient_CPM,width=.30,label=ds,color=['#287f8e','#c47c34'][j])
ax.set_xticks(x,['Epithelial','Fibroblast / stromal','Endothelial']);ax.set_ylabel('Median donor pseudobulk CPM');ax.set_title('CGNL1 is not an endothelial-exclusive transcript');ax.legend(frameon=False);fig.tight_layout();fig.savefig(F/'source_atlas_audit.png',dpi=170);plt.close(fig)
pred=read('EC_incremental_prediction_diagnostic');pred=pred[pred.scope.eq('pooled')].set_index('endpoint').loc[['notch_response','adhesion','historical_junction']];fig,ax=plt.subplots(figsize=(8.3,3.8));ax.barh(['Notch response','Adhesion','Historical junction'],-100*pred.relative_MSE_improvement,color='#bb6248');ax.set_xlabel('Increase in leave-one-patient MSE after adding CGNL1 (%)');ax.set_title('16 previously analyzed patients | study + depth baseline',fontsize=12);fig.tight_layout();fig.savefig(F/'H3_incremental_diagnostic.png',dpi=170);plt.close(fig)

def tbl(t):return '<div class="table">'+t.to_html(index=False,border=0,escape=True,float_format=lambda x:f'{x:.5g}',na_rep='—')+'</div>'
def img(name):return '<img alt="'+name+'" src="data:image/png;base64,'+base64.b64encode((F/name).read_bytes()).decode()+'">'
def link(path,text):return '<a href="'+html.escape(path)+'">'+html.escape(text)+'</a>'
parts=['<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CGNL1 v3.0 逐项可行性核验</title><style>body{font-family:system-ui,"Microsoft YaHei",sans-serif;max-width:1150px;margin:40px auto;padding:0 24px;color:#23343d;line-height:1.65}h1,h2{color:#155d70}h1{font-size:28px}h2{font-size:21px;margin-top:36px}.lead{background:#edf4f5;padding:20px;border-left:5px solid #155d70}.warn{background:#fff5e9;padding:16px;border-left:5px solid #c47c34}.table{overflow:auto;margin:18px 0}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;vertical-align:top;padding:9px;border-bottom:1px solid #dbe3e7}th{background:#edf4f5}img{max-width:100%;height:auto}a{color:#096f90}small{color:#526570}code{overflow-wrap:anywhere}</style><body>']
parts+=['<h1>CGNL1 v3.0：逐项核验与故事定位判断</h1><p>2026-10-08 · 内部分析与审计 · sgs10 · 原结果保留</p>',
 '<div class="lead"><b>结论：技术流程部分成立，内皮机制主线尚未成立，应重新评估CGNL1作为该主线核心基因的定位。</b><p>核验确认了可复现的组织来源差异和弱壁细胞邻域关联，也确认了尺度实现错误、混合细胞来源，以及未得到支持的连接/Notch/免疫功能链。不能将本轮结果概括为“全部无信号”，也不能把有空间关联直接升级为内皮机制。将问题缩减成几个可计算终点，不会解决核心故事的证据缺口。</p></div>',
 '<div class="warn"><b>纠正前次综合判断：</b>GSE294385不是仅登记、尚未落地的输入。本项目已具备24张肝组织切片、11位患者映射、作者三层病理标签及历史初筛。本轮逐一复核。此前“无独立病理区域”的表述应限定于GSE217414四张小试切片，不能用于整个项目。作者区域标签可用，不等于已具备HGP或内皮细胞级定位。</div>',
 '<h2>逐项验收：计算正确性与科学支持分别判断</h2>',tbl(dec),
 '<h2>空间尺度存在真实错误，校正后四张小试切片仍方向不一致</h2>',
 '<p>GSE217414旧代码以55/spot_diameter_fullres换算坐标，而该批文件的字段对应理论65 μm直径。以阵列100 μm间距独立检查，旧换算的最近邻间距为84–85 μm，改为65后为99.35–100.14 μm；旧标称200 μm实际约236.36 μm。依据为<a href="https://www.10xgenomics.com/support/software/space-ranger/2.1/analysis/outputs/spatial-outputs">10x官方字段定义</a>。此修正仅适用于核验的四张文件；GSE294385采用已核查的100 μm阵列坐标，不套用统一65系数。</p>',
 '<p>另将无邻居中心从零暴露改为不可评估。校正200 μm下排除8个无邻居spot后，四位患者系数为+0.05245、−0.03845、−0.02526、+0.01306。固定标记逐一留出及控制中心自身壁细胞分数后，符号仍为两正两负。这些是混合spot描述性系数，本轮不生成spot级P值。LIANA按校正尺度重跑通过；原NicheCompass的物理邻接解释需要重新建立，本轮未将其重新训练包装成生物验证。</p>',tbl(read('visium_geometry')),tbl(vp[['sample','n','beta','isolated','median_neighbors']]),
 '<h2>11患者队列保留弱壁细胞邻域关联，尚不能归为内皮特异效应</h2>',
 '<p>GSE294385的作者RDS来自论文公开仓库，Git blob SHA1与历史记录吻合；341,328条三层标签逐条一致，患者映射与原补充表序列化记录一致。论文说明由两位病理医师完成标注并形成共识。24张肝切片的计数、QC、坐标及历史距离均复现。<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC13390483/">原论文与病理标注方法</a>；<a href="https://github.com/yliuup/CRC_micromets_ST/tree/main/Meta_data">作者标签仓库</a>。</p>',
 '<p>本轮先冻结新增审计分支，再用原v3面板重算：100/200/300 μm、全QC/加入区域截距/仅大转移肿瘤区域。基线非复发切片先按患者等权汇总；9项H1检验均P&lt;0.05，且通过本家族BH-FDR&lt;0.05。200 μm全QC平均标准化系数0.03858，95%CI 0.00380–0.07336，8/11正向，P=0.032997、FDR=0.037121；加入作者区域截距后系数0.03045、P=0.021448、FDR=0.037121。证据强度来自患者重复；小效应、混合spot来源和已探索队列的性质均保留。</p>',
 img('H1_cross_cohort_audit.png'),tbl(ls[ls.scope.eq('baseline')&ls.radius_um.eq(200)][['branch','hypothesis','endpoint','n_patients','positive','mean_beta','CI_low','CI_high','P','FDR']]),
 '<h2>更稳定的已有事实是组织来源差异，不能直接等同肿瘤内皮机制</h2>',
 '<p>历史基线空间结果经复核：肝小叶相对大转移肿瘤的CGNL1伪bulk差异为+1.36089，11/11正向，P=0.000350、FDR=0.001051；肿瘤基质相对肿瘤为+1.11294，11/11正向，P=0.000109、FDR=0.000434。单位均为Δlog2(CPM+1)。预设近300/远600 μm终点为+0.54145、9/11正向、P=0.016866；单一主要终点不需FDR。调整组成代理分数后的系数为+0.18153，4/11正向、P=0.288507、FDR=0.860701。伪bulk与调整后spot模型估计量不同，不将两数相减解释为“被组成解释的百分比”。</p>',
 tbl(hist[hist.scope.eq('baseline')&hist.gene.eq('CGNL1')][['contrast','method','n_patients','positive','mean_effect','P','FDR']]),
 '<p>更早的参考混合模型也存在参考依赖：主参考CGNL1效应+0.26194、P=0.09829、FDR=0.13106；替代肝细胞参考+0.35306、P=0.00784、FDR=0.01568；严格区域下两套参考均未达到原始P&lt;0.05。这些历史阳性和矛盾结果完整保留于下列附表，本轮未重拟合其去卷积模型，因此不将“组成敏感”写成“已证明全是混合伪影”。</p>',
 tbl(pd.read_csv(H/'mixture_primary_reference.tsv',sep='\t').query("gene=='CGNL1'")),tbl(pd.read_csv(H/'mixture_hepatic_reference_sensitivity.tsv',sep='\t').query("gene=='CGNL1'")),
 '<h2>全谱系原始计数复核不支持将CGNL1预先视为内皮专属</h2>',
 '<p>本轮遍历27个原始scRNA矩阵，覆盖185,330个已通过原QC并冻结谱系标签的细胞、24个供者标签；27份内皮缓存与原矩阵逐项相同。GSE245552中内皮482细胞共26个CGNL1计数，上皮1,324、成纤维/基质1,026；GSE285991中内皮1,618细胞共245个计数，上皮620、基质468。计数贡献受捕获的细胞组成影响，不能直接解释为组织内来源比例。</p>',
 '<p>供者中位伪bulk CPM进一步显示：GSE245552基质18.21、内皮3.83、上皮3.48；GSE285991内皮21.12、基质16.35、上皮3.77。支持多个谱系中的目标表达，且相对来源随队列变化。这里的上皮标签不等于经CNV确认的恶性细胞，基质标签也没有区分CAF、周细胞等亚型。</p>',img('source_atlas_audit.png'),tbl(ss),
 '<h2>连接、Notch、黏附及增量诊断没有建立完整功能链</h2>',
 '<p>按v3已有每文库≥20内皮细胞门槛重建16位患者读出，模型固定为研究+平均logUMI+CGNL1。CGNL1目标项与Notch、黏附、历史连接程序的P分别为0.93157、0.68174、0.69754，同类三读出家族FDR均0.93157；HC3结果同样未达到名义显著。新增分析使用已暴露队列，只是探索性核验。11患者空间队列的200 μm H2读出也未达到原始P&lt;0.05。</p>',
 tbl(read('EC_exploratory_axes_models').query("scope=='pooled'")[['endpoint','variance','n','beta','slope','CI_low','CI_high','P','FDR']]),
 '<p>按患者留出，加入CGNL1后Notch、黏附、连接读出的MSE分别增加26.58%、6.93%、9.01%。这是简化研究+深度基线的增量诊断，不包含完整内皮丰度、亚型和应激调整，不能冒充完整H3已被否定；它也没有提供将CGNL1升级为特异状态解释变量的正面依据。</p>',img('H3_incremental_diagnostic.png'),tbl(pred.reset_index()),
 '<p>旧A05/A04/A06的可估计系数、P和各自FDR均复现。18患者主结果β=−0.1134、P=0.3660、FDR=0.4130保留；n=8分支P=0.03231、FDR=0.1667作为名义显著探索性结果保留，HC3 P=0.09019。同亚型匹配7患者差异+0.04139、4/7正向、P=0.647785、FDR=0.777343。早期支持内皮优先级的n=3匹配线索不能取代这些后续结果。</p>',
 '<h2>Methods</h2>',
 '<p>核验读取冻结方案、原计数矩阵、条码/坐标、作者标签、参数、日志、退出码和完成标记。所有新增代码与结果写入独立审计目录。第一轮范围和补充范围分别保存在AUDIT_SCOPE.json、AUDIT_AMENDMENT_01.json；记录已看过历史结果，不宣称前瞻注册。单细胞按原已接受谱系标签聚合，重复基因符号先合并原计数；按原文库资格计算伪bulk，先文库再患者等权。GSE294385按作者区域/条码对齐原矩阵，UMI≥500、检出基因≥200；物理距离使用100 μm间距六角阵列，邻域排除自身。各v3面板至少覆盖两个成员，分数为log1p(CP10k)成员均值。</p>',
 '<p>空间审计使用患者内标准化线性系数，调整预设内皮、技术深度、成纤维、上皮及肝细胞分数；固定区域截距与肿瘤区分支。每位患者等权汇总合格切片，再对患者效应进行双侧t检验。非复发基线和含复发敏感性分别记录；H1、H2、H4读出各自形成检验家族，家族内涵盖固定尺度/区域分支；完整原始P、FDR、方向与区间全部导出。少于10患者按探索性解释；无有效重复不计算推断P。新增16患者功能轴采用OLS与HC3，三读出为同类家族，按研究范围和方差口径分别报告。留一患者预测每次仅在训练患者拟合系数；未以细胞/spot计独立n。</p>',
 '<h2>Discussion：重新评价故事定位，而非用缩减问题保存原结论</h2>',
 '<p><b>已确立的事实：</b>CGNL1在多个细胞来源中表达，组织区域差异明确，GSE294385存在弱的壁细胞邻域关联。<b>当前数据不支持：</b>把这些信号直接归因于内皮内部连接增强、Notch驱动的CGNL1调控链，或CGNL1介导的T细胞排斥/外渗。空间共现、转录程序和RNA层面本身也不等同物理接触、蛋白连接或功能作用。</p>',
 '<p>建议暂停把“CGNL1—内皮连接—Notch/免疫—肝内定植”作为已经获得预分析支持的论文主线。CGNL1目前可保留为多来源组织情境相关候选，但这是一种证据定位，不是自动替换的新论文故事。是否继续以其为核心，应取决于来源明确、疾病情境特异且有独立验证的新增证据；现有阴性后重新挑选区域、亚型或算法得到阳性，不构成重新立项依据。不能因为另一个候选有单项阳性就直接换基因。</p>',
 '<p>本轮未执行完整HGP交互、精细血管亚型功能验证、真实扰动、ATAC调控定位或临床预测，因此不声称否定CGNL1所有生物功能。已有文献支持周细胞接触经Notch调节内皮CGNL1，属于机制背景，不能作为本项目CRCLM患者机制已成立的证据：<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5852532/">原始功能研究</a>。当前核心问题是疾病情境与细胞来源的归属，以及从关联到功能链的证据跨度。</p>',
 '<p>数据边界：这些队列已被项目探索；新审计分支均为探索性，FDR通过不改变这一定位。作者病理区域与HGP是不同属性。现有公开HD元数据在2026-10-08直连复查仍为restricted、无公开文件列表；这只限制该输入，不代表整个空间队列不可用。参考标签、稀疏目标计数、跨平台差异及研究治疗背景持续限制细胞级和因果外推。</p>',
 '<h2>完整证据、代码与运行验收</h2>',
 '<p>四个计算任务退出码均为0；全部原始分析文件保持原状。首次审计有一项M242T计数对齐检查失败，经追溯发现是本次核验器把负对照计数也计入后与gene-only transcript_counts比较；原流程已正确筛选Gene Expression。该核验器问题已纠正并通过，初次记录和修正原因均保留，不登记为原分析错误。</p>',
 '<p>'+link('AUDIT_DECISION.json','本轮判断与边界')+' · '+link('results/ITEMWISE_VERIFICATION.tsv','逐项验收表')+' · '+link('AUDIT_SCOPE.json','初始核验范围')+' · '+link('AUDIT_AMENDMENT_01.json','补充核验范围')+'</p>',
 '<ul>'+''.join('<li>'+link('results/'+p.name,p.name)+'</li>' for p in sorted(O.glob('*.tsv')))+ '</ul>',
 '<p>历史证据快照及来源哈希见historical_evidence与results/HISTORICAL_EVIDENCE_MANIFEST.json。原始大矩阵仍保存在服务器项目原位置；审计包不重复装入大矩阵。</p></body></html>']
content='\n'.join(parts).replace('已暴露队列','已用于项目分析的队列')
content=content.replace(img('H1_cross_cohort_audit.png'),img('H1_cross_cohort_audit.png')+'<p><small>点代表患者；右图黑色菱形为患者均值，误差线为95%置信区间。同队列三个分支使用相同患者，不构成独立验证；左右纵轴范围不同。左图仅描述，不提供患者推断P值。</small></p>')
(R/'CGNL1_itemwise_feasibility_audit.html').write_text(content,encoding='utf-8')
validation={'exit_codes':exit_codes,'author_annotation_provenance_pass':True,'all27_raw_cache_matches':True,'all24_spatial_raw_matches':True,'historical_spatial_FDR_checks':len(fdr),'historical_spatial_FDR_all_match':all(x['passed'] for x in fdr),'initial_validator_issue_resolved':True,'HTML_written':True,'figures':[p.name for p in F.glob('*.png')],'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(R/'DELIVERY_VALIDATION.json').write_text(json.dumps(validation,indent=2))
manifest=[]
for p in sorted(R.rglob('*')):
 if p.is_file() and p.name not in ['SHA256SUMS.txt','CGNL1_audit_bundle.zip']:
  manifest.append(sha(p)+'  '+str(p.relative_to(R)))
(R/'SHA256SUMS.txt').write_text('\n'.join(manifest)+'\n')
with zipfile.ZipFile(R/'CGNL1_audit_bundle.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(R.rglob('*')):
  if p.is_file() and p.name!='CGNL1_audit_bundle.zip':z.write(p,p.relative_to(R))
print(json.dumps(validation),flush=True)
