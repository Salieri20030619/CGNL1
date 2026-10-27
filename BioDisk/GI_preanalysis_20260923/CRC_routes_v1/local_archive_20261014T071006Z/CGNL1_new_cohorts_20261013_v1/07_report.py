from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MPLBACKEND='Agg')
import json,datetime,hashlib,re,html
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1');O=R/'results';F=R/'figures'
S=json.loads((R/'SPECIFICATION.json').read_text());qa=json.loads((R/'NUMERICAL_QA.json').read_text());assert qa['all_pass']
models=pd.read_csv(O/'N1_CGNL1_models.tsv',sep='\t');perf=pd.read_csv(O/'N1_predictive_increment.tsv',sep='\t');v=pd.read_csv(O/'primary_patient_profiles.tsv',sep='\t')
assert set(v.dataset)==set(S['primary_cohorts']),'Both frozen new cohorts must finish before final report'
main=models[models.definition.eq('main')&models.min_cells.eq(20)]
mp=perf[perf.definition.eq('main')&perf.min_cells.eq(20)]
transfer=pd.read_csv(O/'N3_transfer_performance.tsv',sep='\t');decomp=pd.read_csv(O/'N3_transfer_error_decomposition.tsv',sep='\t');liver=pd.read_csv(O/'N4_liver_descriptive.tsv',sep='\t')
coverage=pd.read_csv(O/'primary_cell_coverage.tsv',sep='\t')
registry=json.loads((R/'CANDIDATE_REGISTRY.json').read_text())
for row in registry['candidates']:
 if row['dataset'] in S['primary_cohorts']:
  row.update(status='raw_matrix_and_epithelial_validation_complete',qualified_patients=int(v[v.dataset.eq(row['dataset'])&v.definition.eq('main')].patient.nunique()))
 if row['dataset']=='GSE158692':row['reported_patients']=5
 if row['dataset']=='SCP2079_Tian2023':row['epithelial_pairs']='must qualify separately;23 refers to paired whole-biopsy scRNAseq in paper'
(R/'CANDIDATE_REGISTRY.json').write_text(json.dumps(registry,ensure_ascii=False,indent=2))
atlas=json.loads((R/'metadata/ATLAS_SOURCE_AUDIT.json').read_text());atlas_read=sum(x['status']=='read' for x in atlas['records'])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
colors={'GSE200997':'#137e87','GSE132465':'#b36c39'}
def save(name):
 plt.savefig(F/(name+'.png'),dpi=180,bbox_inches='tight');plt.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close()
fig,axes=plt.subplots(1,2,figsize=(12,4.5))
for ax,out in zip(axes,['actomyosin','cycling']):
 for i,(ds,model) in enumerate([(d,m) for d in S['primary_cohorts'] for m in ['M0','M1']]):
  row=main[main.dataset.eq(ds)&main.model.eq(model)&main.outcome.eq(out)].iloc[0]
  ratio=row.beta/row.slope;lo=row.CI_low*ratio;hi=row.CI_high*ratio
  ax.errorbar(row.beta,i,xerr=[[row.beta-lo],[hi-row.beta]],fmt='o' if model=='M1' else 's',color=colors[ds],capsize=3)
 ax.axvline(0,color='#8b999e',ls='--');ax.set_yticks(range(4),[d+' / '+m for d in S['primary_cohorts'] for m in ['M0','M1']]);ax.invert_yaxis();ax.set_title(out.capitalize());ax.set_xlabel('Standardized CGNL1 coefficient (95% HC3 CI)')
fig.suptitle('New primary-CRC cohorts: depth-only (M0) and differentiation-adjusted (M1)',fontsize=13);fig.tight_layout();save('01_independent_models')
fig,axes=plt.subplots(1,2,figsize=(12,4.6))
for ax,ds in zip(axes,S['primary_cohorts']):
 z=mp[mp.dataset.eq(ds)&mp.model.isin(['M0','M1','M1_VIL1','M1_KRT20'])]
 xx=np.arange(4);w=.34
 for j,out in enumerate(['actomyosin','cycling']):
  q=z[z.outcome.eq(out)].set_index('model').reindex(['M0','M1','M1_VIL1','M1_KRT20']);vals=q.improvement.to_numpy()*100
  bars=ax.bar(xx+(j-.5)*w,vals,w,label=out,color=['#137e87','#b36c39'][j])
  for b,val in zip(bars,vals):
   if np.isfinite(val):ax.annotate(f'{val:+.1f}',(b.get_x()+b.get_width()/2,val),xytext=(0,3 if val>=0 else -4),textcoords='offset points',ha='center',va='bottom' if val>=0 else 'top',fontsize=8)
 ax.axhline(0,color='#66767d');ax.set_xticks(xx,['Depth','+identity/mature','+VIL1','+KRT20'],rotation=15);ax.set_title(ds);ax.set_ylabel('LOPO MSE improvement (%)');ax.margins(y=.18);ax.legend()
fig.suptitle('Positive values indicate lower held-out prediction error');fig.tight_layout();save('02_increment')
fig,axes=plt.subplots(1,2,figsize=(12,4.2))
for ax,ds in zip(axes,S['primary_cohorts']):
 z=coverage[coverage.dataset.eq(ds)&coverage.definition.eq('main')].sort_values('patient')
 ax.bar(np.arange(len(z)),z.cells,color=['#137e87' if x>=20 else '#b36c39' for x in z.cells]);ax.set_yscale('log');ax.axhline(20,color='#a34827',ls='--');ax.set_xticks(np.arange(len(z)),z.patient.str.split(':').str[-1],rotation=90);ax.set_ylabel('Qualified epithelial cells (log scale)');ax.set_title(ds+' patient coverage')
fig.tight_layout();save('03_coverage')
def table(t):return t.to_html(index=False,border=0,float_format=lambda x:f'{x:.5g}',na_rep='不可估计',classes='data')
def f(x):return f'{x:.4g}'
def row(ds,model,out):return main[main.dataset.eq(ds)&main.model.eq(model)&main.outcome.eq(out)].iloc[0]
def cv(ds,model,out):return mp[mp.dataset.eq(ds)&mp.model.eq(model)&mp.outcome.eq(out)].iloc[0]
cohort_results=''
for ds in S['primary_cohorts']:
 n=int(row(ds,'M1','cycling').n)
 cohort_results+='<h3>'+ds+'：'+str(n)+'位合格患者</h3>'
 for out,label in [('actomyosin','actomyosin'),('cycling','增殖')]:
  m0=row(ds,'M0',out);m1=row(ds,'M1',out);p=cv(ds,'M1',out)
  cohort_results+=f'<p>{label}：仅控制深度的M0系数β={f(m0.beta)}，P={f(m0.P)}、FDR={f(m0.FDR)}；控制深度、上皮身份和成熟分化的M1系数β={f(m1.beta)}，P={f(m1.P)}、FDR={f(m1.FDR)}。M1中加入CGNL1后的留一患者预测误差改善为{100*p.improvement:+.1f}%，{int(p.patients_improved)}/{int(p.n)}患者预测改善。负的改善值表示整体误差增加。</p>'
counts={d:int(v[v.dataset.eq(d)&v.definition.eq('main')].patient.nunique()) for d in S['primary_cohorts']}
transfer_note=''
for _,z in decomp[(decomp.train!='old_discovery23')&decomp.outcome.eq('cycling')].iterrows():
 transfer_note+=f'<p>{z.train} → {z.test}：原始MSE改善{100*z.improvement:+.1f}%；去均值误差改善{100*z.centered_improvement:+.1f}%。平方偏差减少量={z.bias_squared_reduction:.5g}，去均值误差减少量={z.centered_error_reduction:.5g}。</p>'
body=f'''<header><div class="tag">CGNL1 · 新队列搜查与独立患者验证 · 2026-10-13</div><h1>增殖跨队列预测出现双向增量，队列内部独立性尚未稳定复现</h1><p class="lead">本轮完成公开来源搜查、患者与细胞类型资格核实，以及GSE200997和GSE132465的固定双终点验证。两队列原报告16和23位患者，本轮分别有{counts['GSE200997']}和{counts['GSE132465']}位达到预设上皮覆盖标准。两者之间的增殖模型迁移出现双向正向误差改善，且去均值误差也下降，形成值得保留的候选预测线索。两队列在身份与成熟分化控制后均未改善内部留一患者预测，目标效应也未达到名义显著；两类证据并列呈现。结果用于评价一般上皮状态信息。</p><p>固定靶点：CGNL1｜固定终点：actomyosin与增殖转录分数｜统计单位：患者｜完整阴性、敏感性及预测结果保留。</p></header>
<section><h2>新患者验证优先使用独立原始计数，并保留转移性CRC配对治疗候选</h2>
<p>对大型CRC图谱的{atlas['loaders']}个来源加载入口进行核查，成功读取{atlas_read}个，其余{atlas['loaders']-atlas_read}个访问失败均保留在审计表。结合GEO、原始论文和作者仓库，形成{len(registry['candidates'])}项候选/重复来源登记；没有将整合图谱的细胞总数或复用研究当作新的患者样本量。</p>
<p>GSE200997为独立美国原发CRC患者队列，GSE132465为韩国原发CRC患者队列。选择这两者依据患者数、计数矩阵和元数据可用性，未查看新队列CGNL1目标关联后择优纳入。原有两个发现队列和Ogden多组学队列保持为历史参照。</p>
<p>本轮新增的高优先级后续候选为SCP2079（Tian等，2023）：论文报告23位转移性BRAFV600E CRC患者的治疗前与第15天配对单细胞数据，包含上皮细胞。该队列能检验状态随治疗的患者内变化，但单个部位、上皮可评估配对数和矩阵尺度仍需单独审核。来源加载代码提示公开epi_matrix为归一化表达；在实际矩阵资格核验前不当作原始UMI，也不直接称为23例肝转移队列。</p>
<p>其他候选包括E-MTAB-12022（6位患者的生长模式队列）、GSE158692（5位患者、6个转移灶，含EPCAM阳性癌细胞与CAF分选）、Stanford/Sathe队列（来源加载清单6份肝转移样本）及GSE97693（多区域单细胞多组学，需要按患者和模态去重）。GSE164522因CD45免疫富集不能承担上皮独立信息验证；GSE274321的3个PDO5处理条件也不构成3位独立患者。</p>
<p><a href="CANDIDATE_REGISTRY.json">完整候选登记及来源链接</a> · <a href="metadata/ATLAS_SOURCE_AUDIT.json">图谱来源与失败访问审计</a></p></section>
<section><h2>新未治疗队列保留低增殖联系，但分化控制后的增量证据减弱</h2>
<p>GSE200997中仅控制测序深度的CGNL1–增殖关联达到校正后显著（β=−0.485，P=0.00164、FDR=0.00328）；加入上皮身份与成熟分化后，方向仍负，但P=0.278、FDR=0.556，患者外误差增加21.2%。这一结果保留了较低增殖状态的联系，当前未支持超出所选身份和分化代理的独立增量。</p>
{cohort_results}
<figure><img src="figures/01_independent_models.png"><figcaption>图1：M0与M1的CGNL1目标系数及HC3区间。标准化效应以各队列、各终点的实际标准差计算。</figcaption></figure>
'''+table(main[['dataset','model','outcome','n','beta','P','FDR']])+'''</section>
<section><h2>两队列均未改善校正后的内部预测，跨队列增殖预测出现正向结果</h2>
<p>所有队列均预设M0、M1、M1+VIL1、M1+KRT20和六普通标志物联合模型；未按显著性选择“最佳”条件。模型资格取决于患者数、满秩及至少5个残差自由度，不足者保留为不可估计。留一患者预测比较相同基线加CGNL1的变化；均方误差与改善患者数同时呈现，避免只按多数患者改善作总体增量判断。</p>
<figure><img src="figures/02_increment.png"><figcaption>图2：正值表示加入CGNL1降低留一患者预测误差，负值表示误差上升。没有将预测误差差异生成细胞层面P值。</figcaption></figure>'''+table(mp[['dataset','model','outcome','n','improvement','patients_improved']])+'''
<p>提高上皮覆盖门槛到50细胞后，GSE200997保留11位患者，M1增殖系数转为正（β=+0.455，P=0.666、FDR=0.666），负向关系对覆盖门槛敏感；GSE132465在两种门槛下均保留22位患者，主结论一致。身份敏感性与全部门槛结果保留在完整统计表。</p><p>跨队列预测完全使用训练来源拟合；旧23患者训练和两新队列互相迁移全部保留。该结果与新队列内部条件模型分别解读；跨来源的均值偏移和平台差别仍能影响整体误差。</p>'''+table(transfer)+'''<p>两新队列之间，增殖预测的原始MSE分别改善5.4%与32.2%，去均值误差也分别降低20.7%与17.0%。因此，正向预测结果不全由整体均值偏移解释。GSE200997到GSE132465方向的总MSE减少主要来自偏差降低，同时也包含患者间误差下降；反向迁移中平均偏差变差，但去均值误差下降抵消了这一变化。</p><p>这些结果支持保留跨队列预测的候选线索；其证据属性是预测误差的描述性改善，不称为名义显著或稳定独立效应。应与两队列内部M1目标效应及留一患者预测结果一起解释。以下为已冻结预测的事后描述性误差分解，恒等式为MSE＝平方平均误差＋去均值误差；没有用验证集重新拟合或校准，也不产生新的P值。</p>'''+transfer_note+table(decomp[decomp.outcome.eq('cycling')])+'''</section>
<section><h2>两组现有肝转移数据可补充方向，但不能替代独立增量检验</h2>
<p>以既有固定、排除CGNL1的谱系标记重评上皮覆盖，GSE178318的6位肝转移患者中5位满足至少20个上皮细胞；GSE231559的8个肝转移患者标签中7个满足门槛，并合并L8的两份肿瘤样本。后者未假定C和L标签相互配对，患者标签独立性继续保留其元数据边界。</p>
<p>合格肝转移群体的actomyosin与CGNL1未调整Spearman相关分别为+0.60（n=5）与+0.18（n=7），增殖相关分别为−0.30与+0.27。它们没有提供一致负方向的补充；没有对小样本强行拟合完整多因素增量模型，也没有合并成一个未经资格确认的“12患者验证队列”。</p>'''+table(liver)+'''</section>
<section><h2>方法</h2>
<p>冻结双终点、基因面板、模型与统计规则后获取新的目标结果。新原发队列采用GEO公开原始UMI计数，逐基因流式读取并验证非负整数性、矩阵列名与元数据细胞的一一对应。患者内上皮原始计数合计，再计算ln(1+CP10k)。CGNL1零表达细胞保留；基因面板不含CGNL1。actomyosin固定为MYH9、MYH10、MYL9、MYL12A、MYL12B、ROCK1、ROCK2；增殖固定为MKI67、TOP2A、PCNA、UBE2C、CENPF。</p>
<p>细胞QC为UMI≥500、检测基因≥200、线粒体比例≤25%；每患者至少20个合格上皮细胞，另保留50细胞门槛。GSE132465使用作者上皮标签并加QC，固定谱系标记交集为敏感性。GSE200997公开注释提供CMS而非细胞类型，因此在新目标关联分析前记录修订，以固定上皮标记最高分、至少2个标记检出及0.4分差为主定义，0.8分差为身份敏感性。两种规则均排除CGNL1和目标程序基因，输出逐细胞QC及完整患者覆盖。</p>
<p>M0控制平均log1p细胞UMI深度；M1另控制EPCAM/KRT8/KRT18/KRT19身份分数及CA1/CA2/GUCA2A/SLC26A3/AQP8/CEACAM7/CLCA4成熟分化分数。六普通标志为CDH1、KRT20、CLDN3、CLDN4、VIL1、CDX2。目标推断使用HC3稳健标准误与t分布，患者数≥10、满秩、残差自由度≥5。原始P&lt;0.05用于名义筛选，BH-FDR按每个队列/模型/细胞定义/门槛内的两个目标终点计算；完整检验家族保留。小肝转移队列本轮仅描述方向和效应。</p>
<p>留一患者时，标准化、常数列判断及拟合均仅使用训练患者；配对基线和增加CGNL1的模型在相同留出患者比较。旧发现队列迁移使用统一患者合计计数的原全基因分母版本，新队列互相迁移亦仅依赖训练来源。未调参或按预测结果重选基因。</p>
<figure><img src="figures/03_coverage.png"><figcaption>图3：每个新队列所有肿瘤患者的上皮覆盖，包括未达到20细胞门槛者。患者数与细胞数分别记录。</figcaption></figure></section>
<section><h2>讨论与下一步</h2>
<p>原23患者中已经建立的上皮状态增量仍然保留。本轮评价其能否跨独立患者来源推广，而不是反向否定原队列事实。M0中的增殖关联、M1条件关联和患者外增量是不同证据；前者成立不能替代后两者。原发CRC与肝转移在组织环境、治疗、平台及肿瘤状态上不同，现有结果不能确定哪一因素造成差异。作者或固定标记定义的上皮群体也不等同于经过独立DNA验证的恶性克隆。</p>
<p>两套新队列与旧来源在研究登记、机构与样本标签上分开；未做跨队列基因型去重。表达程序是RNA代理，不代表肌球蛋白定位、张力、实际增殖率或肝定植因果效应。筛选依据已固定，但整条研究仍属于既有结果驱动的探索性深化。</p>
<p>后续优先核实SCP2079的真实上皮配对数、取样部位和归一化单位，使用患者配对变化检验CGNL1与增殖/骨架状态联系；E-MTAB-12022和GSE158692分别用于病理生长模式和上皮来源明确的方向检查。需要先完成这些数据资格核验，再冻结对应方法；当前没有将这些候选标成已经通过CGNL1验证。</p></section>
<section><h2>内部完整性验收与可复核文件</h2>
<p>GSE200997实际矩阵与注释均为49,859细胞；GEO摘要中的49,589文字值不作为分析计数。文件长度、gzip完整性、跨机SHA256均验证；所有系数、P、FDR、评分与主模型患者外预测使用独立实现复核。六个未成功访问的图谱加载入口、失败传输及修复日志均保留。</p>'''
links=[('冻结主方案','SPECIFICATION.json'),('注释字段修订','ANNOTATION_AMENDMENT.json'),('下载与跨机校验','DOWNLOAD_MANIFEST.json'),('独立数值验收','NUMERICAL_QA.json'),('完整目标统计表','results/N1_CGNL1_models.tsv'),('全部模型资格','results/N1_model_eligibility.tsv'),('逐患者预测','results/N1_heldout_predictions.tsv'),('全部预测增量','results/N1_predictive_increment.tsv'),('普通标志物对照','results/N2_benchmark_increment.tsv'),('跨队列预测','results/N3_transfer_performance.tsv'),('迁移误差分解','results/N3_transfer_error_decomposition.tsv'),('患者覆盖','results/primary_cell_coverage.tsv'),('肝转移描述性结果','results/N4_liver_descriptive.tsv'),('配对部位完整效应','results/N4_paired_site_descriptive.tsv')]
body+='<ul>'+''.join('<li><a href="'+p+'">'+label+'</a></li>' for label,p in links)+'</ul><p>主要原始来源：<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997">GSE200997</a>；<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132465">GSE132465</a>；<a href="https://www.nature.com/articles/s41591-022-02181-8">Tian配对治疗研究</a>；<a href="https://link.springer.com/article/10.1186/s12943-023-01713-1">生长模式队列</a>；<a href="https://crc.i-med.ac.at/">CRC图谱</a>。</p></section>'
css='''body{margin:0;background:#f1f5f6;color:#19333e;font:16px/1.85 "Microsoft YaHei","Noto Sans CJK SC",sans-serif}main{max-width:1140px;margin:auto;padding:28px 22px 55px}header{background:#163f4d;color:white;border-radius:12px;padding:30px}.tag{font-size:14px;color:#b1dfdd}h1{font-size:29px;line-height:1.5}.lead{font-size:18px}section{background:white;border:1px solid #dce5e9;border-radius:9px;padding:25px;margin-top:22px}h2{font-size:22px;border-left:5px solid #137e87;padding-left:13px;line-height:1.5}h3{font-size:18px}table{display:block;overflow:auto;width:100%;font-size:12px;border-collapse:collapse}th,td{border:1px solid #dce5e9;padding:7px;text-align:left}th{background:#e8f2f4}figure{margin:24px 0}img{max-width:100%;height:auto}figcaption{font-size:14px;color:#58727d}a{color:#087f8c}@media(max-width:720px){main{padding:12px}header,section{padding:16px}h1{font-size:23px}}'''
doc='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CGNL1 新队列搜查与验证</title><style>'+css+'</style></head><body><main>'+body+'</main></body></html>'
(R/'CGNL1_new_cohorts_report.html').write_text(doc,encoding='utf-8')
for link in re.findall(r'(?:href|src)="([^"]+)"',doc):
 if not link.startswith('http'):assert (R/link).is_file(),link
decision=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),target='CGNL1',new_primary_cohorts=counts,
 atlas_loaders_attempted=atlas['loaders'],atlas_loaders_read=atlas_read,candidate_registry_entries=len(registry['candidates']),
 primary_M1_results=main[main.model.eq('M1')][['dataset','outcome','n','beta','P','FDR']].to_dict('records'),
 primary_M1_increment=mp[mp.model.eq('M1')][['dataset','outcome','n','improvement','patients_improved']].to_dict('records'),
 interpretation='Positive bidirectional cycling transfer includes centered-error improvement; within-cohort conditional inference and LOPO do not establish stable independence. Retain exploratory predictive clue.',new_primary_to_primary_cycling_transfer=decomp[(decomp.train!='old_discovery23')&decomp.outcome.eq('cycling')].to_dict('records'),liver_small_cohort_descriptive_only=True,next_priority='SCP2079: qualify epithelial paired donors,biopsy sites and normalized expression scale before freezing analysis',
 this_search_and_validation_round_complete=True,whole_CGNL1_project_complete=False,results_directory=str(R))
(R/'INVESTIGATION_DECISION.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2))
(R/'REPORT_VALIDATION.json').write_text(json.dumps(dict(figures=3,relative_links_valid=True,numerical_QA=True),indent=2))
print('REPORT_COMPLETE',flush=True)
