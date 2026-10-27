from pathlib import Path
import json,base64,hashlib,zipfile,datetime
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/multiscale_20261007_v1');A=R.parent/'annotation_panel_20261007_v1'
def read(name):return pd.read_csv(R/'results'/name,sep='\t')
s=read('scale_summary.tsv');pairs=read('scale_comparisons.tsv');edge=read('window_edge_sensitivity.tsv');corr=read('composition_correlations.tsv');focal=read('focal_lineage_composition.tsv');cell=pd.read_csv(R/'results/per_cell_multiscale_metrics.tsv.gz',sep='\t',index_col=0);complete=json.loads((R/'results/MULTISCALE_COMPLETE.json').read_text());assert complete['exit_code']==0
order=['context_knn_15','context_knn_30','context_knn_60','context_radius_100','context_radius_200','context_radius_300'];short=['k=15','k=30','k=60','100 um','200 um','300 um'];palette=['#ce813e']*3+['#247ba0']*3
profiles={c:pd.read_csv(R/f'results/{c}_composition.tsv.gz',sep='\t',index_col=0) for c in s.condition}
focal['matched_all_focal_neighborhood_mean']=[profiles[r.condition][r.neighbor_label].mean() for r in focal.itertuples()];focal['difference_from_matched_neighborhood_pp']=100*(focal.mean_fraction-focal.matched_all_focal_neighborhood_mean);focal.to_csv(R/'results/focal_lineage_composition_with_matched_background.tsv',sep='\t',index=False)
vascular=focal[focal.condition.isin(order)&focal.focal_label.isin(['endothelial','mural'])&focal.neighbor_label.isin(['endothelial','mural'])].copy()
fig,ax=plt.subplots(2,2,figsize=(13,9));ctx=s.set_index('condition').loc[order]
for j,(metric,title,unit) in enumerate([('neighbors','Neighborhood size','Other cells'),('distance','Physical extent','Distance to furthest neighbor (um)')]):
 med=ctx[metric+'_median'+('_um' if metric=='distance' else '')].to_numpy();lo=ctx[metric+'_q25'+('_um' if metric=='distance' else '')].to_numpy();hi=ctx[metric+'_q75'+('_um' if metric=='distance' else '')].to_numpy();ax[0,j].bar(short,med,color=palette);ax[0,j].errorbar(short,med,yerr=[med-lo,hi-med],fmt='none',ecolor='#222222',capsize=4);ax[0,j].set_title(title+' (median and IQR)');ax[0,j].set_ylabel(unit)
compkeys=['context_knn_15_vs_30_TV','context_knn_30_vs_60_TV','context_knn_15_vs_60_TV','context_radius_100_vs_200_TV','context_radius_200_vs_300_TV','context_radius_100_vs_300_TV'];tick=['15-30','30-60','15-60','100-200','200-300','100-300'];box=ax[1,0].boxplot([cell[c].to_numpy()*100 for c in compkeys],tick_labels=tick,showfliers=False,patch_artist=True)
for p,c in zip(box['boxes'],palette):p.set_facecolor(c)
ax[1,0].set_ylabel('Total variation x 100');ax[1,0].set_title('Composition changes across scales');ax[1,0].set_xlabel('Orange: k pairs; blue: radius pairs (um)')
for gene,color in [('endothelial','#247ba0'),('mural','#7859a6')]:ax[1,1].plot(short,ctx['endo_fraction_sd' if gene=='endothelial' else 'mural_fraction_sd']*100,marker='o',label=gene,color=color)
ax[1,1].set_ylabel('SD of local cell fraction (percentage points)');ax[1,1].set_title('Larger neighborhoods average out local variation');ax[1,1].legend();fig.suptitle('Same 6,000 focal cells; neighbors from all 75,094 QC cells',fontsize=14);fig.tight_layout();fig.savefig(R/'figures/scale_comparison.png',dpi=170);plt.close(fig)
vm=[]
for condition in order:
 g=vascular[vascular.condition.eq(condition)].set_index(['focal_label','neighbor_label']);vm.append([g.loc[(a,b),'difference_from_matched_neighborhood_pp'] for a,b in [('endothelial','endothelial'),('endothelial','mural'),('mural','endothelial'),('mural','mural')]])
vm=np.array(vm);fig,ax=plt.subplots(figsize=(8,5));im=ax.imshow(vm,cmap='RdBu_r',vmin=-4.5,vmax=4.5,aspect='auto');ax.set_yticks(range(6),short);ax.set_xticks(range(4),['Endo -> Endo','Endo -> Mural','Mural -> Endo','Mural -> Mural']);ax.set_title('Local vascular composition relative to all focal neighborhoods')
for i in range(6):
 for j in range(4):ax.text(j,i,f'{vm[i,j]:+.2f}',ha='center',va='center',color='white' if abs(vm[i,j])>3 else '#172b3a')
fig.colorbar(im,ax=ax,label='Difference in percentage points');fig.tight_layout();fig.savefig(R/'figures/vascular_scale_comparison.png',dpi=180);plt.close(fig)
pilot=pd.read_csv(A/'results/pilot6000_annotations.tsv.gz',sep='\t',index_col=0).loc[cell.index];xy=pilot[['x_centroid','y_centroid']].to_numpy();fig,axes=plt.subplots(2,3,figsize=(12,9))
for ax,c,title in zip(axes.ravel(),order,short):
 im=ax.scatter(*xy.T,c=profiles[c].loc[cell.index,'endothelial'],s=2,cmap='viridis',vmin=0,vmax=.3,rasterized=True);ax.set_title(title);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([])
fig.suptitle('Endothelial neighbor fraction across scales (same color range)');fig.colorbar(im,ax=axes.ravel().tolist(),label='Endothelial fraction',shrink=.6);fig.savefig(R/'figures/endothelial_spatial_multiscale.png',dpi=170);plt.close(fig)
def table(t,index=False):return t.to_html(index=index,border=0,float_format=lambda x:f'{x:.4f}')
def pic(n):return '<img alt="'+n+'" src="data:image/png;base64,'+base64.b64encode((R/'figures'/n).read_bytes()).decode()+'">'
main=ctx[['neighbors_median','neighbors_q25','neighbors_q75','distance_median_um']].reset_index();main.columns=['条件','邻居数中位数','邻居数Q25','邻居数Q75','最远邻居距离中位数(μm)'];main['条件']=short
contextpairs=pairs[pairs.pool.eq('context')].copy();contextpairs[['mean_TV','median_TV','p90_TV']]*=100;contextpairs=contextpairs.rename(columns={'mean_TV':'mean_TV_x100','median_TV':'median_TV_x100','p90_TV':'p90_TV_x100'})
vp=vascular.copy();vp['邻居比例(%)']=100*vp.mean_fraction;vp['相同中心集合背景(%)']=100*vp.matched_all_focal_neighborhood_mean
report='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CGNL1 双路径多尺度邻域分析</title><style>body{font-family:"Microsoft YaHei",Arial,sans-serif;max-width:1200px;margin:36px auto;padding:0 24px;color:#173043;line-height:1.75}h1,h2{color:#154b69}h2{margin-top:36px}table{border-collapse:collapse;width:100%;font-size:13px;margin:18px 0;line-height:1.5}td,th{border:1px solid #d4e0e7;padding:8px;text-align:left;vertical-align:top}th{background:#edf4f7}img{max-width:100%;height:auto}.lead{background:#edf5f8;padding:20px;border-left:5px solid #247ba0}.note{color:#516777}code{overflow-wrap:anywhere}</style><body><h1>两类邻域均运行成功，真实尺度变化已得到验证</h1><p class="note">20261007_v1 · M242T/GSM9815454 · 同一批 6,000 个中心细胞、冻结注释、1 个独立生物学样本</p>'''
report+='<div class="lead">取消近邻数量上限后，100/200/300 μm 对应的邻居数中位数为 173/687/1,549；k=15/30/60 对应约 29.4/41.5/58.8 μm。两种方案均产生实际不同的邻域。较大邻域使局部组成更加平滑，但这不构成更高准确性的证据。</div>'
report+='<h2>同一中心集合下，两种定义覆盖不同的物理范围</h2><p>原窗口内进行 6 组 QUICHE 原生模块运行，并以全切片 75,094 个 QC 细胞提供邻居，再计算同样 6 组窗口边缘敏感性对照。两组均固定原来的 6,000 个中心、细胞身份和坐标。下表优先展示补足窗口外邻居的结果；原窗口结果完整保存在结果包。</p>'+table(main)+pic('scale_comparison.png')
report+='<h2>组成随尺度改变；200–300 μm 的相似主要体现更强的空间平均</h2><p>补足窗口外邻居后，100→200、100→300、200→300 μm 的平均总变差分别为 0.0573、0.0651、0.0270；15→30、15→60、30→60 近邻分别为 0.1122、0.1574、0.0892。总变差定义为各类别比例绝对差之和的一半，范围 0–1；表格为便于阅读乘以 100。它是组成变化量，不是准确率或显著性。</p><p>内皮邻域比例的跨细胞标准差从 k=15 的 10.15 个百分点降至 k=60 的 5.55 个百分点，并在 100/200/300 μm 下继续降至 3.72/2.17/1.29 个百分点。大邻域稀释局部差异，不能仅凭其更加平稳选择“最佳参数”。</p>'+table(contextpairs)+pic('endothelial_spatial_multiscale.png')
report+='<h2>小尺度存在内皮同类聚集线索，未见稳定的内皮—血管壁共同富集</h2><p>以同一组 6,000 个中心在相同参数下的平均邻域组成为比较背景，内皮中心周围的内皮比例在 k=15/30/60 下分别高出 4.27/2.20/1.58 个百分点；在 100/200/300 μm 下为 +0.78/+0.37/+0.22 个百分点。探索性提示内皮同类空间聚集主要在较小尺度可见。</p><p>内皮中心周围的血管壁细胞比例与上述背景的差值在六种参数下为 −0.09 至 +0.04 个百分点，没有一致的正向变化；血管壁中心周围的内皮比例也随尺度改变方向。当前数据不支持稳定的内皮—血管壁共同富集。这里的“周围”指质心距离邻域，不能直接解释为细胞膜接触。</p>'+pic('vascular_scale_comparison.png')+table(vp[['condition','focal_label','focal_cells','neighbor_label','邻居比例(%)','相同中心集合背景(%)','difference_from_matched_neighborhood_pp']])
report+='<h2>窗口边界影响大半径邻居数量，但没有解释掉尺度变化</h2><p>100/200/300 μm 下，29.47%/55.17%/75.27% 的中心原本遗漏窗口外邻居。300 μm 时，每个中心平均补入 315 个邻居；补足前后组成的平均总变差为 0.0150。k=15/30/60 的组成受窗口补足影响的中心比例为 6.58%/11.10%/16.40%。固定 k 的邻居数不变，但边缘细胞的具体邻居可以替换。</p><p>另取距所有窗口外 QC 细胞至少 300 μm 的 1,484 个内部中心：窗口内与全切片背景的半径结果完全一致，100→300 μm 的平均总变差仍为 0.0548。说明尺度差异在不受窗口截断影响的区域仍存在。</p>'+table(edge)
report+='<h2>Methods</h2><p>读取既有注释对象和原窗口 cell ID，验证冻结标签 SHA256，保留 uncertain 类别参与分母；所有邻域排除中心自身。固定半径通过未修改 QUICHE graph.py 的 compute_spatial_neighbors，设置 n_neighbors=None、delaunay=False，分别构造 100/200/300 μm 邻域。固定数量通过 spatial_niches_khop 设置 k=15/30/60、khop=1、radius=∞。按实际质心坐标计算欧氏距离，最少邻居数设为 3。全切片背景对照用 cKDTree 构造相同定义的邻域，仅改变可用邻居池，中心集合不变。</p><p>各邻域以等权细胞数计算类型比例。记录邻居数、最大距离、比例标准差、尺度间总变差和类别比例 Spearman 相关。按中心身份汇总邻居比例，并同时保留全邻居池背景及相同中心集合的邻域平均背景；正文使用后者描述局部组成差异。后者是描述性比较背景，不是置换零模型，也未调整空间混杂。未进行差异丰度或其他生物学推断检验，没有以细胞数代替生物学样本量，也未生成 P/FDR。</p>'
report+='<h2>后续参数选择与解释边界</h2><p>针对局部血管微环境，可继续以 k=30 为工作尺度，k=15 和 k=60 作为敏感性分析；100 μm 提供固定物理距离的补充视角，200/300 μm 描述更广的组织背景。所有正式窗口分析应允许全切片 QC 细胞提供邻居。这个安排依据本样本显示的物理跨度与平均效应，并不表示已确定生物学最优参数。</p><p>本次完成的是单样本的邻域构建和描述性尺度评估，没有运行完整 QUICHE 差异丰度推断。注释、分割和跨谱系转录本背景仍沿用上一版已记录的适用边界。当前面板不包含 CGNL1；本结果不能作为 CGNL1 机制验证。窗口外邻居补足只消除人为裁剪边界，不能补回切片之外的组织。</p>'
report+='<h2>验证与完整结果</h2><p>12 个条件均保留全部 6,000 个中心；检查比例有限且和为 1、无自身邻居，并各抽取固定随机种子的 20 个中心，以全量欧氏距离独立重算核对。原窗口 k=30 与上次注释后 QUICHE 结果逐项一致。完成标记、退出状态和文件校验均保留。</p><ul><li><a href="results/scale_summary.tsv">全部尺度参数和邻域大小</a></li><li><a href="results/focal_lineage_composition_with_matched_background.tsv">完整中心身份—邻居组成表</a></li><li><a href="results/composition_correlations.tsv">各类别尺度间相关和比例变化</a></li><li><a href="results/window_edge_sensitivity.tsv">窗口边界敏感性</a></li><li><a href="results/per_cell_multiscale_metrics.tsv.gz">逐细胞指标</a></li><li><a href="standards/MULTISCALE_SPECIFICATION.json">冻结参数</a></li><li><a href="results/NUMERICAL_VALIDATION.json">数值核查</a></li></ul><p>邻接矩阵另存于 sgs10：<code>/mnt/BioData/bio/wangmingyuan/CGNL1_multiscale_20261007_v1/graphs/</code>。结果包含 12 份组成矩阵、全部统计表、图件和运行脚本；旧结果未改写。本次未新增下载外部分析数据。</p></body></html>'
(R/'CGNL1_multiscale_report.html').write_text(report,encoding='utf-8')
files=[p for p in R.rglob('*') if p.is_file() and p.suffix!='.zip' and p.name!='SHA256SUMS.txt'];(R/'SHA256SUMS.txt').write_text('\n'.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(R)) for p in sorted(files))+'\n')
with zipfile.ZipFile(R/'CGNL1_multiscale_results.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in R.rglob('*'):
  if p.is_file() and p.suffix!='.zip':z.write(p,p.relative_to(R))
print('REPORT_BUNDLE_COMPLETE',flush=True)
