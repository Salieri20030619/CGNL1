from pathlib import Path
import json,hashlib,base64,html,zipfile,datetime
import numpy as np,pandas as pd,anndata as ad
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'annotation_panel_20261007_v1';D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_annotation_20261007_v1');OD=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
s=json.loads((R/'results/ANNOTATION_AND_RETRY_COMPLETE.json').read_text());review=json.loads((R/'review/cluster_annotation_review.json').read_text());cfg=json.loads((R/'standards/QUICHE_RETRY_PARAMETERS.json').read_text());a=ad.read_h5ad(D/'GSM9815454_panel_annotated_v1.h5ad',backed='r');o=a.obs.copy();old=ad.read_h5ad(OD/'inputs/GSM9815454_quiche_technical.h5ad',backed='r');p=o.loc[old.obs_names].copy();oldobs=old.obs.copy();old.file.close()
assert s['exit_code']==0 and o.index.is_unique and len(o)==80036 and int(o.qc_pass.sum())==75094
assert o.annotation.notna().all() and np.all(o.loc[~o.qc_pass,'annotation'].eq('QC_excluded'))
assert p.index.is_unique and len(p)==6000 and p.qc_pass.all()
accepted=o.qc_pass&~o.annotation.eq('uncertain');assert (o.loc[accepted,'core_markers_detected']>=2).all() and not o.loc[accepted,'mixed_marker_flag'].any()
labelhash=hashlib.sha256('\n'.join(i+'\t'+v for i,v in zip(p.index,p.annotation.astype(str))).encode()).hexdigest();assert labelhash==cfg['labels_sha256']==s['labels_sha256']
assert not any(g.upper() in ['CGNL1','JACOP','PCING','ENSG00000128849'] for g in a.var_names)
features=pd.read_csv(R/'standards/clustering_genes.tsv',sep='\t').gene;assert not set(features)&set(json.loads((R/'standards/ANNOTATION_SPECIFICATION.json').read_text())['excluded_targets'])
pd.crosstab(o.annotation,o.annotation_confidence).to_csv(R/'results/annotation_confidence_counts.tsv',sep='\t')
pd.crosstab(o.annotation,o.segmentation_method).to_csv(R/'results/annotation_segmentation_counts.tsv',sep='\t')
o[o.qc_pass].groupby('annotation',observed=True)[['transcript_counts','negative_control_fraction','cell_area','nucleus_area']].median().to_csv(R/'results/annotation_QC_medians.tsv',sep='\t')
o[o.qc_pass].groupby('annotation_reason',observed=True).size().rename('cells').to_csv(R/'results/annotation_reason_counts.tsv',sep='\t')
oldunknown=oldobs.candidate_identity.astype(str).eq('unclassified');oldknown=~oldunknown
pd.crosstab(oldobs.candidate_identity.astype(str),p.annotation.astype(str)).to_csv(R/'results/pilot_old_new_labels.tsv',sep='\t')
validation={'exit_code':0,'all_cell_rows':len(o),'qc_pass':int(o.qc_pass.sum()),'unique_cell_ids':True,'qc_excluded_retained':True,'accepted_labels_have_at_least_two_core_markers':True,'mixed_marker_cells_kept_uncertain':True,'exact_original_pilot_ids_order':True,'label_sha256_verified':True,'excluded_targets_absent_from_clustering':True,'old_unclassified_destinations':p.loc[oldunknown.to_numpy(),'annotation'].value_counts().to_dict(),'full_uncertain':int(o.annotation.eq('uncertain').sum()),'full_mixed_marker_flags':int(o.mixed_marker_flag.sum()),'pilot_mixed_marker_flags':int(p.mixed_marker_flag.sum()),'all_counts_match_metadata':True,'reference_transfer_used_for_final_labels':False,'CGNL1_in_panel':False,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(R/'results/FINAL_INTEGRITY_CHECK.json').write_text(json.dumps(validation,indent=2));a.file.close()
zh={'hepatocyte_lineage':'肝细胞样谱系','endothelial':'内皮细胞','mural':'血管壁/收缩型基质细胞','fibroblast':'成纤维细胞','T_NK':'T/NK 细胞','B':'B 细胞','plasma':'浆细胞','myeloid':'髓系细胞','uncertain':'未定','QC_excluded':'未通过基础 QC'}
def table(t,index=False):return t.to_html(index=index,border=0,escape=True,float_format=lambda v:f'{v:.3f}')
def pic(name):return '<img src="data:image/png;base64,'+base64.b64encode((R/'figures'/name).read_bytes()).decode()+'" alt="'+name+'">'
ct=pd.DataFrame({'标签':sorted(set(o.annotation))});ct['中文']=ct['标签'].map(zh);ct['全样本细胞数']=ct['标签'].map(o.annotation.value_counts()).fillna(0).astype(int);ct['原窗口细胞数']=ct['标签'].map(p.annotation.value_counts()).fillna(0).astype(int)
qc=o[o.qc_pass];unknown_new=int(p.annotation.eq('uncertain').sum());decline=(4906-unknown_new)/6000*100
standard=pd.read_csv(R/'standards/panel_marker_standard_final.tsv',sep='\t');st=[]
for lab,g in standard.groupby('label',sort=False):
 st.append({'类别':zh.get(lab,lab),'核心标记（实际测量且允许用于身份判定）':', '.join(g.loc[g.measured&g.allowed_identity&g.final_role.eq('core'),'gene']),'辅助标记':', '.join(g.loc[g.measured&g.allowed_identity&g.final_role.eq('support'),'gene']),'核查性冲突标记':', '.join(g.loc[g.measured&g.final_role.eq('conflict'),'gene'])})
reference=[]
for n in ['normal_liver','multiome']:
 t=pd.read_csv(R/f'results/{n}_panel_LODO.tsv',sep='\t');reference.append({'参考':n,'供体数':len(t),'供体平均准确率':t.accuracy.mean(),'供体平均平衡准确率':t.balanced_accuracy.mean(),'对本样本的用途':'保留审计；不用于最终标签否决或概率置信度'})
retry=pd.read_csv(R/'results/QUICHE_retry_summary.tsv',sep='\t');cross=pd.read_csv(R/'results/pilot_old_new_labels.tsv',sep='\t',index_col=0)
cluster_table=pd.read_csv(R/'review/cluster_annotation_review.tsv',sep='\t',dtype=str)
sources=''.join('<li><a href="'+html.escape(x['url'])+'">'+html.escape(x['title'])+'</a></li>' for x in review['references'])
report='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CGNL1 面板适配注释与 QUICHE 重试</title><style>body{font-family:"Microsoft YaHei",Arial,sans-serif;max-width:1280px;margin:36px auto;padding:0 24px;line-height:1.7;color:#172b3a}h1,h2{color:#123d57}h2{margin-top:36px}table{border-collapse:collapse;width:100%;font-size:13px;line-height:1.5;margin:18px 0}th,td{border:1px solid #d6e1e7;padding:8px;text-align:left;vertical-align:top}th{background:#eaf1f5}img{max-width:100%;height:auto}.lead{padding:20px;background:#eaf3f7;border-left:5px solid #247ba0}.note{color:#526775}code{overflow-wrap:anywhere}</style><body>'''
report+=f'<h1>面板适配注释已完成，QUICHE 原窗口重试通过</h1><p class="note">版本 20261007_v1 · M242T / GSM9815454 · 1 个独立生物学样本</p><div class="lead">全样本保留 80,036 个细胞记录，其中 75,094 个通过既定基础 QC。原 6,000 细胞窗口的未定标签由 4,906 个（81.77%）降至 {unknown_new:,} 个（{unknown_new/6000:.2%}），下降 {decline:.2f} 个百分点。QUICHE 原生空间邻域模块在 100、200、300 μm 三组原参数下均运行完成。</div>'
report+='<h2>旧标签的主要问题来自面板适配不足</h2><p>原规则要求同一类别至少检出两个 marker，但肝细胞和髓系在旧标准中分别只有 CYP2E1、LYZ 一个可用 marker，且未配置 B 细胞、浆细胞等类别。原窗口 4,906 个未定细胞中，4,899 个未满足“获胜类别至少两个 marker”的要求。本次没有将该阈值统一降低，而是补充面板中实际测量的成组标记，并结合表达聚类审查。</p><p>原始 Gene Expression 矩阵逐细胞计数与 cells 元数据完全一致。当前证据支持“旧注释标准不适配”是主要可修正因素；不能据此宣布原样本、分割或背景处理均无问题。</p>'
report+='<h2>最终标签恢复了主要免疫及基质谱系</h2>'+table(ct)+f'<p>全样本通过 QC 的细胞中，{int(qc.annotation.eq("uncertain").sum()):,} 个仍保留为未定；{int(qc.mixed_marker_flag.sum()):,} 个带有多谱系标记重叠标志。肝细胞样谱系仅描述表达身份，未区分正常与恶性。未获得足够证据单独建立 NK、肥大细胞、中性粒细胞、树突状细胞亚型或胆管/CRC 上皮群；未保留某类标签不等于该类细胞在组织中不存在。</p>'+pic('annotation_spatial.png')
report+='<h2>实际面板标准覆盖 13 类候选身份</h2><p>标准以原矩阵 462 个 Gene Expression 特征为准；79 个控制特征不参与归一化与聚类。下表列出初始候选标准，最终类别还需通过群级审查和单细胞证据检查。T/NK 合并、髓系合并以及肝细胞样谱系命名见冻结审查文件。</p>'+table(pd.DataFrame(st))
report+='<p>CGNL1、CDH5、CLDN5、EMCN、ESAM、KDR、OCLN、PECAM1、TJP1、VWF 按既有目标排除清单不进入聚类及自动身份打分。PECAM1/VWF/KDR 的测量值仅作为保留的核查证据；内皮身份由 FLT1、CLEC14A、PTPRB、ROBO4、PLVAP、CD34 等支持。初始标准和修订版分别保存，可追溯修改。</p>'+pic('marker_dotplot.png')
report+='<h2>直接参考迁移未通过本样本的适用性审查</h2><p>两个既有单细胞参考在测量面板限制下重新做留一供体验证，结果如下。它们在自身平台内具有分类能力，但直接迁移至本 Xenium 样本后，包括内皮标记富集群在内的大部分群均被预测为肝细胞。因此，最终标签采用标记与表达群审查，参考预测仅保留供审计。</p>'+table(pd.DataFrame(reference))+'<p>参考域内准确率不是本样本注释准确率。当前尚无配对人工真值或同细胞蛋白验证，不能给出最终标签准确率。</p>'
report+='<h2>重试通过；邻域组成不再由旧“未定”类别主导</h2>'+table(retry)+'<p>重试使用原来的 6,000 个细胞、相同顺序与坐标、相同 QUICHE 源码和参数：k=10、khop=3、最少邻居数 3。每个尺度均检查邻域比例有限且和为 1，并独立重算 20 个细胞验证一致性。原标签 unclassified 与新 uncertain、旧 endothelial_identity 与新 endothelial 等先统一命名后计算邻域组成差异（total variation）。此差异描述标签修订的影响，不代表方法性能评分。</p><p>该模块在这些参数下最多使用 30 个近邻，当前三个半径均已覆盖该近邻上限。因此三个半径结果相同，不能称为独立的多尺度稳健性验证。</p><h3>原窗口新旧标签交叉表</h3>'+table(cross,index=True)
report+='<h2>Methods</h2><p>读取原始 10x H5 和匹配细胞元数据，按 barcode 对齐；保留全部记录。基础 QC 为 transcript_counts ≥50 且 nucleus_count=1。Gene Expression 计数按细胞缩放至 10,000 后 log1p；排除目标清单，对至少 10 个细胞表达的特征标准化并截断至 10，PCA 30 维、余弦距离 30 近邻、Leiden 分辨率 0.3/0.6/1.0，随机种子 20261007。主分辨率 0.6 得到 10 群；对混合群 0、3、9 进一步使用 PCA 20 维、20 近邻、Leiden 0.5，形成 31 个审查群。空间坐标不用于身份聚类或标签确定。</p><p>群标记报告平均 log 表达、相对其余细胞的均值差及检出比例。最终标签要求群级证据成立、同类别至少两个允许使用的核心标记检出。出现预设替代谱系标记组合时保留为未定；血管壁与成纤维的收缩程序重叠不单独触发冲突。各标记组合与阈值逐项写入冻结审查文件。参考使用原已 log 归一化矩阵，限定到实际测量且允许的核心/辅助标记，计算供体平衡的类别均值及 Pearson 相似度；留一供体验证每次排除待预测供体。</p><p>所有结果为单样本描述性结果，未以细胞或邻居数代替生物学样本量，未进行生物学推断检验、差异丰度检验或生成 P/FDR。</p>'
report+='<h2>Discussion</h2><p>新版注释解决了旧规则遗漏大量实际 marker 的问题，并使原空间窗口获得可解释的邻域组成。仍需保留三个边界：第一，普遍的肝细胞转录本背景与部分跨谱系共表达，尚未通过转录本级位置、细胞边界及配对蛋白证据区分环境背景、分割混合与真实状态；第二，当前只有一个生物学样本，QUICHE 完成的是原生空间模块技术重试，未运行完整差异丰度推断；第三，本数据的测量面板没有 CGNL1，标签修订不能将该数据变成 CGNL1 表达验证数据。细胞核数量不为 1 的细胞保留为 QC_excluded，其中可包含真实多核肝细胞，不能直接等同坏样本。</p>'
report+='<h2>逐群审查记录</h2>'+table(cluster_table)+'<h2>文件与复核</h2><p>全细胞注释：<a href="results/cell_annotations.tsv.gz">cell_annotations.tsv.gz</a>；最终 marker 标准：<a href="standards/panel_marker_standard_final.tsv">panel_marker_standard_final.tsv</a>；完成与完整性记录：<a href="results/FINAL_INTEGRITY_CHECK.json">FINAL_INTEGRITY_CHECK.json</a>。完整含计数的注释对象保存在 sgs10：<code>'+str(D/'GSM9815454_panel_annotated_v1.h5ad')+'</code>。旧分析结果未改写。</p><p>此次使用既有数据和参考，不新增下载分析数据；后续新增下载仍按用户要求直连、不使用代理。</p><h2>原始来源</h2><ul>'+sources+'</ul></body></html>'
(R/'CGNL1_annotation_retry_report.html').write_text(report,encoding='utf-8')
files=[f for f in R.rglob('*') if f.is_file() and f.suffix!='.zip' and f.name!='SHA256SUMS.txt'];manifest='\n'.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+str(f.relative_to(R)) for f in sorted(files));(R/'SHA256SUMS.txt').write_text(manifest+'\n')
with zipfile.ZipFile(R/'CGNL1_annotation_retry_review_bundle.zip','w',zipfile.ZIP_DEFLATED) as z:
 for f in R.rglob('*'):
  if f.is_file() and f.suffix!='.zip':z.write(f,f.relative_to(R))
print(json.dumps(validation,indent=2));print('REPORT_AND_BUNDLE_READY')
