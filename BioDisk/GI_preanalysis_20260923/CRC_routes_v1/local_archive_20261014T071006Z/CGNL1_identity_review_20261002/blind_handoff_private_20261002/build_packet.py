from pathlib import Path
import json, random, hashlib, zipfile, datetime, csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
SOURCE=Path(r'C:\Users\31937\CGNL1_identity_review_20261002')
PRIVATE=SOURCE/'blind_handoff_private_20261002'
DOWNLOAD=Path(r'C:\Users\31937\Downloads')
OUT=DOWNLOAD/'单细胞身份_匿名盲审_20261002_v1'
OUT.mkdir(exist_ok=False)
(OUT/'figures').mkdir()
rng=random.Random(2026100201)
datasets=['GSE245552','GSE285991'];rng.shuffle(datasets)
dsmap={ds:f'D{i+1:02d}' for i,ds in enumerate(datasets)}
cluster_keys=[];data={}
for ds in datasets:
 d=SOURCE/'results'/ds
 o=pd.read_csv(d/'cell_identity_evidence.tsv.gz',sep='\t',dtype={'cluster05':str})
 m=pd.read_csv(d/'cluster_markers_descriptive.tsv',sep='\t',dtype={'cluster':str})
 data[ds]=(o,m)
 cluster_keys.extend((ds,c) for c in sorted(o.cluster05.unique()))
rng.shuffle(cluster_keys);clmap={k:f'C{i+1:03d}' for i,k in enumerate(cluster_keys)}
assert len(clmap)==44
blocked=['CGNL1','CLDN5','ESAM','TJP1','OCLN','GSE245552','GSE285991','GSM','one_marker_removed','two_plus_retained','multiome_prediction','normal_liver_prediction','both_immune','both_endothelial']
overview=[];markers=[];donors=[];private=[];template=[];card=[];donormaps={}
plt.rcParams.update({'font.size':9,'svg.fonttype':'none'})
for ds in datasets:
 o,m=data[ds];dcode=dsmap[ds];donorids=sorted(o.donor.unique());rng.shuffle(donorids);dm={d:f'{dcode}_P{i+1:02d}' for i,d in enumerate(donorids)};donormaps[ds]=dm
 o['anon_cluster']=[clmap[(ds,c)] for c in o.cluster05];o['anon_donor']=o.donor.map(dm)
 origpacket=pd.read_csv(SOURCE/'results'/ds/'cluster_review_packet.tsv',sep='\t',dtype={'cluster':str}).set_index('cluster')
 for cl,g in o.groupby('cluster05'):
  cc=clmap[(ds,cl)];dist=g.donor.value_counts();fr=dist/len(g)
  overview.append(dict(dataset=dcode,cluster=cc,cells=len(g),donors=g.donor.nunique(),median_UMI=float(g.UMI_QC.median()),median_detected_genes=float(g.genes_QC.median()),median_mito_pct=float(g.mito_QC.median()),largest_donor_fraction=float(fr.max())))
  for donor in sorted(dm,key=lambda x:dm[x]):donors.append(dict(dataset=dcode,cluster=cc,donor=dm[donor],cells=int(dist.get(donor,0))))
  gm=m[m.cluster==cl].sort_values('mean_difference',ascending=False)
  assert len(gm)==30
  for rank,row in enumerate(gm.itertuples(),1):markers.append(dict(dataset=dcode,cluster=cc,rank=rank,gene=row.gene,mean_log_expression=float(row.mean_log),mean_difference_vs_other_clusters=float(row.mean_difference),fraction_detected=float(row.fraction_detected),cells=int(row.cells)))
  pred=origpacket.loc[cl]
  private.append(dict(anonymous_dataset=dcode,anonymous_cluster=cc,original_dataset=ds,original_cluster=cl,current_disease_reference=pred.disease_reference,current_normal_reference=pred.normal_reference))
  template.append(dict(dataset=dcode,cluster=cc,broad_identity='',optional_subtype='',supporting_genes='',conflicting_evidence='',confidence='',mixed_or_doublet_concern='',additional_evidence_needed='',reasoning=''))
 # Only anonymous labels are used in figure titles, legends and metadata.
 fig,axs=plt.subplots(1,2,figsize=(15,6))
 for ax,column,title in [(axs[0],'anon_cluster',f'{dcode}: anonymous clusters'),(axs[1],'anon_donor',f'{dcode}: anonymous donors')]:
  names=sorted(o[column].unique());colors=plt.cm.turbo(np.linspace(.05,.95,len(names)))
  for color,name in zip(colors,names):
   v=o[o[column]==name];ax.scatter(v.UMAP1,v.UMAP2,s=.7,color=color,alpha=.65,rasterized=True,label=name)
  ax.set_title(title);ax.set_xlabel('UMAP 1');ax.set_ylabel('UMAP 2');ax.set_xticks([]);ax.set_yticks([]);ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),fontsize=6,markerscale=5,frameon=False)
 fig.tight_layout();fig.savefig(OUT/'figures'/f'{dcode}_anonymous_umap.png',dpi=160,metadata={'Title':f'{dcode} anonymous evidence'});plt.close(fig)
overview=pd.DataFrame(overview).sort_values(['dataset','cluster']);markers=pd.DataFrame(markers).sort_values(['dataset','cluster','rank']);donors=pd.DataFrame(donors).sort_values(['dataset','cluster','donor']);template=pd.DataFrame(template).sort_values(['dataset','cluster'])
overview.to_csv(OUT/'01_cluster_overview.tsv',sep='\t',index=False)
markers.to_csv(OUT/'02_top30_markers.tsv',sep='\t',index=False)
donors.to_csv(OUT/'03_cluster_donor_counts.tsv',sep='\t',index=False)
template.to_csv(OUT/'05_审阅结果模板.tsv',sep='\t',index=False)
for v in overview.itertuples():
 card.extend([f'[{v.dataset} / {v.cluster}]',f'细胞数：{v.cells}；供者数：{v.donors}；最大单一供者占比：{v.largest_donor_fraction:.1%}',f'UMI中位数：{v.median_UMI:.1f}；检出基因数中位数：{v.median_detected_genes:.1f}；线粒体比例中位数：{v.median_mito_pct:.2f}%', '排名\t基因\t簇内平均log表达\t相对其余细胞平均差\t簇内检出比例'])
 for row in markers[(markers.dataset==v.dataset)&(markers.cluster==v.cluster)].itertuples():card.append(f'{row.rank}\t{row.gene}\t{row.mean_log_expression:.4f}\t{row.mean_difference_vs_other_clusters:.4f}\t{row.fraction_detected:.3f}')
 card.append('')
(OUT/'04_逐簇标记证据.txt').write_text('\n'.join(card),encoding='utf-8')
instructions='''单细胞身份注释：匿名盲审说明
版本：2026-10-02 v1

一、交接方式
请在一个未接触本项目既往结论的新对话中上传本证据包，再粘贴“06_可直接粘贴的审阅指令.txt”。不要同时提供项目名称、原始数据集编号、已有细胞标签、自动参考预测、原先筛选分组、目标基因或关联统计结果。若不能读取压缩包，可解压后上传其中的说明、标记表、逐簇证据、供者表和图片。
本包不包含解盲对照。请先保存独立审阅结果，再进行后续解盲比较。

二、审阅任务
这是两组人类组织单细胞RNA表达数据的细胞身份盲审，共44个匿名簇。请依据提供的多基因标记组合、簇内检出比例、表达差异、供者分布与聚类位置，独立判断各簇的主要细胞谱系。研究目的、原始队列和先前答案均被隐藏；不得猜测或反推这些信息。
本轮只判断身份，不评价任何目标基因的功能、关联或显著性，也不选择可支持某条研究结论的细胞。

三、证据来源与处理方式
1. 两组数据分别处理，使用各自全部通过宽松基础QC的细胞，没有先按预期细胞类型筛出某个谱系再聚类。
2. 基础QC：每细胞原始计数总量至少200、检出基因数至少100、线粒体计数比例不超过35%。本轮没有另做双细胞模型剔除，混合簇/疑似双细胞应作为待核查问题记录。
3. 预先指定的一小组研究相关基因已从归一化后的分析对象、聚类特征和标记证据中屏蔽，具体名单不向审阅者揭示。请勿猜测屏蔽名单，亦勿把某基因未列出当作该基因未表达。
4. 屏蔽后的计数按细胞总量归一化至10000后取log1p；按样本选择约3000个高变基因，标准化截断值10，PCA 30维，20邻居，余弦距离，Leiden聚类。此次证据包对应分辨率0.5的主簇。图中展示相同分析的UMAP，仅用于呈现邻域，不能单凭UMAP距离判定身份。
5. 标记表按“本簇平均log表达减去同队列其余细胞平均log表达”排序，每簇列前30个基因；同时给出簇内检出比例。它是描述性标记排序，不是差异表达显著性检验，也没有提供细胞级P值/FDR。
6. 供者、队列和簇编号已匿名化。细胞数不等于独立生物学样本量。最大供者占比高的簇需留意样本驱动，但不能仅因此判为错误簇。
7. 文件不提供任何现有注释、参考预测或原先筛选标签。审阅时请仅使用本包证据及一般细胞生物学知识；如确需查文献，可引用原始研究，但不要检索或反推原始数据集身份。

四、判断要求
1. 每个簇都要给出判断，包括“不确定”或“混合/需要进一步拆分”；不要遗漏小簇。
2. 先给主要谱系，再在证据足够时给亚型。身份与激活/增殖等状态分开，不能以状态词替代主要谱系。
3. 使用多基因组合及相反证据，不以单个标记或固定“至少两个标记”规则自动决定身份。前30名中未出现某个经典标记，不构成可靠阴性证据。
4. 仅凭这些RNA标记，不应把某一上皮簇直接定为恶性，也不应确证双细胞、环境RNA污染或特定功能。
5. 自由选择最符合证据的细胞类型，不预设需要找出哪一类细胞，也不预设各类数量。
6. 提供高/中/低置信度及理由。置信度来自证据一致性，不来自簇大小、显著性或是否符合预期。
7. 如果材料不足以判断，明确需要补充哪些基因的表达、哪些图或原始信息；不要补造未提供的表达、标签、P值或实验结果。

五、返回格式
请使用“05_审阅结果模板.tsv”的全部44行编号，返回一行对应一个簇的表：
匿名队列｜匿名簇｜主要身份｜可选亚型｜支持基因｜相反/缺失证据｜置信度｜混合或双细胞疑点｜需要补充的证据｜简短理由。
最后另列：需优先复核的簇；可合并或需拆分的簇及理由；整体证据局限。
审阅结果完成后，请声明本轮是否接触过原标签、自动预测或既往项目结论。先冻结这份答案，再接收任何解盲材料。

六、文件列表
01_cluster_overview.tsv：44个簇的细胞数、供者数和QC描述。
02_top30_markers.tsv：各簇前30个描述性标记及表达/检出比例。
03_cluster_donor_counts.tsv：匿名簇×匿名供者的细胞数。
04_逐簇标记证据.txt：同一证据的纯文本阅读版。
05_审阅结果模板.tsv：待填写的44行审阅表，未预填任何身份。
06_可直接粘贴的审阅指令.txt：简短启动指令。
figures/：两组匿名聚类及供者UMAP图。
MANIFEST.json：文件完整性校验清单，不含原始身份对照。
'''
(OUT/'00_盲审说明_先读.txt').write_text(instructions,encoding='utf-8')
prompt='''请对我上传的匿名单细胞证据进行独立细胞身份盲审。先阅读“00_盲审说明_先读.txt”，然后逐一审阅全部44个簇，使用“05_审阅结果模板.tsv”的编号返回主要谱系、可选亚型、支持与相反证据、置信度及需要补充的材料。
不要预设目标谱系，不要猜测项目、屏蔽基因或既往结论；不依赖单个标记或固定标记个数作自动分类。标记表仅包含前30名，未列出不代表未表达。身份与细胞状态分开；证据不足时明确标为不确定。
只依据本包的无标签证据作判断，不要求计算关联显著性。先完成并冻结全部簇的独立判断，之后才接收解盲信息。请在结尾说明你本轮是否接触过原标签、自动参考预测或既往项目结论。
'''
(OUT/'06_可直接粘贴的审阅指令.txt').write_text(prompt,encoding='utf-8')
# Private restoration keys and generation code must never enter the upload ZIP.
pd.DataFrame(private).to_csv(PRIVATE/'解盲对照_禁止上传.tsv',sep='\t',index=False)
(PRIVATE/'匿名供者对照_禁止上传.json').write_text(json.dumps(donormaps,ensure_ascii=False,indent=2),encoding='utf-8')
checks=[]
for p in OUT.rglob('*'):
 if not p.is_file():continue
 if p.suffix in ['.txt','.tsv','.json']:
  text=p.read_text(encoding='utf-8')
  leaks=[v for v in blocked if v in text]
  assert not leaks,(p.name,leaks)
 assert not any(v in str(p.relative_to(OUT)) for v in blocked)
 checks.append(dict(path=str(p.relative_to(OUT)).replace('\\','/'),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
assert set(template.cluster)==set(overview.cluster)==set(markers.cluster)==set(donors.cluster)
assert len(markers)==44*30 and len(template)==44
for cc,g in donors.groupby('cluster'):assert int(g.cells.sum())==int(overview[overview.cluster==cc].cells.iloc[0])
(OUT/'MANIFEST.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
zip_path=DOWNLOAD/'单细胞身份_匿名盲审_20261002_v1.zip'
assert not zip_path.exists()
with zipfile.ZipFile(zip_path,'w',compression=zipfile.ZIP_DEFLATED) as z:
 for p in sorted(OUT.rglob('*')):
  if p.is_file():z.write(p,arcname=str(p.relative_to(OUT)))
with zipfile.ZipFile(zip_path) as z:
 assert z.testzip() is None
 assert len(z.namelist())==len(checks)+1
qa={'generated':datetime.datetime.now().astimezone().isoformat(),'clusters':44,'marker_rows':len(markers),'upload_files':len(checks)+1,'source_labels_in_text_or_names':False,'private_key_in_zip':False,'counts_reconcile':True,'zip_sha256':hashlib.sha256(zip_path.read_bytes()).hexdigest(),'upload_zip':str(zip_path),'upload_folder':str(OUT),'private_key':str(PRIVATE/'解盲对照_禁止上传.tsv')}
(PRIVATE/'PACKET_QA.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(qa,ensure_ascii=False,indent=2))