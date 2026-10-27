import pathlib,json,csv,gzip,shutil,zipfile,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8');O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');D=O/'delivery';D.mkdir(exist_ok=True)
def read(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
a=read(O/'results/route1/program_summary.tsv');b=read(O/'results/route2/program_summary.tsv');plan=json.loads((O/'ANALYSIS_PLAN.json').read_text());pool=json.loads((O/'results/route1/AUDIT_COMPLETE.json').read_text())['sensitivity_passes']
report='''# 结肠癌两方向首轮预分析报告

启动：2026-09-23｜完成：2026-09-24｜本地与 sgs10 同步归档

## 两条路线均可继续筛选，方向①的结肠癌纳入依据更明确

本轮完成实际数据下载、患者及部位核查、患者级表达比较、空间病理区域比较和细胞来源复核。两条路线均有可执行的数据基础；目前未确定可直接进入客户湿实验的新基因。

| 项目 | ① 肝转移相关肿瘤上皮状态 | ② 原发肿瘤侵袭前沿与基质界面 |
|---|---|---|
| 已实际分析 | 5位经临床表确认的结肠癌患者，原发—肝转移配对单细胞数据 | 5位患者、8张原发空间切片，含作者侵袭前沿注释 |
| 最明确的线索 | 严格QC后侵袭/ECM面板4/5升高，分化面板4/5下降 | 分化面板5/5下降；深部肿瘤对照下4/4仍下降 |
| 当前主要未决问题 | 作者肿瘤标签尚未完成独立恶性身份确认；肝相关RNA变化的来源 | 逐患者结肠/直肠信息不全；肿瘤区域spot中混有基质成分 |
| 创新性判断 | 需要落到研究不足的具体调控机制；“肝转移适应”本身不是新发现 | 需要落到明确的细胞互作机制；常规CAF标记和空间共存不足以构成创新 |
| 后续实验难度初评 | 中高：癌细胞功能及救援可行，完整肝转移模型要求较高 | 中等至中高：共培养、空间蛋白验证和救援；须明确作用细胞 |
| 推进建议 | 作为目前更明确的纯结肠癌筛选主线 | 保留为空间机制探索线，补齐纯结肠纳入依据后再定题 |

当前可以进入进一步生信筛选，不应把本报告写成已建立两个新机制或已经满足目标期刊要求。药物残存/再生长、胃癌和FABP4方向均未恢复。

## 数据清单已从论文规模落实到实际可用患者

**方向①：GSE225857。** 原文有6位单细胞患者，实际公开的非免疫细胞矩阵覆盖其中5位。临床附表Table S2（补充PDF第13页）确认P0107、P0115、P0920、P1231为左半结肠，P0813为右半结肠；P1125未进入公开非免疫矩阵。五位均有原发CCT和肝转移LCT数据，全部接受术前化疗，其中P0813为FOLFOXIRI，其余为FOLFOX。实际读取41,892个非免疫细胞、17,515个基因，原始计数逐细胞总和与作者nCount_RNA完全一致。[1]

**方向②：GSE294385＋作者GitHub病理注释。** 从公开资源取得8张原发切片：M-ST-13、15、29、30、31、32、33、34，对应患者1、6、11、11、9、9、10、10。共匹配60,801个病理注释spot，矩阵包含18,085个基因。M-ST-33没有侵袭前沿标注，不能形成该对比；主要比较实际为7张切片、5位患者。深部肿瘤对照要求至少20个spot，实际为6张切片、4位患者。[2–3]

该队列病理字段使用“Primary Colon”，临床样本表使用“Primary CRC”，未给出这5位逐人的结肠/直肠明细。因此本轮方向②按CRC探索处理，未将其算作已经确认的纯结肠癌验证。另已核实GSE225857空间样本C1–C4均为结肠原发，但本轮未下载分析这四张空间矩阵，它们是后续来源明确的补充入口。[1]

**访问状态修正。** 2026年空间多组学论文对应Zenodo记录19043057及最新版本21431802在本轮API核查中均为restricted，文件列表为空。作者已经写明存储位置，不等于当前可匿名下载。本轮只使用实际开放取得的GEO数据及GitHub注释，没有把受限的HD、WGS或蛋白成像算作已取得数据，也没有提交访问申请。[2–4]

## 方向①显示患者级状态变化，但肝相关RNA不能直接当作新靶点

下表为肝转移−原发的程序分数差，中位数来自5位独立患者。分数为面板成员平均log2(CPM＋1)。严格QC敏感性排除Tu10_COL3A1群，并要求线粒体比例≤25%、检出基因数≥200。

| 肿瘤细胞口径 | 面板 | 中位数差 | 升高/下降患者 |
|---|---|---:|---|
'''
labels={'invasion_ECM':'侵袭/ECM','intestinal_differentiation':'肠上皮分化','proliferation':'增殖','epithelial_proxy':'上皮标记'}
for ct,ctname in [('Tumor_author','作者肿瘤标签'),('Tumor_strict_QC','严格QC')]:
 for prog in labels:
  r=next(r for r in a if r['celltype']==ct and r['program']==prog);report+=f"| {ctname} | {labels[prog]} | {float(r['median_delta']):+.3f} | {r['positive']}/{r['negative']} |\n"
report+='''
全基因组按表达覆盖、效应大小和患者方向一致性作描述筛选；进一步要求宽群与严格QC两种归一化方向一致，并检查可配对的Tu01_AREG和Tu04_RGMB亚群，每个亚群覆盖4位患者。最终保留126个描述性筛选条目，完整列表附后。它们是待查新的候选池，不是126个已通过统计显著性或独立验证的靶点；这些敏感性仍来自原队列，不能追加为独立样本。

未经来源核查的前列基因包含ALB、FGA、FGG等肝脏相关转录本。本轮另建11基因肝相关RNA代理，严格QC肿瘤群中5/5患者升高，中位数差+2.971；成纤维细胞仅2/5升高，中位数差−0.124；内皮细胞2/4升高，中位数差+0.011。肝相关变化没有在所有细胞群中普遍出现，不能简单归结为统一的样本污染；但目前也不足以认定为肿瘤细胞的肝样功能重编程。

![方向①患者级变化](figures/route1_patient_effects.png)

图中每点为一位患者的配对差，红线为中位数。肿瘤及成纤维比较为5对，内皮比较为4对；无细胞级显著性检验。

## 方向②显示侵袭前沿分化降低，部分升高信号主要指向基质来源

采用作者病理标签进行切片内比较，同一患者的多张切片差值等权合并后，再计算患者间中位数。

| 对比 | 面板 | 患者数 | 中位数差 | 升高/下降患者 |
|---|---|---:|---:|---|
'''
for ct,cn in [('tumor_front_vs_other_tumor','前沿肿瘤−其他肿瘤'),('tumor_front_vs_deep_tumor','前沿肿瘤−深部肿瘤')]:
 for prog in ['intestinal_differentiation','invasion_ECM','stromal_proxy','epithelial_available3']:
  r=next(r for r in b if r['contrast']==ct and r['program']==prog);lab={'intestinal_differentiation':'肠上皮分化','invasion_ECM':'侵袭/ECM','stromal_proxy':'基质标记','epithelial_available3':'可测上皮三基因代理'}[prog];report+=f"| {cn} | {lab} | {r['n']} | {float(r['median_delta']):+.3f} | {r['positive']}/{r['negative']} |\n"
report+='''
分化面板的下降在改变比较区域后仍保持方向一致，支持继续检查侵袭前沿的上皮分化状态。侵袭/ECM面板与基质标记同时升高，提示解释时应考虑肿瘤—基质混合。前沿基质相对其他基质的侵袭/ECM面板中位数差为−0.198，2/5升高、3/5下降；因此当前数据不支持将结果叙述为“前沿所有基质细胞的侵袭程序整体增强”。

前沿肿瘤区域的描述筛选得到104个基因，其中103个可在GSE225857非免疫单细胞参考中定位。SULF1、ITGBL1、THBS2、COL10A1、FBLN2等升高基因，在5位参考患者中均显示成纤维细胞表达高于作者肿瘤群。它们说明了来源核查的重要性，不能因其出现在病理肿瘤区域就归为癌细胞内在表达，也不能仅凭本轮结果列为创新靶点。

![方向②患者级变化](figures/route2_patient_effects.png)

![作者病理区域映射](figures/route2_pathology_map.png)

空间图为患者11的M-ST-29，显示作者标注的前沿肿瘤与前沿基质位置。图示与距离计算只描述空间邻接，不证明细胞通信或因果作用。

## 两条路线的创新问题需要进一步收窄

方向①更合适的问题是：在恶性身份可靠、且匹配同类细胞状态后，哪些调控因子与结肠癌肝转移的分化/侵袭变化稳定相关？优先核查126条目的来源、已有研究和独立队列方向，再考虑RNA＋ATAC等调控证据。现有多组学静态资料不能代替转移过程中的动态调控证明。

方向②更合适的问题是：侵袭前沿的低分化上皮与哪一种明确的基质状态相邻，是否存在独立可复现的分子联系？应先补纯结肠纳入信息、分解spot中的细胞成分，再从上皮和基质各自的表达状态筛选机制。原GSE225857论文已研究F3/NRG1–ERBB3等互作，2026年论文已研究微转移干性、静息及空间微环境，不能把这些已发表主线当作本项目的新发现。[1–2]

本轮没有完成126/104个条目的逐基因全面查新，也没有宣布任何单基因为最终推荐。后续湿实验应服务于已收窄并得到生信支持的机制，而非用来解决尚未确认的细胞身份或样本来源问题。

## Methods

下载GEO处理后矩阵、公开作者细胞注释、GitHub空间病理注释及原始论文临床附表。逐条记录URL、文件大小、下载退出码、传输退出码和SHA256。9个矩阵输入共1,684,437,570字节，本地与服务器哈希全部一致；受限Zenodo资源仅记录访问状态。

方向①按作者patients与organs字段进行原发CCT、肝转移LCT配对，临床P编号与矩阵s编号按相同日期编号映射。细胞条码中的连字符/点差异经过唯一性和全集匹配核查；原始计数为非负整数，逐细胞总计数与元数据完全一致。肿瘤、成纤维和内皮分别采用Tu、F、E前缀标签；每患者每部位每群至少20个细胞才进行比较。计算聚合log2(CPM＋1)及平均单细胞CP10K，并检查排除Tu10_COL3A1、技术QC和作者同亚群的敏感性。

方向①初筛要求至少3位患者达到表达覆盖、主效应绝对中位数≥1、至少80%患者方向一致，排除线粒体和核糖体符号。来源核查后增加描述性敏感性筛选：主肿瘤群与严格QC均覆盖5对、绝对中位数差≥1、各至少4/5同向，并与平均单细胞CP10K方向一致；Tu01_AREG与Tu04_RGMB各至少3/4同向；排除记录的肝相关代理成员及VTN、UGT2B15。这是筛选后的稳健性步骤，不是独立验证。完整筛选前后结果均保留。

方向②通过sample~barcode将病理注释与每张切片矩阵、坐标一一匹配。主要对比为“Tumor cells at tumor invasive front/border”对黏膜/黏膜下、肌层、浆膜下肿瘤区域；敏感性只保留肌层及浆膜下对照；另比较前沿基质与其他肿瘤基质。两侧均至少20个spot。先将区域计数聚合至log2(CPM＋1)，计算同切片差值，再等权合并到患者。单spot平均CP10K作为表达加权敏感性。

程序按以下固定成员平均log2(CPM＋1)计分：

'''
for name,gs in plan['programs'].items():report+='- '+name+'：'+', '.join(gs)+'。\n'
report+='''
空间FFPE探针集缺少LDHA和KRT18：缺氧面板使用可测的5/6基因；原4基因上皮代理仅覆盖3/4，未达预设80%覆盖，保留为不可计算。另标记加入EPCAM/KRT8/KRT19三基因诊断代理，不以此替代原面板或恢复显著性。肝来源代理为ALB、FGA、FGB、FGG、APOA1、APOA2、TTR、HP、AGT、GC、CRP，来源异常出现后作为诊断性核查添加。

空间候选细胞来源采用另一队列每患者各细胞群的平均聚合表达，原发与肝转移在患者内汇总，用于定位而非转移方向的独立复现。所有结果以患者为单位；本轮有效独立患者少于10，故仅报告效应、方向和一致性，不报告推断性P值或FDR，不将细胞、spot、切片或同患者亚群当作新增生物学重复。

## Discussion

方向①的优势是纯结肠临床身份、配对关系和可追溯作者标签明确。其边界是仅5位术前治疗患者，无法直接推广至未经治疗的自然转移；作者肿瘤标签主要由标记支持，本轮没有完成逐细胞CNV/突变证据确认。肝相关RNA集中于Tu群的现象既可能包含肿瘤状态，也可能受注释及环境RNA影响，当前结果不能作单一归因。

方向②的优势是采用独立于待筛基因的病理区域定义，且同一患者可进行区域比较。其边界是逐患者结肠/直肠部位未齐全，空间spot并非纯细胞，基质标记升高不能直接解释为癌细胞内的调控。患者级分化下降在深部对照中一致，是探索性事实；小样本不允许把它升级为已验证机制，也不应仅因无法形成强FDR证据而抹去其方向。

两方向目前只有转录与空间定位层面的首轮实证，尚未形成新的多组学因果链。下一阶段的价值在于筛掉来源错误、重复已有研究或不能独立复现的候选。是否能支撑目标期刊，取决于后续具体机制与验证质量。

## 来源与复核位置

1. [GSE225857原始研究：Single-cell and spatial transcriptome analysis…](https://pmc.ncbi.nlm.nih.gov/articles/PMC10275599/)；临床Table S2见补充PDF第13页，治疗Table S3见第14页。 [GEO数据](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225857)。
2. [2026年空间多组学原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC13390483/)；样本映射见Table S1 Tab 1。 [GSE294385](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294385)。
3. [作者病理注释和分析代码](https://github.com/yliuup/CRC_micromets_ST)。
4. [Zenodo原记录](https://zenodo.org/records/19043057)、[最新版本](https://zenodo.org/records/21431802)：本轮API核查为restricted，原始JSON已归档。

服务器：`/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/`。
本地：`C:/Users/31937/CRC_routes_20260923/`。

证据包包含固定计划、数据状态、临床映射、完整患者级效应与筛选表、所有阴性/敏感性结果、代码、图件和校验记录；原始大矩阵保留在上述raw目录，不重复装入证据包。图件为描述性科研图，均已核查标签和布局。
'''
(D/'结肠癌两方向首轮预分析报告.md').write_text(report,encoding='utf-8');shutil.copytree(O/'results/figures',D/'figures',dirs_exist_ok=True)
clinical=[{'matrix_id':k,'clinical_id':'P'+k[1:],**v,'source':'GSE225857 supplementary PDF Table S2 page13 and Table S3 page14','paired_nonimmune_matrix':True} for k,v in plan['clinical_map'].items()];(O/'metadata/clinical_eligibility.json').write_text(json.dumps(clinical,ensure_ascii=False,indent=2),encoding='utf-8')
(D/'README.txt').write_text('Read the Markdown report first. Tables contain descriptive patient-level results, not inferential P/FDR. Original matrices are retained in the project raw directory. No novel target or wet-lab mechanism has been confirmed.\n',encoding='utf-8')
print('REPORT_READY',len(report),str(D))
