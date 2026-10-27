from pathlib import Path
import json, re, hashlib, html
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT

ROOT=Path(__file__).parent
OUT=ROOT/'output'/'pdf'; OUT.mkdir(parents=True,exist_ok=True)
pdfmetrics.registerFont(TTFont('Song',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0))
pdfmetrics.registerFont(TTFont('Hei',r'C:\Windows\Fonts\simhei.ttf'))
pdfmetrics.registerFontFamily('Song',normal='Song',bold='Hei',italic='Song',boldItalic='Hei')
W,H=595.276,841.89
NAVY=colors.HexColor('#173E50');TEAL=colors.HexColor('#157A83');GRAY=colors.HexColor('#52636D')
styles={
 'body':ParagraphStyle('body',fontName='Song',fontSize=11,leading=18.5,textColor=NAVY,wordWrap='CJK',spaceAfter=8),
 'small':ParagraphStyle('small',fontName='Song',fontSize=8.6,leading=13.5,textColor=GRAY,wordWrap='CJK'),
 'h2':ParagraphStyle('h2',fontName='Hei',fontSize=13,leading=20,textColor=TEAL,wordWrap='CJK'),
 'table':ParagraphStyle('table',fontName='Song',fontSize=9.6,leading=15.3,textColor=NAVY,wordWrap='CJK'),
 'th':ParagraphStyle('th',fontName='Hei',fontSize=9,leading=14,textColor=colors.white,wordWrap='CJK'),
 'title':ParagraphStyle('title',fontName='Hei',fontSize=23,leading=32,textColor=NAVY,wordWrap='CJK')}

def esc(t):
 t=html.escape(t).replace('\n','<br/>')
 return re.sub(r'\[R(\d{2})\]',lambda m:f'<a href="{REFS[int(m[1])-1][3]}" color="#157A83">[R{m[1]}]</a>',t)

REFS=[
('R01','2008 | Mol Biol Cell','Paracingulin Regulates the Activity of Rac1 and RhoA GTPases by Recruiting Tiam1 and GEF-H1 to Epithelial Junctions','https://pmc.ncbi.nlm.nih.gov/articles/PMC2555940/'),
('R02','2015 | J Cell Biol','ZO-1 controls endothelial adherens junctions, cell-cell tension, angiogenesis, and barrier formation','https://pubmed.ncbi.nlm.nih.gov/25753039/'),
('R03','2017 | Cardiovascular Research','Cgnl1, an endothelial junction complex protein, regulates GTPase mediated angiogenesis','https://pmc.ncbi.nlm.nih.gov/articles/PMC5852532/'),
('R04','2023 | J Cell Biol','Cingulin and paracingulin tether myosins-2 to junctions to mechanoregulate the plasma membrane','https://pmc.ncbi.nlm.nih.gov/articles/PMC10202830/'),
('R05','2023 | J Cell Sci','Paracingulin recruits CAMSAP3 to tight junctions and regulates microtubule and polarized epithelial cell organization','https://pmc.ncbi.nlm.nih.gov/articles/PMC10184829/'),
('R06','2025 | PMID 40636102','Paracingulin controls the junctional accumulation of nonmuscle myosin-2 in endothelial cells in vivo','https://pmc.ncbi.nlm.nih.gov/articles/PMC12238879/'),
('R07','2014 | Molecular & Cellular Proteomics','MicroRNA-205 Targets Tight Junction-related Proteins during Urothelial Cellular Differentiation','https://pmc.ncbi.nlm.nih.gov/articles/PMC4159652/'),
('R08','2021 | Int J Mol Sci','Identification of Tumor Suppressive Genes Regulated by miR-31-5p and miR-31-3p in Head and Neck Squamous Cell Carcinoma','https://www.mdpi.com/1422-0067/22/12/6199'),
('R09','2013 | PLoS ONE','The Junctional Proteins Cingulin and Paracingulin Modulate the Expression of Tight Junction Protein Genes through GATA-4','https://pmc.ncbi.nlm.nih.gov/articles/PMC3567034/'),
('R10','2019 | Cancer Research','Hepatic Endothelial Notch Activation Protects against Liver Metastasis by Regulating Endothelial-Tumor Cell Adhesion Independent of Angiocrine Signaling','https://pubmed.ncbi.nlm.nih.gov/30530502/'),
('R11','2013 | Cancer Cell','Endothelial Cells Promote the Colorectal Cancer Stem Cell Phenotype Through a Soluble Form of Jagged-1','https://pmc.ncbi.nlm.nih.gov/articles/PMC3574187/'),
('R12','2016 | Nature Medicine','Vessel co-option mediates resistance to anti-angiogenic therapy in liver metastases','https://pmc.ncbi.nlm.nih.gov/articles/PMC5104270/'),
('R13','2022 | British Journal of Cancer','Histopathological growth patterns of liver metastasis: updated consensus guidelines for pattern scoring, perspectives and recent mechanistic insights','https://www.nature.com/articles/s41416-022-01859-7'),
('R14','2021 | Communications Biology','Runt related transcription factor-1 plays a central role in vessel co-option of colorectal cancer liver metastases','https://www.nature.com/articles/s42003-021-02481-8'),
('R15','2023 | Gut；在线发表于2022','Novel TCF21high pericyte subpopulation promotes colorectal cancer metastasis by remodelling perivascular matrix','https://gut.bmj.com/content/72/4/710'),
('R16','2025 | Nature Communications','Pericytes promote metastasis by regulating tumor local vascular tone and hemodynamics','https://www.nature.com/articles/s41467-025-62475-6'),
('R17','2021 | Journal of Hepatology','Endothelial GATA4 controls liver fibrosis and regeneration by preventing a pathogenic switch in angiocrine signaling','https://doi.org/10.1016/j.jhep.2020.08.033'),
('R18','2024 | Nature Medicine','An organotypic atlas of human vascular cells','https://pmc.ncbi.nlm.nih.gov/articles/PMC11645277/'),
('R19','2023 | Science Advances','Single-cell and spatial transcriptome analysis reveals the cellular heterogeneity of liver metastatic colorectal cancer','https://pmc.ncbi.nlm.nih.gov/articles/PMC10275599/'),
('R20','2024 | Clinical Cancer Research','Spatial Transcriptome Mapping of the Desmoplastic Growth Pattern of Colorectal Liver Metastases by In Situ Sequencing Reveals a Biologically Relevant Zonation of the Desmoplastic Rim','https://pmc.ncbi.nlm.nih.gov/articles/PMC11443209/'),
('R21','2025 | Cancer Cell','Phenotypic heterogeneity and plasticity in colorectal cancer metastasis','https://pmc.ncbi.nlm.nih.gov/articles/PMC12278653/'),
('R22','2026 | Cancer Cell；DOI 10.1016/j.ccell.2026.06.009','Spatial multi-omics landscape of colorectal cancer macro- and micrometastases','https://www.sciencedirect.com/science/article/pii/S1535610826002965'),
('R23','2026 | Molecular Cancer；2026-07-20','Evidence of endothelial anergy and T-cell exclusion in microsatellite stable colorectal liver metastases revealed through human hemi-liver perfusion and spatial transcriptomics','https://link.springer.com/article/10.1186/s12943-026-02732-4'),
('R24','2024 | Nature','In vivo interaction screening reveals liver-derived constraints to metastasis','https://www.nature.com/articles/s41586-024-07715-3'),
('R25','2022 | Scientific Reports','Endothelial Notch signaling directly regulates the small GTPase RND1 to facilitate Notch suppression of endothelial migration','https://pubmed.ncbi.nlm.nih.gov/35102202/'),
('R26','2024 | J Mol Cell Cardiol','NCoR1 limits angiogenic capacity by altering Notch signaling','https://pubmed.ncbi.nlm.nih.gov/38359551/'),
('R27','2024 | Genome Biology','Niche-DE: niche-differential gene expression analysis in spatial transcriptomics data identifies context-dependent cell-cell interactions','https://pubmed.ncbi.nlm.nih.gov/38217002/'),
('R28','2023 | Nature Methods','SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks','https://www.nature.com/articles/s41592-023-01938-4'),
('R29','2024 | Cell Reports','Integrative multi-omic profiling of adult mouse brain endothelial cells and potential implications in Alzheimer’s disease','https://pmc.ncbi.nlm.nih.gov/articles/PMC10843806/'),
('R30','2024 | PMID 38297362','Role of a novel circRNA-CGNL1 in regulating pancreatic cancer progression via NUDT4-HDAC4-RUNX2-GAMT-mediated apoptosis','https://pmc.ncbi.nlm.nih.gov/articles/PMC10829403/')]

def P(t):return ('p',t)
def H2(t):return ('h',t)
def T(headers,rows,widths=None):return ('table',headers,rows,widths)
def C(t):return ('callout',t)
def S(t):return ('small',t)

AUDIT=[]; PLAN=[]
def audit(title,*blocks):AUDIT.append((title,list(blocks)))
def plan(title,*blocks):PLAN.append((title,list(blocks)))

audit('文献调查与 v1.1 方案改进审阅',
 C('核心判断：保留 CGNL1 选题及既有统计框架，将研究重心调整为“空间血管情境 + 上游调控证据”，并把连接转录评分恢复到它能回答的层面。'),
 P('审阅对象：用户指定《CGNL1纯生信研究计划_客户版_v1.1.pdf》，2026年9月29日，共23页。原文件保持不变。本调查截至2026年10月5日，纳入30项核心原始研究、方法或共识来源，并定向核查公开数据入口。'),
 H2('这次调查形成的五项决定'),
 P('1. 已有直接机制支持 CGNL1 参与连接骨架组织；四基因转录相关性不是这些机制的充分替代测量。\n2. 周细胞接触及 Notch 是有文献依据的上游起点；直接 TF-顺式元件-CGNL1 关系仍须逐边核查。\n3. 内皮 Notch 与癌细胞 Notch 的方向不能混写；原发灶入血、肝内黏附和已建立转移灶分开解释。\n4. 血管共选、血管生成、肝窦毛细血管化和血管成熟分别建模。\n5. CRC肝转移空间研究已经存在，增量应落在 CGNL1 位点、细胞类型、邻域和独立复现的交汇处。'),
 H2('交付结构'),
 P('本文提供证据分层、竞争研究、逐项修改对照和当前状态纠正；配套《CGNL1空间血管情境与上游调控研究方案_v2.0》提供新的可执行设计。此次交付为研究设计与审阅，没有新增患者统计结果。'),
 S('源文件SHA256：3f55cc579ae5647d46e234271cf802d64201f5e9dafe2db8ed6a2cff0ba9650e。文献编号[R01]-[R30]可点击，并在文末给出完整标题与链接。'))

audit('调查范围、证据等级与查阅边界',
 H2('定向范围综述'),
 P('检索使用 CGNL1、paracingulin、JACOP、cingulin-like 1，并组合 endothelial、Notch、pericyte、promoter、miRNA、spatial transcriptomics、colorectal liver metastasis、vessel co-option、histopathological growth pattern、perturbation。沿关键原始研究追溯参考文献，并核对2024-2026年的邻近工作和官方数据记录。'),
 P('来源优先级：原始论文全文/出版商页面、PubMed摘要、GEO/ArrayExpress等官方记录。综述和聚合网页仅作发现线索；本交付不把数据库预测、供应商说明或自动生成摘要作为机制依据。不是PRISMA系统综述，也不报告未经完整记录的总检索命中数。'),
 T(['层级','本项目中的含义'],[
 ['D1 直接机制','直接涉及线性CGNL1蛋白/表达的干预、结合、定位或救援证据；同时标出实验系统。'],
 ['D2 邻近机制','血管、周细胞、Notch或肝转移机制成立，但论文未直接证明CGNL1参与。'],
 ['D3 人体空间/组学','观察性定位、程序关联或模型推断；区分患者、切片、细胞及模态。'],
 ['D4 方法/资格','统计或技术方法、数据公开记录；不算额外生物学验证。']],[92,415]),
 H2('访问与结论边界'),
 P('部分PMC页面触发访问限制。对应条目使用已取得的原始文献摘要、检索返回的原文段落或出版商页面，未宣称逐篇补充材料已全部审完。D1重点核查操纵对象、实验读出和方向；D2/D3用于确定竞争假说和设计要求。'),
 P('“本轮未发现专门证据”只描述此次定向调查；不能据此断言全球首次。尚未核验的细胞面板、供者独立性和公开矩阵覆盖，均留在新方案资格任务中。'))

audit('直接机制支持：测量层级比面板大小重要',
 T(['证据','模型与实测结论','方案含义'],[
 ['[R01] 2008','上皮细胞；CGNL1招募Tiam1/GEF-H1，影响Rac1/RhoA活性与连接形成。','活性和亚细胞定位不等于相应基因mRNA。'],
 ['[R02] 2015','内皮；ZO-1参与JACOP及p114RhoGEF的连接招募和局部张力组织。','定位调控与表达上游分开。'],
 ['[R03] 2017','HUVEC共培养与小鼠视网膜；CGNL1影响血管管腔稳定，周细胞接触诱导其表达。','优先构建壁细胞-内皮邻域假说。'],
 ['[R04] 2023','直接结合及敲除/救援；CGNL1结合NM2A/NM2B并影响连接组织。','MYH9/MYH10表达只作背景读出。'],
 ['[R05] 2023','上皮细胞与小鼠肠上皮；CGNL1结合并招募CAMSAP3，影响微管与极性。','独立的微管机制模块，不能直接外推肝内皮。'],
 ['[R06] 2025','小鼠主动脉和bEnd.3；敲除改变连接处NM2定位，基础条件下连接外观未明显受损。','不得预设“高CGNL1=强屏障=促转移”。']],[70,245,192]),
 H2('对既有四基因评分的判断'),
 P('CLDN5、ESAM、TJP1、OCLN评分测量连接相关转录程序。它不直接测量蛋白连接定位、肌球蛋白招募、GTP结合状态、血管通透性或血流。建议原样保留作为历史预设端点，并另设有明确文献来源的新假设。'),
 P('后续RNA中将连接转录、肌动球蛋白、微管、内皮炎症黏附分开。不同模块不强行同向合成一个“功能分数”；是否存在蛋白层读出取决于公开资料实际覆盖。'))

audit('上游证据：三条已知入口与一个未闭合环节',
 T(['入口','已知证据','仍需回答'],[
 ['细胞接触/Notch [R03]','周细胞直接共培养提高内皮CGNL1；干扰DLL4或NOTCH1/4阻断诱导，非接触条件未重现同等诱导。','CRLM内是否存在相应邻域；具体配体来源、接收内皮状态及转录执行者。'],
 ['miRNA [R07][R08]','miR-205在尿路上皮、miR-31-3p在头颈鳞癌具有CGNL1 3′UTR靶向支持。','肝内皮中的表达、成熟链覆盖及作用方向，不能借用其他细胞系直接确认。'],
 ['连接招募 [R02][R05]','ZO-1/PLEKHA7相关招募框架及CGNL1-CAMSAP3作用已有研究。','表达上升是否伴随功能定位变化；RNA不单独回答。'],
 ['转录/顺式元件','本轮未核实CRLM内皮中某一TF直接结合并调控CGNL1的完整证据链。','候选TF活性、位点可及性、motif、结合及干预逐层衔接。']],[88,219,200]),
 H2('优先级与防止方向误读'),
 P('第一优先为 Notch应答及RBPJ相关候选；NCOR1作为调控复合体背景，不将其当作独立motif的序列特异TF。[R25][R26] GATA4/MYC可用于比较肝窦身份和毛细血管化解释。[R17] ERG/ETS等身份相关候选须有独立位点证据后再提升，不因表达相关直接命名为上游。'),
 P('[R09]研究的是CGN/CGNL1对GATA4及连接蛋白表达的影响，不能倒置成“GATA4已直接调控CGNL1”；也不能把CGN的证据自动转给CGNL1。circCGNL1是另一RNA分子，其胰腺癌机制不作为线性CGNL1内皮机制。[R30]'))

audit('邻近机制揭示：方向必须绑定细胞与阶段',
 T(['文献情境','发现','设计改变'],[
 ['肝内皮Notch [R10]','内皮特异Notch激活降低ICAM1相关黏附，并减少实验性肝转移。','保留保护性/适应性解释，不能预设促转移。'],
 ['内皮至癌细胞 [R11]','内皮来源可溶Jagged-1在研究系统中促进癌细胞Notch与干性。','与“Notch作用于内皮”分别建模；JAG1 RNA不证明剪切分泌。'],
 ['血管共选 [R12][R14]','既有肝窦可供替代型生长使用；癌细胞RUNX1相关运动参与共选。','纳入病理生长方式；RUNX1不能借此成为内皮CGNL1上游。'],
 ['周细胞状态 [R15][R16]','原发CRC周细胞的基质和血流调节与转移有关。','周细胞数量、基质状态、收缩状态分开，原发证据不等于肝转移灶机制。'],
 ['肝窦身份 [R17]','内皮GATA4缺失、MYC/PDGFB改变与肝窦毛细血管化和纤维化相关。','将器官身份变化作为CGNL1差异的竞争解释。']],[99,210,198]),
 P('[R15]原文还报告肝转移瘤周细胞中TCF21未检出，因此尤其不能把原发灶TCF21high周细胞标签直接移植给肝转移灶。2023与2025周细胞论文复用GSE199726的部分资料，不能算作两套独立单细胞验证。'),
 P('肝细胞来源PLXNB2的体内相互作用筛选提供“宿主细胞影响肝定植”的设计参照，但操纵对象是肝细胞PLXNB2，不是内皮CGNL1。[R24] 跨器官内皮研究也提示CGNL1定位具有组织背景，需实际核对肝内来源，而非套用脑内皮特异性。[R29]'),
 C('新的三种竞争解释：H-A局部接触诱导；H-B器官/血管亚型组成改变；H-C应激或适应性重塑。三者可以同时存在，效应方向由患者数据决定。'))

audit('空间邻近文献与真实增量',
 T(['研究','已经覆盖','本项目应补充'],[
 ['[R18] 2024 血管图谱','器官特异血管细胞状态；使用空间技术，CGNL1出现在相关空间标记面板中。','不能写“CGNL1首次被空间测量”；聚焦CRLM特定问题。'],
 ['[R19] 2023 CRLM图谱','单细胞与空间联合描述CRC肝转移异质性。','泛细胞图谱不是新增价值。'],
 ['[R20] 2024 ISS','纤维增生边缘的分区及血管生成相关邻域。','界面几何、纤维环内外及血管类别应进入设计。'],
 ['[R21] 2025 multiome','癌细胞状态、AP-1/NF-κB、免疫抑制邻域；配套RNA/ATAC与空间资料。','借助数据回答内皮CGNL1调控，但不能把癌细胞TF结果移植。'],
 ['[R22] 2026 空间多组学','宏/微转移灶、克隆及免疫生态；已用GSE294385属于该研究。','新HD模态先查重复；同患者多模态增强测量，不增加独立n。'],
 ['[R23] 2026 内皮失能','肝灌流、空间及多重成像，描述区域内皮差异与T细胞排斥。','只再现“肿瘤内皮与免疫排斥”竞争力不足；需CGNL1特异调控和邻域增量。']],[106,201,200]),
 P('目前可守住的创新表述：以已知CGNL1连接骨架机制为依据，在CRC肝转移内检验壁细胞邻域、内皮转录状态及候选顺式调控的对应关系，并在患者层面分离细胞组成与组织背景。该表述是研究目标，尚不是已获得的新机制。'))

audit('原方案逐项对照：A01-A08',
 T(['原页/模块','保留与不足','具体改进'],[
 ['p4 A01','患者登记与终点冻结正确；模态版本需更新。','新增研究原始来源、分析暴露史、同一患者跨库重复和新旧终点版本。'],
 ['p4 A02','配对bulk及组成敏感性保留；组织升高不足以独自定义主线。','引言先给机制依据；组织差异作为人体入口，明确原发CRC/肝转移对比。'],
 ['p5 A03','去目标身份判定与全谱系来源图谱正确。','加强周细胞、血管平滑肌、肝星状细胞与CAF区分；不预设内皮特异。'],
 ['p5 A04','亚型/数量分离已写入，但病理血管背景弱。','增加肝窦身份、毛细血管化、血管生成与壁细胞状态；精细标签稳定才采用。'],
 ['p6 A05','四基因患者级检验可复核，但被赋予“最核心”地位。','冻结保留历史结果；新核心改为邻域-CGNL1及上游-CGNL1两个问题。'],
 ['p6 A06','检出匹配及技术对照合理。','作为技术解释模块；保留≥10对原规则，已有宽阈值分析另标探索版本。'],
 ['p7 A07','新增患者与时间限定搜索合理。','增加每供者目标信息量、空间可测性及上游扰动的资格，不只追求人数。'],
 ['p7 A08','冻结分数和覆盖检查正确；单分数跨模态负担过重。','拆分连接转录、骨架/微管、Notch应答、炎症黏附；分别冻结、分别解释。']],[82,201,224]),
 P('A03/A04此前未完成的发现队列与GSE231559:L6复核仍是输入任务。已审核两队列的粗亚型结果无需重复运行以寻找阳性；新的机制模块由文献决定，并明确标为方案修订后分析。'))

audit('原方案逐项对照：A09-A16',
 T(['原页/模块','保留与不足','具体改进'],[
 ['p8 A09','独立病理区域、物理尺度和计算分区标签正确。','新增HGP及治疗后形态证据等级；局部切片不足时仅标界面特征。'],
 ['p8 A10','解卷积和参考敏感性完整；残余不等于纯内皮状态。','保留测量、预测、残余三层；优先用单细胞空间验证来源，传播组成不确定性。'],
 ['p9 A11','主端点偏向连接分数与肿瘤距离；未突出周细胞。','主问题转为内皮CGNL1与壁细胞邻域；肿瘤/肝距离为背景，分尺度建模。'],
 ['p9 A12','独立患者及平台可比性正确。','加入高分辨率候选与靶向面板覆盖；同研究HD不当作外部验证。'],
 ['p10 A13','两大样本队列与组成调整保留。','验证冻结机制程序和组织组成，不称bulk证明单细胞内皮调控。'],
 ['p10 A14','真实分型、治疗和交互设计保留。','补充可靠HGP、取样部位及病灶阶段；未知字段不补造。'],
 ['p11 A15','连续OS模型、左截断/landmark审查合理。','保持条件增强，降低投入；预后不是主线验收门槛。'],
 ['p11 A16','简洁模型与训练内处理合理。','默认暂缓，待真实结局、事件量及外部验证均可用再启动。']],[82,201,224]),
 P('p21“来源与状态通过后推进空间解释”应修改为：身份与测量资格通过后启动空间及上游研究，不以旧四基因相关P值为共同闸门。每条假设单独验收，避免一次转录阴性阻断不同测量层的问题。'))

audit('原方案逐项对照：A17-A22 与可选分支',
 T(['原页/模块','调整','改进内容'],[
 ['p12 A17','提升为核心','从“连接状态相关TF”转为“CGNL1目标调控”；固定Notch候选与有限竞争程序，估计TF活性时剔除目标泄漏。'],
 ['p12 A18','资格审查前置；分析条件执行','先做配对条码与内皮ATAC小试。每个候选列可及性、motif、peak-gene、结合/扰动证据。'],
 ['p13 A19','限定轴的通讯升为核心','优先壁细胞-内皮；来源、受体、响应、空间与患者一致性逐边登记。低分辨率不等于接触。'],
 ['p13 A20','优先级提升','公开真实上游扰动先于虚拟敲除；明确DLL4刺激、NCOR1干预和直接CGNL1扰动的区别。'],
 ['p14 A21','延后','不因RNA网络可运行就模拟非TF支架蛋白功能；仅保留有稳定网络时的假说生成。'],
 ['p14 A22','保留加强','增加证据边和矛盾证据登记，区分独立患者、技术复现与机制类比。'],
 ['p15 B01','暂缓','不与A16重复投入，不以复杂模型数量增加论文完整性。'],
 ['p16 B02','保留条件入口','优先相关内皮真实扰动；癌细胞药敏不能填补血管机制。'],
 ['p17-18 B03','保留可行性储备','一般CRC风险与肝转移不等同；没有合格cis工具及共定位不强行展开。']],[87,96,324]),
 P('p1权重建议从30/25/20/10/15调整为：来源与状态20%、空间30%、上游与公开干预30%、组织临床10%、稳健性复现10%。比例是资源分配，不是统计权重。p23参考资料需由泛称升级为可追溯标题、年份、链接和具体支持边。'))

audit('当前证据必须进入修订起点',
 T(['项目记录','实际状态','修订含义'],[
 ['审阅后A05/A06','18患者主连续β=-0.1134，P=0.3660，FDR=0.4130；同7患者深度匹配差异+0.1011，P=0.2869。','现有主检验不支持稳定正相关；保留完整结果，不重命名为新验证。'],
 ['预设小样本分支','≥50细胞、研究/深度校正n=8，β=+0.1867，P=0.03231，FDR=0.1667。','保留名义显著探索信号；不能取代主分析。'],
 ['A04粗亚型','18患者双比例校正β=-0.1167，P=0.4046；同亚型匹配7人差异+0.0414，P=0.6478。','亚型粗分未恢复总体正向；不是整个A04已完成。'],
 ['multiome RNA','15患者，作者内皮QC后489核；259核通过旧严格交集；作者≥50核患者3，严格交集0。','旧交集不是身份真值；审核作者与目标排除证据，不靠降阈值造阳性。'],
 ['ATAC原始下载','304/304完成并有官方MD5校验完成记录、退出码0。','下载完成不等于内皮ATAC合格；先小试，不再把“下载未完成”作为障碍。']],[89,233,185]),
 S('内部依据：review_consensus_retest_20261002_v1、A04_subtype_check_20261002_v1、formal_v1_20261002/02_multiome/COMPLETE.json；2026-10-05已读取。原v1.1第2页是9月旧快照，不应无日期地沿用。'),
 P('GSE294385官方记录为37个组织样本/11患者的Visium系列；2026正式文章扩展为49肿瘤/19患者的多模态研究。二者范围不同，不能将19直接替换为本项目的空间n。先做发表版本-模态-患者映射。[R22]'))

audit('新增数据入口与调查后的行动优先级',
 T(['入口','已核实到的范围','下一步资格'],[
 ['GSE163568 [R25]','官方元数据含HUVEC DLL4包被刺激6h：对照3个rep、处理3个rep；另有HREC、EGTA及小鼠子实验。','核实供者/培养批次独立性；各子实验分开。原始count表已列出，尚未分析CGNL1结果。'],
 ['GSE237653 / GSE242671 [R26]','同一NCOR1研究的HUVEC RNA与ATAC；分别有对照/干预，标注n=3。RNA附件为NormCounts。','验证重复、样本配对及矩阵尺度；不可将归一化表当整数count输入NB模型。'],
 ['E-MTAB-13651/13652/13655 [R21]','RNA/ATAC与空间同研究；Mendeley数据yd7bb3shn5.1。','本地已使用部分RNA；再次核查既有/新增模态及患者重复。'],
 ['GSE294385 / Zenodo [R22]','论文提供Visium及其他模态入口；Zenodo DOI 10.5281/zenodo.19043057。','HD数据文件、肝转移覆盖、CGNL1检出和同患者重叠均待审核。'],
 ['2026内皮空间研究 [R23]','论文报道CosMx及多重蛋白验证。','尚未确认可直接获取的完整分析包；先查面板CGNL1覆盖，不预设可作目标验证。'],
 ['GSE199726 [R15][R16]','原发CRC周细胞参考，被两篇研究复用。','仅作发送细胞状态参考；不计新增CRLM内皮患者。']],[114,203,190]),
 C('优先顺序：保存旧结果 → 补齐来源和周细胞身份 → 同时审核空间与上游数据 → 冻结新端点 → 空间/上游分别执行 → 独立验证 → 组织临床补充。'),
 P('这次调查不启动大规模ATAC重算、不替换服务器原方案或统计状态。新方案独立版本化；已完成结果、旧P值及检验家族原样保留。'))

plan('CGNL1空间血管情境与上游调控研究方案',
 C('v2.0 | 文献驱动修订版 | 2026年10月5日\n研究范围：结直肠癌肝转移；不依赖新增湿实验。'),
 H2('研究目标'),
 P('利用公开单细胞、空间转录组、合格染色质与真实扰动资料，确定CGNL1相关内皮状态的细胞来源和局部血管情境，检验壁细胞邻域、Notch应答及候选顺式调控是否构成一致的关联证据，并评价组织层面的可迁移性。'),
 H2('建议工作题目'),
 P('结直肠癌肝转移中CGNL1相关内皮状态的空间血管情境与上游调控线索。\n若获得多层一致证据，可进一步具体到某条候选轴；若只有空间组成证据，则采用来源与血管邻域题目。'),
 T(['模块','投入','核心产出'],[
 ['来源、身份与测量','20%','患者级来源、周细胞/肝星状细胞区分、目标信息量'],
 ['空间血管情境','30%','壁细胞邻域、组织界面、组成分离和外部患者效应'],
 ['上游与公开扰动','30%','受体/响应、候选TF-元件-基因边、真实刺激响应'],
 ['组织与临床扩展','10%','冻结程序的组织复现及条件临床关联'],
 ['稳健性与复现','10%','版本、患者去重、统计表和完整证据链']],[160,55,292]),
 S('资源比例不是统计权重。本文件为执行设计，尚未产生新患者分析结果。旧v1.1结果继续有效保存；本版本新增假设不能改称旧分析的预设验证。'))

plan('机制依据与可检验问题',
 P('CGNL1已有连接处GTPase调节、肌球蛋白和微管组织的实验基础；已有周细胞接触/Notch依赖表达诱导证据。[R01]-[R06] 新研究围绕这些功能发生的组织条件及调控来源展开。'),
 T(['编号','预设问题','主要测量'],[
 ['H1 空间主问题','同一患者内，壁细胞更邻近的内皮单元是否具有不同CGNL1表达？','内皮归属CGNL1与壁细胞邻域暴露的患者级效应；双侧检验，正向为文献预期。'],
 ['H2 上游主问题','内皮Notch应答与CGNL1表达是否在独立患者中关联？','患者级CGNL1与不含CGNL1的独立Notch响应评分/调控活性之关联。'],
 ['H3 调控解释','候选TF是否有指向CGNL1的顺式调控和干预支持？','RNA活性、CGNL1位点可及性、motif、peak-gene连接、结合与真实干预分层记录。'],
 ['H4 情境修饰','H1/H2是否依赖血管类别、病理界面、生长方式或治疗？','预设少量交互及分层效应；条件数据合格后执行。'],
 ['L1 历史问题','既有CGNL1与四基因连接分数关系是否可复现？','原CLDN5/ESAM/TJP1/OCLN定义和旧模型保留，单独登记历史家族。']],[77,209,221]),
 C('H1、H2分别属于空间邻域和患者内皮调控的假设家族，不跨家族合成“总显著性”。H3是逐边证据任务；单个motif不作直接调控结论。'),
 P('同一患者、同一血管状态可以兼有局部接触诱导、肝窦身份变化和应激适应。Notch作用于内皮或癌细胞、肝内黏附或原发灶入血必须分别叙述。[R10][R11][R15][R16]'))

plan('版本、既有结果与新分析边界',
 H2('从已完成结果续接'),
 P('两个审阅队列的主连续分析为18位患者，β=-0.1134、P=0.3660、FDR=0.4130；A04粗亚型校正未恢复总体正向；严格同亚型匹配≥10对只剩2位患者，按描述性呈现。n=8深度校正分支的P=0.03231、FDR=0.1667保留为探索性线索。'),
 P('这些结果继续进入论文证据表；不通过增加面板、改变细胞阈值或重命名亚群来替换旧主分析。新端点专门回答空间邻域与上游调控问题。'),
 H2('三个冻结时点'),
 P('F0：仅依据文献和输入资格冻结本方案版本、目标、对比、方向、统计单位和家族。\nF1：在查看新目标关联前，冻结身份、病理区域、空间尺度、基因面板、TF集合、候选调控窗口和模型。\nF2：发现阶段结束后冻结候选轴、全部参数和验证集；验证阶段只运行一次主版本，敏感性预先列出。'),
 T(['登记项','要求'],[
 ['分析暴露史','已看过CGNL1目标结果的队列记为已探索；新命名不重置独立性。'],
 ['供者/病灶/切片','以patient_id串联组织、模态、处理、时间点与原始来源；多模态不新增患者。'],
 ['旧/新规则','原≥10对匹配、既有距离端点和检验家族保留；新设计另设analysis_id。'],
 ['变更记录','记录原因、日期、是否已看目标结果及影响范围；事后诊断不替代冻结模型。']],[111,396]),
 S('数据阶段只做资格与精度评估，不以显著性确定主面板、阈值或验证队列。'))

plan('数据分工：复用资料与空间候选',
 T(['资料','角色','进入条件'],[
 ['GSE225857；GSE231559:L6','发现/已探索补齐','优先完成同口径身份与来源复核；核对治疗和目标原始计数。'],
 ['GSE245552；GSE285991','已审阅机制探索','保留接受的主要谱系和粗亚型结果；新上游问题仍标已探索队列。'],
 ['GSE234804等未用目标资料','候选新增验证','核实人源、组织、患者与既有来源重叠；不能仅按GEO名称认定独立。'],
 ['GSE294385','空间发现','官方Visium系列与2026文章多模态版本逐患者映射；已有11人不能扩写为19人。'],
 ['GSE206552/254364/217414、Niche-DE来源','外部空间资格池','核对源病例、H&E、坐标、目标覆盖与重复；已看目标的部分属已探索支持。'],
 ['2026宏/微转移HD资料 [R22]','高分辨率候选','先检验处理数据可得性、肝转移覆盖、细胞分割与CGNL1检出；同患者仅技术支持。'],
 ['2026内皮失能空间资料 [R23]','独立候选/情境参考','实际可获取性和靶向面板先审；面板无CGNL1时只验证邻域背景。'],
 ['正常肝、血管图谱 [R18][R29]','身份和器官背景','与疾病参考并用；正常参考不能承担疾病调控真值。']],[140,110,257]),
 P('每个端点分别列患者数、病灶数和切片数。阶段性搜索限定5个工作日形成资格表；未就绪的新数据可留为后续资源，不无限延长主体任务。患者数不足的端点仍可交付探索性效应。'))

plan('数据分工：调控、扰动与组织扩展',
 T(['资料','本轮确认','资格与用途'],[
 ['CRCLM multiome [R21]','RNA已用；ATAC原始304文件下载完成。','RNA作者内皮QC后489核，≥50核供者仅3。先恢复配对条码与ATAC内皮覆盖；不默认可做患者级网络推断。'],
 ['GSE163568 [R25]','HUVEC对照3rep/DLL4包被6h 3rep；另含HREC和小鼠实验。','首选真实上游响应；确认独立培养/供者及批次。各细胞类型、时间与物种分开。'],
 ['GSE237653、GSE242671 [R26]','NCOR1干预的RNA/ATAC，同研究、每组标注3rep。','RNA附件为NormCounts；先寻合适原始计数或变换表达模型。ATAC不是额外独立队列。'],
 ['GSE199726、GSE200064/200065 [R15]','原发CRC周细胞及TCF21相关资料。','用于发送者状态参考和方法对照，不验证肝转移内皮CGNL1。'],
 ['GSE252596、GSE248500 [R16]','NKX2-3周细胞RNA/ChIP入口。','仅在发送者状态有需求时审核；不预先扩大为第二个主机制课题。'],
 ['GSE159216；GSE73255','原计划171/167患者组织队列。','核实患者去重、探针、冻结分数覆盖；OS等只在实际字段支持时使用。']],[131,177,199]),
 P('文件已公开、已下载、可读取、端点可测、有效重复充分是五个不同状态。所有入口输出逐项资格表；本次文献调查尚未读取这些新扰动的CGNL1差异结果，后续方向由实际分析决定。'))

plan('M01 来源与测量：先保证目标属于谁',
 H2('目的与流程'),
 P('构建患者级CGNL1来源图谱，覆盖内皮、恶性/非恶性上皮、肝细胞、胆管、周细胞/平滑肌、肝星状细胞/CAF及免疫群。复用既有QC并记录版本；疾病参考与正常肝参考分别提供身份支持。'),
 P('身份聚类、参考特征和组织分区排除CGNL1及随后要检验的目标模块。周细胞不能仅凭PDGFRB或ACTA2命名；结合RGS5、CSPG4、MCAM、NOTCH3等组合、平滑肌/CAF/肝星状细胞程序及可用血管定位。标记名单在F1阶段定稿，不要求每个细胞都检出所有标记。'),
 H2('新端点测量规则'),
 P('新患者级RNA主分支暂定每病灶≥20个经审核内皮细胞/核，≥10及≥50作测量敏感性；这是v2新端点规则，不覆盖旧A05。先病灶汇总再患者等权合并。保留CGNL1测得0的患者，分别报告总计数、检出率、文库量和变异；不只纳入CGNL1阳性细胞。'),
 P('原始整数计数用于检测率与计数模型；作者decontX等调整矩阵仅在明确尺度后用于适当连续表达或注释。缺失特征、测得0、QC失败、全组无变异分别记录。snRNA与scRNA按平台分层，技术零不能写成生物学缺失。'),
 T(['产出','验收'],[
 ['患者×谱系来源矩阵','包含作者标签、接受标签、支持证据、分歧、目标计数及有效n。'],
 ['周细胞/CAF/星状细胞核查','不清楚的细胞保留“壁细胞样/间质未定”标签，不强制用于接触主分析。'],
 ['测量精度表','固定深度下重复抽样仅评价测量稳定性，不增加生物学重复或产生患者效应P值。']],[134,373]),
 S('对应原A03/A04/A07；图1来源与有效信息量，完整QC置补充。'))

plan('M02 分开定义转录程序与功能读出',
 T(['测量模块','预设内容','解释'],[
 ['历史连接面板L1','CLDN5、ESAM、TJP1、OCLN；不含CGNL1。','连接相关转录分数，原定义原家族保留。'],
 ['Notch响应H2','优先使用固定版本、有人内皮刺激支持的独立响应集合；HES1、HEY1/2、NRARP等为核查种子。','最终成员/符号方向由文献和参考资源冻结；不把受体丰度等同信号活性。'],
 ['骨架与微管','MYH9/MYH10、TIAM1/ARHGEF2/ARHGEF18、CAMSAP3等分模块列出。','候选表达读出；不命名为GTPase活性或肌球蛋白定位。'],
 ['炎症黏附/血管状态','ICAM1/VCAM1/SELE等响应与肝窦、血管生成、收缩/基质程序分开。','辅助检验竞争解释，不合成为任意“促转移分数”。'],
 ['蛋白或原位成像','仅限公开面板实际包含CGNL1或相关定位蛋白且图像可定量。','直接定位读出单列；没有CGNL1蛋白不能用RNA图代替。']],[105,208,194]),
 P('每个面板记录版本、来源、权重、基因方向和平台覆盖。参考标记、CGNL1暴露、TF活性网络和结局程序之间做重叠表；主检验剔除CGNL1泄漏，重叠基因敏感性单列。小面板以固定等权log尺度均值为主，秩评分为敏感性；不能按最小P选择。'),
 P('冻结文件至少包含gene_symbol、program_id、direction、weight、evidence_source、platform_coverage、exclude_reason。历史四基因面板要求全部可测；新Notch程序在看目标结果前确定最小覆盖和主算法，未达到覆盖的平台不通过随意删基因补齐。'),
 P('新扩展程序在发现资料以外的文献/公共扰动证据定义；若用发现数据训练，则在外部冻结验证。程序存在不依赖所有成员同向，机制种子表达图与经验证的程序评分分开输出。'),
 S('对应原A05/A08/A17；图2为模块定义、跨平台覆盖及患者效应。'))

plan('M03 空间资格与血管结构背景',
 H2('区域与生长方式'),
 P('在不查看CGNL1和目标程序的前提下，标注肿瘤、肝组织、间质、坏死区及可用纤维环。HGP优先复用作者可靠病理判读；如有完整且足够的H&E界面，按2022共识记录各生长方式占比。[R13] 无专业判读的自动分区称计算分区，不把肝窦基因高直接等同replacement型。'),
 P('局部取样只能支持局部界面特征，不能据一小块切片判定整个病灶“纯desmoplastic”。治疗前后、宏/微转移灶、正常肝与纤维化背景分别记录。血管共选、血管生成和肝窦毛细血管化作为不同机制轴。[R12][R14][R17]'),
 H2('尺度预案'),
 T(['平台','主尺度/单位','可解释范围'],[
 ['标准Visium','壁细胞邻域半径150 μm；100/250 μm为敏感性。去除中心spot，优先使用连续邻域量。','区域共邻近；不解释为细胞直接接触。'],
 ['可靠细胞分割的原位平台','内皮至最近壁细胞的边界距离，20 μm内比例为预设邻近指标；10/30 μm敏感性。','近距离相邻；仅有可靠膜/边界形态时另评直接接触。'],
 ['HD bin未可靠分割','按实际bin大小汇总到预定血管ROI，报告bin与ROI两个尺度。','不将2 μm bin直接当一细胞。']],[99,234,174]),
 P('以上为v2计划值，F1前用平台几何、覆盖和分割误差作无目标资格核查；若平台不支持，端点标不可评估而非换尺度追阳性。旧固定距离和相对25%端点仍单独保留。'),
 P('空间可估计预案：每切片至少30个合格内皮索引spot/细胞或血管ROI，且邻域暴露的低/高三分位各至少10个；连续暴露变异和区域覆盖同时通过。按聚集性折算有效信息量，30个位置仅为技术门槛，不是30个生物学样本；20/50索引单元为预设资格敏感性。'),
 S('对应原A09/A11/A12；图3展示区域、几何、血管背景，全部切片列补充。'))

plan('M04 空间主分析：壁细胞邻域与内皮CGNL1',
 H2('高分辨率资料'),
 P('以独立注释的内皮细胞/血管ROI为索引，量化周围壁细胞密度、最近距离和覆盖关系。主暴露为预定邻域壁细胞量，主结局为内皮归属CGNL1表达；控制测序深度/细胞面积、区域和局部总细胞密度。身份不确定或分割边界疑似串扰的单元作预设敏感性。'),
 H2('spot资料'),
 P('使用不含目标面板的参考特征估计组成；主方法加一种不同原理方法作敏感性。先评估混合重构和参考供者留一，再用细胞类型特异的空间表达模型检验邻域效应，Niche-DE可作为候选实现。[R27] 原始表达、组成预测表达及模型归属效应分别展示。'),
 P('每切片分别估计效应，患者内按切片-病灶层级等权汇总；主群体统计单位为患者。记录含足够内皮与壁细胞邻域变异的有效患者。以空间块置换/保留区域结构的零模型检查切片内关联；不得将随机打乱spot的P值作为患者重复证据。'),
 H2('预设调整与对照'),
 P('基础模型调整技术量、肿瘤/肝距离、内皮及肝/间质组成；壁细胞暴露不再同时以同义变量重复调整。Notch响应可能位于机制路径中，不进入H1基础调整模型；加入Notch的递进模型只作分解描述，不称因果中介。'),
 P('对照包含：非邻近同片区域、丰度与检出率匹配基因、参考标记排除敏感性及区域内空间零模型。比较周细胞、血管平滑肌和星状细胞/CAF相邻效应，避免把所有基质统一叫周细胞。'),
 S('H1报告患者效应、方向一致性、区间、原始P和家族。纯spot模型只给模型归属证据；单细胞空间为独立测量加强。对应A10-A12。'))

plan('M05 上游主分析：Notch响应与候选调控',
 H2('患者级主模型H2'),
 P('在合格内皮内按病灶伪bulk计算CGNL1连续表达及独立Notch响应量，再合并至患者。主模型为CGNL1标准化表达 ~ Notch响应 + 研究；有足够可辨识信息时加入内皮粗亚型和技术量。小样本采用最简模型，各研究效应单列；研究完全混同治疗或平台时不声称统计校正已分离二者。'),
 P('细胞与核不直接拼接为统一计数模型。目标/应答全零或无变异时登记不可估计；缺失协变量不虚构。若同队列贡献H1和H2，两个模态的结果属于同患者多层支持，不称两次独立复现。'),
 H2('候选上游排序'),
 T(['层级','候选入口','升级要求'],[
 ['首选文献轴','DLL4/NOTCH1/4相关接触与RBPJ响应。','H2目标效应及发送/接收表达；随后位点证据。'],
 ['调控复合体','NCOR1及Notch相关转录调控背景。','真实扰动下CGNL1响应；NCOR1本身不作motif扫描对象。'],
 ['竞争解释','GATA4/MYC相关肝窦身份；有限ETS/应激TF候选。','各自独立网络和CGNL1位点连接，不据身份相关直接认定因果。'],
 ['转录后储备','miR-205、miR-31-3p。','细胞来源匹配的成熟miRNA测量和重复；普通3′scRNA缺少此测量时不运行ceRNA链。']],[94,173,240]),
 P('所有TF在预先固定同一家族内评价；有有效目标检验时原始P<0.05进入名义候选，完整FDR保留。按方向、效应、供者一致性和位点证据排序；数据库网络没有CGNL1覆盖时不能靠默认插补补出调控边。'),
 S('对应A17/A19；图4展示发送者-接收者-响应及患者模型。'))

plan('M06 ATAC与CGNL1顺式调控：先小试再扩展',
 H2('数据资格与计算顺序'),
 P('先整理已下载multiome的样本索引、RNA/ATAC条码映射、参考基因组和作者QC；选两个内皮覆盖较高且身份明确的文库作技术小试，不根据CGNL1效应选择。优先已有fragment/peak文件；只有原始reads时才比对。先记录单文库运行时间、内存和临时空间，再决定扩展。'),
 P('按文库评估有效片段、TSS富集、FRiP、核小体周期、双细胞及RNA身份一致性。阈值由平台规范与不含CGNL1的QC分布冻结。每供者内皮核/片段覆盖单列；作者内皮与旧严格交集作为不同身份敏感性，不把两者之一自动当真值。'),
 H2('预设位点与证据边'),
 P('固定同一参考版本下CGNL1转录本/TSS清单：主启动子TSS前2 kb至后500 bp；候选远端元件先审TSS上下游500 kb，窗口外仅在独立染色质连接证据支持时单列扩展。不同TSS分别登记，不择最显著启动子。'),
 T(['证据','实际检验/输出'],[
 ['位点可及性','各供者覆盖轨迹、峰可重复性及相关内皮亚型分布。'],
 ['motif与TF','开放峰内motif、候选TF表达/活性；共享motif先报告TF家族。'],
 ['元件-基因','同核配对RNA/ATAC、供者内连接及留一稳定性；稀疏下仅描述。'],
 ['直接结合/干预','匹配细胞背景的ChIP/CUT&Tag及扰动；SCENIC+输出是推断边。[R28]']],[109,398]),
 P('少量供者不使用细胞数替代n做患者组间显著性。可以交付可靠基因座描述，不能将失败的ATAC网络包装为机制阴性。未获得合格内皮ATAC时，正文定位为多尺度转录组及公开扰动研究。'),
 S('对应A18；图5基因座证据与候选边，资源消耗写入技术日志。'))

plan('M07 公开真实扰动：把上游关系落实到处理效应',
 H2('第一批任务'),
 P('优先GSE163568的HUVEC DLL4包被6h对照，HREC子实验作为不同内皮背景的支持；EGTA激活、抑制剂和小鼠脑实验分别建模，不合并为一个处理。其次核查GSE237653/242671的NCOR1干预RNA/ATAC，以检验更间接的调控背景。[R25][R26]'),
 P('正式分析前检查每个rep对应的供者、独立培养/动物和批次；同供者技术重复先合并。独立培养实验单位与独立人供者分别报告。处理方案、时长、对照、siRNA序列及RNA与ATAC重复配对关系写入manifest。'),
 H2('分析与验收'),
 P('CGNL1作为预设单靶点处理效应；独立Notch阳性响应基因用于确认处理确实激活/抑制预期程序。整数count采用合适的计数差异模型；若只有NormCounts且无原始count，则核查其归一化定义并采用适用连续模型或作者完整统计，不补造原始计数。'),
 P('单一预设处理对比的CGNL1无需FDR；多个细胞类型/处理/时间的CGNL1对比组成单独家族。全基因差异、程序富集、开放峰分别校正。有效独立样本少于10时定位探索性；没有有效重复只报告变化方向和效应。'),
 T(['情形','交付与解释'],[
 ['处理阳性对照响应，CGNL1同向','支持该实验条件下上游响应；再与人体H2及位点证据对照。'],
 ['处理有效，CGNL1不变或反向','保留条件依赖或阴性；不得删除研究来拼接预期轴。'],
 ['目标不测/处理未有效/无重复','分别记为测量不足、操作响应不足或仅描述，不生成目标推断P值。']],[190,317]),
 S('公开上游刺激不是直接CGNL1功能扰动；不能据响应推断其介导整个处理效应。对应A20，优先于A21。'))

plan('M08 外部验证与跨层证据整合',
 H2('验证独立性'),
 P('验证集须与发现病例无重叠，且不参与新面板、阈值和模型选择。若只能使用已探索队列，称冻结流程复测或跨队列支持。同一患者的Visium、HD、RNA与ATAC分别增加测量层，不增加独立样本量。'),
 P('外部H1优先验证相同邻域效应；若平台分辨率不同，分别验证区域共邻近或细胞近距离，不把二者数值直接合并。外部H2使用冻结面板、模型和相同效应方向，若关键基因未覆盖则该端点不可评价。'),
 H2('证据链逐边验收'),
 T(['边','最低可报告证据','更强支持'],[
 ['发送细胞-内皮','两类细胞及配体/受体可测。','独立患者空间邻域及接收者响应。'],
 ['Notch响应-CGNL1','患者级H2效应及完整统计。','真实DLL4等处理下目标响应一致。'],
 ['TF-元件-CGNL1','候选活性、开放元件、motif与连接。','匹配细胞结合和直接扰动；逐项标明是否实测。'],
 ['CGNL1-血管表型','程序或空间关系。','公开CGNL1直接扰动及实际定位/功能读出。'],
 ['血管情境-临床','患者级观察关联。','独立队列、可靠治疗/随访信息；仍不等于治疗效应。']],[120,192,195]),
 P('不将不同家族P值合成“总显著性”。同一表达矩阵派生的通讯、TF和程序分数具有依赖性，证据矩阵同时列出底层数据来源。仅在效应定义和平台可比时作预设meta分析；研究少时异质性估计作为描述。'),
 S('对应A07/A12/A22；图6区分直接实验、人体关联、模型推断和未支持的边。'))

plan('M09 组织和临床：验证群体意义',
 P('保留GSE50760/GSE14297已核实的配对组织证据，并明确18对CRC与13对严格结肠来源。肝背景调整前后效应并列；bulk差异用于来源研究起点，不直接证明癌细胞或内皮中的变化。'),
 H2('大样本组织关联'),
 P('GSE159216与GSE73255分别预处理；技术重复合并、同患者多病灶主分析等权。按注释/质量预定探针汇总，不按目标P挑探针。检验CGNL1与冻结Notch/连接等程序的联系，并递进加入不含目标的内皮、肝、间质组成量。'),
 P('若需评估壁细胞信号，在组成设计中避免同时放入高度同义分数；报告共线性、残余及平台差异。组织关联只称组织程序支持，不当作内皮特异上游验证。'),
 H2('临床条件分析'),
 P('分型、治疗、同步/异时转移、可靠HGP等按实际字段分析；LMS与CRIS不直接合并为同一类别。OS为条件主临床终点，CSS等为次要，先核实事件与表达取样相对随访起点。日期充分时采用适当landmark或左截断；没有字段不补造。'),
 P('暴露按连续1 SD建模；临床协变量由先验和事件数决定。多个预设暴露的OS效应在同一临床家族BH校正。报告HR、区间、P、FDR、事件数与比例风险检查；KM采用预设展示分组，不优化截点。'),
 C('临床生存阳性不是空间/调控研究的启动或完成条件。A16/B01预测、B02药物预测、B03遗传分析保留资格储备，本阶段默认不投入正式建模。'),
 S('对应原A02/A13-A16/B01-B03；主文以人体机制情境为主，临床证据合格后再组织。'))

plan('统计方案：估计对象、独立样本与模型',
 H2('样本单位'),
 P('人体分析以独立患者为单位，病灶、切片、ROI、spot、细胞和多模态均为嵌套测量。体外按有效独立培养/供者、动物按独立动物记录；技术重复合并。无有效重复、全零或无变异模型不作推断。少于4名有效患者的总体关联仅描述；4-9名可估计时仍定位探索性。'),
 H2('主要估计量'),
 T(['端点','效应量','主要模型'],[
 ['H1 邻域','每1 SD壁细胞邻域量对应的内皮CGNL1标准化变化；单细胞与spot分开。','切片/病灶模型先估效应，再患者等权汇总；研究分层。'],
 ['H2 上游','每1 SD Notch响应对应CGNL1变化，报告标准化β及原尺度斜率。','患者级研究调整线性模型；分队列与小样本稳健区间/研究内置换。'],
 ['处理响应','CGNL1 log2FC和区间，完整处理对比。','符合数据尺度的计数/连续模型；配对与批次按真实设计。'],
 ['TF/区域','目标关联β、可及性变化和位点支持。','患者级推断；细胞内网络或motif分数标模型统计。']],[78,236,193]),
 P('采用双侧检验；预期方向与实际方向分别记录。多因素模型用目标效应或目标对比的P，交互以交互项P评价，不比较两个亚组显著与否。技术深度、研究、区域是预设调整来源；可能的中介变量不机械塞入基础模型。'),
 P('缺失处理按字段缺失机制及样本量选定；高缺失变量退出主模型并报告，不能在全数据先做有结局监督的插补。置信区间和重采样在患者层进行。'))

plan('统计方案：检验家族、功效与决策',
 T(['家族','检验范围','口径'],[
 ['F-SPATIAL','冻结H1在主分析层的一项主要总体效应；其他细胞类别/尺度另作敏感性家族。','单个预设主要终点不需FDR；各外部队列及合并若均作推断，另定义复制家族。'],
 ['F-NOTCH','冻结H2主总体效应；受体、发送者和竞争程序单列家族。','不把H1和H2混入一个池；每队列/亚型对比提前列表。'],
 ['F-TF / F-PEAK','全部预设且可估计TF目标关联；可及性峰各自家族。','在全体受检成员校正，不只对P<0.05成员重算FDR。'],
 ['F-PERT / F-DE / F-PATH','目标处理对比、全基因差异、通路分别成家族。','模型P、富集P与生物学处理P分别标记。'],
 ['F-CLIN / F-SENS','临床终点与各范围敏感性分析。','旧A05/A06家族保持原值；新家族不能覆盖旧FDR。']],[91,264,152]),
 P('具备有效统计检验时，目标原始P<0.05为名义候选筛选硬标准；FDR未通过仍可进入后续，但明确标为名义显著/探索性。涉及多重检验同时报告FDR，只有通过对应阈值才称校正后显著。不将无有效统计属性的网络分数补成P值。'),
 H2('人数目标不是功效保证'),
 P('以≥10名有效独立患者作为报告分层界线；不承诺10例足够。F1阶段在不查看新效应的条件下，用可得患者数、测量误差和ICC范围模拟区间宽度及可检测效应；报告多个合理效应场景，不只选择乐观场景。小样本离散检验最小P限制单列。'),
 P('是否继续由数据可测性与问题相关性决定，不能用FDR未通过否定明确方向；也不能仅凭一次P<0.05宣布机制链成立。方向、区间、患者一致性、独立性与模型适用性共同解释。'))

plan('执行顺序、资源预算与阶段验收',
 T(['阶段','工作量估算','完成标准'],[
 ['S0 资格与冻结','3-5工作日','文献/数据矩阵、去重、患者覆盖、H1/H2规则、正式统计家族；此文献调查已完成，具体输入资格仍待执行。'],
 ['S1 来源和测量','4-6工作日','补齐未审队列、壁细胞身份、原始/调整计数及目标信息量。'],
 ['S2 空间邻域','6-9工作日','独立分区、组成诊断、各患者H1效应、尺度/空间零模型；有资格即执行。'],
 ['S3 上游与真实扰动','5-8工作日','H2患者关联、DLL4/NCOR1处理响应、候选边；不等待旧L1转正。'],
 ['S4 ATAC技术小试','2-4工作日','两个文库条码/QC/基因座覆盖、单文库资源估计和扩展判断。'],
 ['S5 独立复现及组织整合','4-6工作日','冻结验证、证据冲突清单、组织程序；临床仅字段合格时。'],
 ['S6 图件与稿件','4-6工作日','约6张主图、完整统计、Methods/Results/Discussion及复现包。']],[117,101,289]),
 P('基础人工工作量估计28-44工作日；S2/S3可在输入独立且资源允许时并行，实际日历工期由任务依赖与资源决定。ATAC全量比对、受控数据申请和大规模图像重分割不含在基础工期，按S4实测另估。'),
 P('当前BioDisk约剩51 GB，BioData项目盘约有5.6 TB空闲；ATAC原始文件位于BioData。临时比对/排序和缓存先规划到已核实的可写大盘，按试运行峰值预留；不删除旧结果换空间。'),
 S('估算是工作量范围，非交付日期承诺。阶段结束核对退出码、完成标记与输出完整性；不以进程消失当成功。'))

plan('图件蓝图、叙事顺序与交付',
 T(['主图','回答的问题','核心内容'],[
 ['Fig.1','CGNL1来自谁、测量是否充分？','机制依据简图、队列去重、患者来源与信息量。'],
 ['Fig.2','血管状态与旧连接结果如何对应？','内皮/壁细胞身份、独立程序、旧主分析与探索分支。'],
 ['Fig.3','CGNL1位于什么血管邻域？','独立病理分区、距离/邻域、H1患者效应与外部验证。'],
 ['Fig.4','何种上游信号与CGNL1相联系？','发送/接收及Notch响应、H2患者效应、竞争解释。'],
 ['Fig.5','调控是否获得正交支持？','合格ATAC位点、真实扰动、TF-元件-基因证据；不合格面板不强做。'],
 ['Fig.6','哪些跨层结论成立？','组织复现、条件临床、证据链支持/冲突/未知与最终模型。']],[68,177,262]),
 P('主文按实际结果采用结论性标题；先给方向、效应、患者一致性及统计证据。Methods完整描述数据和方法。样本、测量、因果与外推边界集中在Discussion；完整QC、全切片和敏感性放补充。会改变主线的阴性与验证失败保留正文。'),
 P('交付：患者/病灶/切片/模态清单；基因面板与模型冻结文件；逐患者效应；全部P/FDR/家族；空间掩膜和坐标；调控证据边；脚本/环境/日志/哈希；可编辑图件与source data。'),
 C('预期产出是一个可复核的疾病情境研究：确定CGNL1信号属于何种来源、邻域及候选上游。是否能写出具体机制轴，取决于跨层与独立证据的实际结果。'))

plan('证据不足时的交付路径与解释边界',
 T(['实际结果','继续交付','不能升级的表述'],[
 ['仅组成解释稳定','CGNL1相关血管/肝背景组成与空间分布；完整报告状态阴性。','不称细胞内激活或新内皮亚群。'],
 ['H1成立，H2/调控不支持','空间邻域关联论文，保留上游阴性与公开刺激结果。','不称Notch驱动CGNL1。'],
 ['H2/扰动支持，空间资格不足','调控响应及人体来源研究；明确组织空间证据不足。','不声称特定直接接触界面。'],
 ['位点模型支持，缺少直接结合','候选顺式调控；结合公开扰动作优先排序。','不把motif或SCENIC+边写成直接结合。'],
 ['H1/H2在独立资料均无支持','完整结题与范围收窄，整理已有机制和人体测量差距。','不继续换阈值/面板制造主要阳性。']],[141,207,159]),
 H2('Discussion需集中说明的边界'),
 P('纯生信的最大可达结论由测量决定：空间共邻近、转录程序、染色质可及性和公开实验响应可以形成一致支持，但不能自动证明人体肝转移中的CGNL1必要性、充分性或临床可干预性。血管成熟、屏障完整、免疫外渗与肿瘤入血/定植不是同一终点。'),
 P('原发灶周细胞文献、正常肝/脑内皮参考及其他癌种miRNA证据各自保留背景标签；它们用于提出与反驳假说，不能替代CRLM内皮直接验证。CGNL1与CGN、线性转录本与circRNA分开。'),
 P('本版本完成标准为上述资格、分析及交付任务被真实执行并记录，或因客观缺失形成明确不可评估结论；不以获得指定阳性、某一TF或某一预后方向作为验收要求。'))

def add_refs(pages,selected=None):
 refs=REFS if selected is None else [r for r in REFS if r[0] in selected]
 for i in range(0,len(refs),6):
  b=[]
  for ident,yr,title,url in refs[i:i+6]:
   b.extend([H2(ident+' | '+yr),P(title),S(url)])
  pages.append(('参考文献与原始来源 '+str(i//6+1),b))

EVIDENCE_NOTES=[
 ('D1','上皮细胞','CGNL1与Tiam1/GEF-H1及Rac1/RhoA调节','与CRLM内皮非同一模型'),
 ('D1','培养内皮','ZO-1/JACOP/p114RhoGEF连接定位','定位上游不等于表达上游'),
 ('D1','HUVEC/周细胞及小鼠视网膜','周细胞接触、Notch依赖诱导和血管稳定','未证实CRLM中的直接TF-位点边'),
 ('D1','细胞生化及敲除救援','CGNL1与非肌肉肌球蛋白直接结合','不以MYH转录本替代定位'),
 ('D1','上皮及小鼠肠组织','CGNL1-CAMSAP3和微管组织','不能直接移植到肝内皮'),
 ('D1','小鼠主动脉/bEnd.3','敲除影响NM2连接定位','基础条件下连接外观未明显破坏'),
 ('D1','尿路上皮等','miR-205对CGNL1的3UTR靶向','与内皮miRNA调控不同背景'),
 ('D1','头颈鳞癌细胞','miR-31-3p靶向CGNL1','需成熟链及相同细胞背景验证'),
 ('D1','MDCK','CGN/CGNL1维持GATA4相关连接蛋白表达','不能倒置为GATA4直接调控CGNL1'),
 ('D2','小鼠肝内皮/实验性转移','Notch增强降低黏附和肝转移','不是CGNL1依赖性证据'),
 ('D2','内皮-CRC细胞/动物及组织','内皮可溶Jagged1与癌细胞干性','信号接收者为癌细胞'),
 ('D2','CRLM组织与体内模型','血管共选与抗血管生成反应','不将所有血管归为新生血管'),
 ('D4','病理共识','HGP的组织形态及评分','不能用表达评分替代病理判读'),
 ('D2','CRLM癌细胞','RUNX1/ARP2/3与血管共选','不能当作内皮CGNL1上游'),
 ('D2','原发CRC周细胞','TCF21与血管周基质重塑','原文肝转移周细胞TCF21未检出；留意已发布更正'),
 ('D2','原发CRC周细胞','NKX2-3、血管张力和血流与转移','与R15复用部分单细胞来源'),
 ('D2','肝窦内皮/肝损伤模型','GATA4/MYC/PDGFB及毛细血管化','不证明CGNL1直接上游'),
 ('D3','人跨器官血管','器官特异图谱和空间测量','已有CGNL1空间面板相关展示'),
 ('D3','人CRC及肝转移','单细胞与空间异质性图谱','已覆盖泛空间整合'),
 ('D3','人CRLM纤维环','ISS空间分区及血管相关邻域','定向面板目标覆盖待核'),
 ('D3','人原发CRC/肝转移multiome','癌细胞可塑性、TF与空间邻域','作者RNA已处理；内皮有效供者另算'),
 ('D3','人CRC宏微转移','空间多组学和病灶生态','GEO与发表扩展版范围不同'),
 ('D3','人MSS-CRLM','肝灌流、区域内皮及T细胞排斥','完整数据可得性/CGNL1面板资格待核'),
 ('D2','肝细胞相互作用筛选','宿主PLXNB2与肝定植','不属于CGNL1/内皮操纵'),
 ('D2','HUVEC/HREC及小鼠内皮','动态Notch响应与RND1直接靶点','公开数据可再检CGNL1；当前未查结果'),
 ('D2','HUVEC','NCOR1扰动与Notch/RBPJ相关作用','RNA/ATAC同研究；非CGNL1直接干预'),
 ('D4','空间统计方法','细胞类型特异邻域关联','模型输出不直接证明细胞接触'),
 ('D4','多组学网络方法','TF-元件-基因关系推断','motif/网络不等于实测结合'),
 ('D3','小鼠脑/肺/肝内皮比较','器官相关蛋白表达与定位差异','物种和器官迁移需标注'),
 ('排除直接转用','胰腺癌circRNA','circCGNL1相关机制','非线性CGNL1蛋白的机制证据')]

def write_evidence():
 rows=[]
 for r,n in zip(REFS,EVIDENCE_NOTES):
  row=dict(zip(['id','year_venue','title','url'],r));row.update(dict(zip(['evidence_level','experimental_context','supported_scope','boundary'],n)))
  row['review_depth']='原始来源的摘要/已获取正文段落定向核查；不宣称全文补充材料逐项复核'
  rows.append(row)
 (OUT/'evidence_register.json').write_text(json.dumps({'review_date':'2026-10-05','review_type':'targeted_scoping_review','entries':rows},ensure_ascii=False,indent=2),encoding='utf-8')

def draw_document(name,pages,subtitle):
 p=OUT/name
 c=canvas.Canvas(str(p),pagesize=(W,H));c.setTitle(subtitle);c.setAuthor('CGNL1 project research synthesis')
 qa=[]; text_pages=[]
 for index,(title,blocks) in enumerate(pages,1):
  c.setFillColor(TEAL);c.rect(0,H-10,W,10,fill=1,stroke=0)
  c.setFont('Hei',9);c.setFillColor(GRAY);c.drawString(44,H-35,'CGNL1 | '+subtitle)
  c.setStrokeColor(colors.HexColor('#DAE5E8'));c.line(44,H-46,W-44,H-46)
  y=H-64
  h=Paragraph(esc(title),styles['title']);_,hh=h.wrap(W-88,H);h.drawOn(c,44,y-hh);y-=hh+18
  texts=[title]
  for b in blocks:
   typ=b[0]
   if typ in ('p','small','h'):
    st=styles['h2' if typ=='h' else 'small' if typ=='small' else 'body']
    pa=Paragraph(esc(b[1]),st);_,ph=pa.wrap(W-88,H)
    if y-ph<56:raise ValueError(f'Overflow {name} page {index} {title}: {y-ph}')
    pa.drawOn(c,44,y-ph);y-=ph+(10 if typ=='h' else 9);texts.append(b[1])
   elif typ=='callout':
    pa=Paragraph(esc(b[1]),styles['body']);_,ph=pa.wrap(W-114,H)
    c.setFillColor(colors.HexColor('#EAF4F4'));c.roundRect(44,y-ph-21,W-88,ph+21,5,fill=1,stroke=0)
    pa.drawOn(c,57,y-ph-9);y-=ph+35;texts.append(b[1])
   elif typ=='table':
    _,heads,rows,widths=b
    data=[[Paragraph(esc(s),styles['th']) for s in heads]]+[[Paragraph(esc(str(s)),styles['table']) for s in row] for row in rows]
    tab=Table(data,colWidths=widths or [(W-88)/len(heads)]*len(heads),hAlign='LEFT')
    tab.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),NAVY),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),8),('RIGHTPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#EFF5F6'),colors.white]),('LINEBELOW',(0,0),(-1,-1),0.3,colors.HexColor('#D9E4E8'))]))
    _,th=tab.wrap(W-88,H)
    if y-th<56:raise ValueError(f'Table overflow {name} page {index} {title}: {y-th}')
    tab.drawOn(c,44,y-th);y-=th+12
    texts.extend([' | '.join(heads)]+[' | '.join(map(str,row)) for row in rows])
  qa.append({'page':index,'title':title,'bottom_y':round(y,2)})
  c.setFont('Song',8);c.setFillColor(GRAY);c.drawString(44,31,'2026-10-05 | 研究设计与证据分层 | v2.0')
  c.drawRightString(W-44,31,f'{index} / {len(pages)}')
  c.showPage();text_pages.append('\n'.join(texts))
 c.save()
 (OUT/(p.stem+'.txt')).write_text('\n\n\f\n\n'.join(text_pages),encoding='utf-8')
 return {'file':name,'pages':len(pages),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'layout':qa}

if __name__=='__main__':
 add_refs(AUDIT);add_refs(PLAN)
 results=[draw_document('CGNL1_literature_and_v1.1_review.pdf',AUDIT,'文献调查与原方案改进审阅'),draw_document('CGNL1_research_plan_v2.0.pdf',PLAN,'空间血管情境与上游调控方案')]
 (OUT/'references.json').write_text(json.dumps([dict(zip(['id','year_venue','title','url'],r)) for r in REFS],ensure_ascii=False,indent=2),encoding='utf-8')
 (OUT/'BUILD_QA.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
 write_evidence()
 print(json.dumps([{k:v for k,v in r.items() if k!='layout'} for r in results],ensure_ascii=False))
