from pathlib import Path
import json,hashlib,datetime
import pandas as pd,numpy as np
from reportlab.lib.pagesizes import A4
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle,Image
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.colors import HexColor,white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from xml.sax.saxutils import escape
B=Path(__file__).resolve().parent;R=B/'remote_results/results';O=B/'output/pdf';O.mkdir(parents=True,exist_ok=True)
pdfmetrics.registerFont(TTFont('Chinese',r'C:\Windows\Fonts\simsun.ttc',subfontIndex=0))
body=ParagraphStyle('body',fontName='Chinese',fontSize=10.5,leading=16.5,spaceAfter=11,wordWrap='CJK',textColor=HexColor('#223547'))
title=ParagraphStyle('title',parent=body,fontSize=19,leading=27,spaceAfter=18,textColor=HexColor('#084b61'))
small=ParagraphStyle('small',parent=body,fontSize=8.5,leading=12.5,spaceAfter=7)
th=ParagraphStyle('th',parent=small,textColor=white)
S=[]
def p(t,style=body):S.append(Paragraph(escape(t).replace('\n','<br/>'),style))
def page(t):
 if S:S.append(PageBreak())
 p(t,title)
def table(rows,widths):
 data=[[Paragraph(escape(str(x)),th if i==0 else small) for x in row] for i,row in enumerate(rows)]
 t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#0b5268')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7),('LINEBELOW',(0,0),(-1,-1),.3,HexColor('#dce4e9')),('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f1f6f8')])]))
 S.append(t);S.append(Spacer(1,12))
def img(name,w):
 from PIL import Image as PI
 path=R/name
 with PI.open(path) as im:h=w*im.height/im.width
 S.append(Image(str(path),width=w,height=h));S.append(Spacer(1,10))
def footer(c,d):
 c.setStrokeColor(HexColor('#cad7df'));c.line(44,44,A4[0]-44,44);c.setFont('Chinese',8);c.setFillColor(HexColor('#526b7d'))
 c.drawString(44,30,'CGNL1 | v3 初步数据与方法验证 | 2026-10-06')
 c.drawRightString(A4[0]-44,30,str(d.page))
inp=pd.read_csv(R/'input_qualification.tsv',sep='\t')
comp=pd.read_csv(R/'TECHNICAL_BASELINE_COMPARISON.tsv',sep='\t')
liana=json.loads((R/'LIANA_COMPLETE.json').read_text())
nc=json.loads((R/'NICHECOMPASS_COMPLETE.json').read_text())
qc=json.loads((R/'QUICHE_COMPLETE.json').read_text())
ec=json.loads((R/'ACCEPTED_EC_MEASUREMENT_SUMMARY.json').read_text())
eff=pd.read_csv(R/'mixed_spot_effects_descriptive.tsv',sep='\t')
page('01 新方案已冻结，核心方法完成分级小试')
p('v3.0方案于2026-10-06 22:28（UTC+8）冻结，随后启动新方法小试。本轮完成两项原生方法运行及一项原生模块测试；未把技术运行成功解释为CGNL1机制成立。')
p('方案主线转为血管界面情境、多细胞邻域、连接与免疫黏附读出的分离，以及CGNL1的额外解释信息。旧A05/A06和A04结果完整保留。')
table([['方法','本轮实际状态','决策'],['基础空间/界面','4张患者标识切片、9,353个QC spot；无独立病理掩膜','可作测量；界面模型待标注'],['LIANA','两张切片原生空间模块完成','保留作候选信号测量'],['NicheCompass','3,000 spots，两种种子各30轮训练完成','保留技术候选；未证优于基线'],['QUICHE','官方邻域模块完成；未运行差异丰度','先解决标签与患者重复'],['SMURF','目标HD记录受限，未取得可用输入','本轮不运行'],['FlowSig / OSDR / SpaCEy','条件输入审核；未运行模型','按资格保留，不能记为验证通过']],[110,270,125])
p('主要发现：CGNL1在现有Visium中可测；混合spot小试未显示一致正向的壁细胞邻域关联。完整结果、阴性与输入失败均保留。',small)
page('02 CGNL1可测与未被提供必须分开记录')
p('GSE217414四张切片均在公开表达矩阵中包含CGNL1且有非零计数。下表为UMI≥500且位于组织内的QC结果；重复符号合并后有17,941个基因符号。')
table([['样本','QC spots','CGNL1总计数','阳性spots','阳性比例']]+[[r['sample'],int(r['n_spots_QC']),int(r['CGNL1_counts']),int(r['CGNL1_positive_spots']),f"{100*r['CGNL1_positive_spots']/r['n_spots_QC']:.1f}%"] for _,r in inp.iterrows()],[110,85,105,100,105])
p('“未测/未提供”不是“表达为零”。本轮用CGNL1、JACOP、PCING、ENSG00000128849及Entrez 84952复核名称映射。')
table([['来源','核查事实','用途'],['BJC 2026 ISS','公开探针附表无目标符号/别名和Entrez 84952','组织背景或方法测试'],['GSE335552','四份541-feature矩阵的name/id均无CGNL1；独立panel.json未取得','按公开输入缺少目标测量处理'],['接受内皮缓存','27文库；17文库满足≥20细胞，对应16位患者（两研究各8位）','新面板可测性，尚非新关联结论']],[100,305,100])
p('接受内皮缓存中8个文库的CGNL1测量值为0，已保留。该数字是文库数，不是患者数。重复基因符号先汇总计数再计算分数，避免一个面板成员被重复加权。',small)
page('03 NicheCompass可运行，基线比较未证实优势')
p('按预设几何规则从前两张切片各选择1,500 spots，共812个训练基因、8个先验及8个附加程序。CGNL1及Notch、黏附、历史连接目标面板均排除于训练。两次训练输出均为有限数值的3,000×16表征。')
table([['方法','跨种子ARI','相邻同簇比例','按簇大小期望同簇比例'],['标记-邻域简单基线',f"{comp.iloc[0].cross_seed_ARI:.3f}",f"{comp.iloc[0].neighbor_same_cluster_fraction:.3f}",f"{comp.iloc[0].cluster_size_expected_same_fraction:.3f}"],['NicheCompass小模型',f"{comp.iloc[1].cross_seed_ARI:.3f}",f"{comp.iloc[1].neighbor_same_cluster_fraction:.3f}",f"{comp.iloc[1].cluster_size_expected_same_fraction:.3f}"]],[155,95,115,140])
img('nichecompass_baseline_maps.png',490)
p('NicheCompass产生更连续的空间簇，但简单基线的随机种子稳定性更高。空间平滑与ARI均不等于生物学正确率；目前没有审核后的邻域真值，也没有患者外准确性比较，因此不宣称复杂模型优于基线。图中颜色在不同运行间不对应。',small)
page('04 LIANA通过空间测试，QUICHE仅完成邻域模块')
table([['LIANA输入','保留的预设信号对','打乱位置后的平均绝对分数变化'],['GSM6716963：1,500 spots','8 / 8','0.1254'],['GSM6716964：1,500 spots','3 / 8','0.0968']],[175,125,205])
p('LIANA 1.10.0采用原生SpatialBivariate余弦分数。八对候选来自DLL4/JAG1-Notch、VEGF、PDGF、CXCL12及TGFβ相关配体受体；第二张切片有五对未通过预设表达覆盖。未运行置换P值。')
p('QUICHE使用作者未修改的graph.py中的spatial_niches_khop，在GSM9815454的6,000个QC细胞上运行。100、200、300微米均完成，且20个细胞通过独立数值核对。')
img('quiche_technical_maps.png',490)
p('粗标签有4,906/6,000（81.8%）未定；其余标签包括790个内皮标记候选、131个壁细胞候选等，均未经正式身份审核。三个半径输出完全相同，符合30近邻上限在100微米内已饱和，不能据此称“尺度稳健”。本轮未运行QUICHE全套差异丰度，也未把这份不含CGNL1的数据当目标验证。',small)
page('05 混合spot试算未显示一致正向关联')
p('调整内皮标记量、文库深度及成纤维、上皮、肝细胞标记量后，CGNL1与邻域壁细胞分数的描述性标准化系数在四张切片中呈两正两负。三个尺度下方向保持这一格局。')
wide=eff.pivot(index='sample',columns='radius_um',values='beta_neighbor_mural_adjusted')
table([['样本','β：100 µm','β：200 µm','β：300 µm']]+[[idx]+[f'{r[k]:+.4f}' for k in [100,200,300]] for idx,r in wide.iterrows()],[160,115,115,115])
p('在200微米尺度，CGNL1与邻域T细胞分数的描述性Spearman相关依次为0.0243、0.0303、-0.0005和0.0719。这些数值较小，不能支持普遍强关联。')
p('上述系数来自混合spot测量；没有独立病理分区、细胞归属或患者级正式检验，因此不将其解释为内皮特异效应或机制阴性。本轮没有计算新的患者关联P值/FDR，也没有依据这些效应选择训练样本、算法、尺度或基因面板。')
p('这项结果决定了扩展顺序：先建立可信的病理坐标和细胞归属，随后按v3的预设情境比较开展正式检验。不能通过改变半径或挑选正向切片替换当前结果。')
p('旧证据仍保留：18位患者A05主结果β=-0.1134、P=0.3660、FDR=0.4130；n=8分支P=0.03231、FDR=0.1667仍属探索性。新方案不将旧队列重新称为独立验证。',small)
page('06 当前输入限制决定下一步实施顺序')
p('本轮结果支持继续建设“病理界面+基础空间测量+可解释邻域”的主线。LIANA可用于空间信号测量，NicheCompass保留为与简单基线并列的候选；QUICHE需先获得合格标签和患者重复。')
table([['待解决项','本轮证据','明确处理'],['病理区域','现有小试未取得独立区域掩膜/HGP标签','优先整理作者区域和H&E审核；不从CGNL1反推区域'],['细胞身份','单细胞空间小试81.8%粗标签未定','重新做多谱系、壁细胞/CAF/星状细胞审核'],['HD输入','Zenodo 19043057及最新21431802均为restricted，公开清单无文件','SMURF未运行；待合格公开或获授权输入'],['ISS输入','查看器URL返回HTML，而非HDF5；未用于分析','保留失败记录；不记为数据下载成功'],['方向/动力学/临床模型','尚无合格内皮信号链、动力学假设和足量空间蛋白结局队列','FlowSig、OSDR、SpaCEy保留资格阶段']],[85,225,195])
p('讨论边界：小试的不同种子不是独立生物重复；3,000 spots和6,000细胞不增加患者n。平滑的空间簇不证明真实邻域，空间共现不证明接触、信号传递或T细胞穿越。CGNL1 RNA也不能单独证明蛋白连接定位或屏障功能。')
p('正式模型继续采用用户口径：有效目标原始P<0.05作为名义候选标准；同类家族报告完整FDR，不将FDR作为唯一门槛。无有效重复不作推断检验。')
page('07 方法、直连记录与复现交付')
p('数据和方法：Visium按公开H5计数与条码坐标对齐，以55微米spot直径和scalefactors转换物理尺度；QC为组织内UMI≥500。邻域按切片构建。面板分数为成员log1p(CP10k)均值；目标与训练输入分离。')
p('LIANA：预设八对配体受体、200微米线性邻接核、非零比例≥1%，余弦空间分数；位置打乱仅作技术对照。NicheCompass：0.3.3、CPU、两种固定种子各30轮、8先验+8附加程序。QUICHE：官方源文件SHA256固定，仅调用原生邻域模块；单个样本，不运行差异丰度。')
p('全部新下载走直连：清除代理变量、NO_PROXY=*；requests trust_env=False；curl禁用配置及代理。服务器和本地直连获取的软件与公开输入通过SSH传输。未启用代理，未申请或绕过受限数据权限。')
p('新增细胞坐标：GSM9815454_cells.parquet.gz，1,572,536字节；SHA256为 b76fd61efe1d90c0e7184ce705b3562f8cb45642938405e67fee0ae09c9ed038。GEO官方下载入口成功，gzip/Parquet可读取，80,036个细胞与现有矩阵条码对齐。',small)
p('技术修复记录：LIANA文档与1.10.0命名空间差异已按实际接口修正；隔离环境补齐匹配的PyTorch图扩展；接受内皮缓存中的重复基因符号已合并。原失败日志与修复前测量表保留，未改写旧项目结果。')
p('复现位置：sgs10项目CRC_routes_v1/v3_pilot_20261006；大数据、独立环境与模型权重位于BioData的CGNL1_v3_pilot_20261006。PLAN_FREEZE.json记录方案哈希，方法参数、脚本、退出码、完成JSON、完整TSV和图件均归档。')
p('主要官方来源：QUICHE github.com/jranek/quiche；NicheCompass nichecompass.readthedocs.io；LIANA liana.readthedocs.io；10x Xenium输出文档；NCBI Gene 84952。完整文献和端点定义见v3方案及前轮20项文献矩阵。',small)
out=O/'CGNL1_v3_initial_validation_report.pdf'
SimpleDocTemplate(str(out),pagesize=A4,leftMargin=44,rightMargin=44,topMargin=43,bottomMargin=60).build(S,onFirstPage=footer,onLaterPages=footer)
import pymupdf
doc=pymupdf.open(out)
qa={'pages':len(doc),'characters_per_page':[len(x.get_text()) for x in doc],'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'new_patient_inference_run':False,'visual_QA':'pending'}
(B/'REPORT_BUILD_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
print(json.dumps(qa))
