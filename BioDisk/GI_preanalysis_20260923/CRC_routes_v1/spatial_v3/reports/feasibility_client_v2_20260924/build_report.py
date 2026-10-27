from pathlib import Path
import json,html
from reportlab.platypus import SimpleDocTemplate,Paragraph,Table,TableStyle,Spacer,PageBreak,Image
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image as PI
from pypdf import PdfReader
S=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3');O=S/'reports/feasibility_client_v2_20260924';D=S.parent/'deep_v2'
pdfmetrics.registerFont(TTFont('yh',r'C:\Windows\Fonts\msyh.ttc',subfontIndex=0));pdfmetrics.registerFont(TTFont('yhb',r'C:\Windows\Fonts\msyhbd.ttc',subfontIndex=0));pdfmetrics.registerFontFamily('yh',normal='yh',bold='yhb',italic='yh',boldItalic='yhb')
ink=colors.HexColor('#143A4B');teal=colors.HexColor('#137F85');pale=colors.HexColor('#ECF4F5')
st={k:ParagraphStyle(k,fontName='yhb' if k in ['title','h','sub'] else 'yh',fontSize=size,leading=lead,spaceAfter=gap,textColor=ink if k in ['title','h'] else (teal if k=='sub' else colors.black),wordWrap='CJK') for k,size,lead,gap in [('title',23,34,18),('h',17,25,12),('sub',12,19,7),('body',10.2,17,9),('small',8.5,13.5,6),('cell',9,14.5,0)]}
story=[];md=[]
def p(t,k='body'):story.append(Paragraph(t,st[k]));md.append(t.replace('<b>','**').replace('</b>','**').replace('<br/>','\n'))
def h(t):p(t,'h');md[-1]='## '+t
ndef=0
def sub(t):p(t,'sub');md[-1]='### '+t
ndef=0
def page():story.append(PageBreak());md.append('\n---\n')
def tab(rows,widths):
 t=Table([[Paragraph(str(v),st['cell']) for v in row] for row in rows],colWidths=widths,hAlign='LEFT',repeatRows=1);t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),pale),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),8),('RIGHTPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8),('LINEBELOW',(0,0),(-1,0),.8,teal),('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#D7E3E8'))]));story.append(t);story.append(Spacer(1,10));z=['| '+' | '.join(str(v).replace('<br/>','；') for v in r)+' |' for r in rows];z.insert(1,'|'+'|'.join(['---']*len(rows[0]))+'|');md.append('\n'.join(z))
def fig(path,w=505):
 iw,ih=PI.open(path).size;story.append(Image(str(path),width=w,height=w*ih/iw));story.append(Spacer(1,6));md.append('图件：'+str(path))
p('结肠癌肝转移多组学研究<br/>立项可行性报告','title');p('CGNL1与肝转移局部微环境  |  客户版 · 2026年9月24日','small')
p('<b>立项建议：启动以CGNL1为重点候选的多组学深化分析，围绕肝转移局部微环境中的表达来源与功能状态，建立可衔接湿实验的研究主线。</b>')
p('本项目已经具备三个立项基础：两组独立配对数据重复验证CGNL1升高；空间数据将其定位到近肝区域；单细胞分析提供癌细胞相关的探索性线索。研究已从全基因组筛选推进到有明确候选、有组织定位、有后续验证路径的阶段。')
tab([['评估维度','立项判断','实施依据'],['创新性','具有明确切入点','将CGNL1连接/骨架功能与肝转移局部状态联系，解析特定细胞来源及功能条件。'],['可行性','支持启动深化分析','两独立队列重复验证，已有空间、单细胞与RNA+ATAC数据可继续整合。'],['验证难度','可分阶段实施','先完成生信定位与模型选择，再依次开展定位、扰动救援和体内验证。']],[70,140,295])
sub('为什么值得继续投入')
p('18对CRC中16对、13对严格结肠癌中12对呈升高，全基因组FDR均显著，逐一剔除患者后方向稳定。近肝空间主分析9/11患者升高，P=0.0169。这些结果共同支持对CGNL1开展有针对性的深化研究。')
sub('下一阶段将交付什么')
p('明确主要表达细胞与适用模型；拆分细胞组成和细胞状态的贡献；完成外部空间验证；据此筛选上游条件、关键表型和救援方案。已完成36份外部人类单细胞样本矩阵基础检查，为深化分析提供数据基础。')
sub('面向约10分目标的研究设计')
p('以“临床配对关联 - 空间定位 - 细胞来源与状态 - 条件依赖的功能 - 机制救援及体内验证”组织完整证据链。当前证据支持启动生信阶段，后续实验围绕其收敛的主线开展。')
page();h('两组独立队列确认CGNL1在肝转移组织升高')
p('两个公开队列均采用患者配对分析。GSE14297进一步按解剖部位排除直肠和直乙交界病例，提供明确的结肠癌验证。[1-2]')
tab([['指标','GSE50760 / CRC','GSE14297 / 严格结肠子集'],['患者对数 / 升高','18对 / 16对','13对 / 12对'],['平均差与95%CI','+1.031（0.571-1.491）','+1.251（0.894-1.607）'],['效应单位','log2(FPKM+1)差','作者VST标准化表达差'],['P / 全基因组FDR','1.48e-4 / 0.00532','3.12e-6 / 0.000997'],['固定126候选内FDR','0.00222','0.000173'],['逐患者剔除的均值范围','+0.877至+1.105','+1.160至+1.367']],[125,190,190])
fig(D/'results/figures/lead_independent_replication.png')
p('图1｜每条线连接同一患者的原发灶与肝转移灶。两平台尺度不同，效应不可相加；VST差值不能解释为表达倍数。第二队列18对CRC全集包含上述13对，未重复记作独立验证。','small')
p('加入间质与上皮表达代理评分后，两队列关联仍存在，全基因组FDR分别为0.00245和0.000365。组成代理评分并非精确的细胞比例。','small')
page();h('空间分析将研究问题聚焦到肝-肿瘤邻域')
p('近肝区域CGNL1升高将组织层面的差异进一步定位到局部邻域，为研究肝环境相关状态提供入口。GSE294385共11名患者、24张切片；主分析采用19张非复发切片，按患者比较肿瘤区域近肝≤300微米与远端≥600微米的表达。[4]')
tab([['比较','患者一致性','平均差 / 95%CI','统计证据'],['预设空间主分析','9/11升高','+0.541 / 0.120-0.963','P=0.0169'],['排除混合边缘：5/6邻居为肿瘤','7/10升高','+0.231 / -0.193-0.655','P=0.249'],['更严格：6/6邻居为肿瘤','5/8升高','+0.164 / 仅描述','探索性，无推断检验']],[157,86,150,112])
p('前两项分别为预先指定的主要终点，不另作FDR校正；严格结果属于同队列稳健性分析，不是新的独立验证。效应单位均为log2(CPM+1)差。','small')
import csv
from reportlab.graphics.shapes import Drawing,Line,Rect,String
with (S/'results/screen/patient_effects.tsv').open(encoding='utf-8-sig') as ff:
 rows=[r for r in csv.DictReader(ff,delimiter='\t') if r['gene']=='CGNL1' and r['contrast']=='tumor_near300_vs_far600' and r['method']=='pseudobulk' and r['scope']=='baseline']
rows.sort(key=lambda r:int(r['patient']))
assert len(rows)==11 and sum(float(r['effect'])>0 for r in rows)==9
assert abs(sum(float(r['effect']) for r in rows)/11-.541452)<.00001
dr=Drawing(505,250)
dr.add(String(252,234,'近肝区域较远端升高：9/11位患者',fontName='yhb',fontSize=11,textAnchor='middle',fillColor=ink))
yy=lambda v:52+(v+.4)/2.2*159
for tick in [-.25,0,.5,1,1.5]:
 dr.add(Line(46,yy(tick),489,yy(tick),strokeColor=colors.HexColor('#D7E3E8'),strokeWidth=.5))
 dr.add(String(38,yy(tick)-3,str(tick),fontName='yh',fontSize=8,textAnchor='end'))
for j,r in enumerate(rows):
 val=float(r['effect']);xx=56+j*39
 dr.add(Rect(xx,min(yy(0),yy(val)),25,abs(yy(val)-yy(0)),fillColor=teal if val>0 else colors.HexColor('#899FAB'),strokeColor=None))
 dr.add(String(xx+12.5,34,'P'+r['patient'],fontName='yh',fontSize=8,textAnchor='middle'))
dr.add(String(252,12,'患者级近肝 - 远端差值，log2(CPM+1)',fontName='yh',fontSize=9,textAnchor='middle'))
story.append(dr);story.append(Spacer(1,6));md.append('图2：11名患者空间主分析，近肝-远端差值；逐患者数值见归档patient_effects.tsv。')

p('图2｜每个柱代表一位独立患者的空间主分析效应。严格边缘结果见上表；限于共同10名患者，原始均值+0.470，排除边缘后+0.231。患者8因严格规则下spot不足而排除，效应减弱并非仅由患者数变化造成。','small')
p('200/400微米敏感性支持同向富集（10/11升高，差+0.529，P=0.00219，家族FDR=0.00876）。严格边缘排除后效应减弱且未显著，表明来源与边界混合是下一步必须拆解的问题；这一结果将深化分析聚焦到具体细胞和局部状态。')
page();h('多层数据明确细胞来源验证的研究重点')
p('肝小叶与转移基质均呈稳定的CGNL1富集；RNA-CNV支持的上皮细胞也提供升高线索。多层结果支持以细胞来源为主线，解析哪些细胞贡献组织差异，以及其中是否存在转移相关状态。')
tab([['已完成分析','结果','支持到的层级'],['肝小叶区域 vs 肿瘤区域','11/11更高；差+1.361；P=0.000350，FDR=0.00105','肝区域表达富集'],['转移基质区域 vs 肿瘤区域','11/11更高；差+1.113；P=0.000109，FDR=0.000434','间质区域表达富集'],['空间主距离，组成评分调整','系数+0.182；P=0.289，敏感性FDR=0.861','当前不支持独立剩余关联显著'],['严格边缘后，再作组成调整','系数+0.112；P=0.484，敏感性FDR=1.000','优先开展来源拆分']],[163,215,127])
p('前两项为区域伪bulk差；后两项为spot层面log2(CPM+1)的调整系数，不能与伪bulk效应直接比较大小。来源比较的FDR在4候选×3来源对比的12项家族内计算。','small')
sub('肝相关背景是深化分析的关键变量')
p('GSE50760加入肝表达评分后，效应为+0.016，P=0.972、全基因组FDR=0.99975；严格结肠13对为+0.795，P=0.0837、FDR=0.624。肝评分与组织部位高度重叠（发现队列VIF=5.11），需通过细胞分辨率分析区分组成贡献与状态变化。调整后结果不支持独立于肝评分的显著关联。')
sub('单细胞提供可继续追踪的癌细胞相关线索')
p('发现单细胞队列在排除候选基因参与RNA-CNV推断后，95%参考阈值子集保留4对，CGNL1为3升1降、中位差+0.900；99%阈值为4升0降、中位差+0.749。这一小样本探索性结果为独立恶性细胞验证提供了明确方向。[3]')
p('P001组织与类器官上皮样细胞均检出CGNL1 RNA及启动子邻近ATAC可及性，但只来自同一患者，也不是原发-肝转移比较。这为后续模型筛选及调控分析提供了可检测性依据。[5]')
p('功能轴仍待筛选：当前连接、肌球蛋白及微管程序未见CGNL1校正后显著关联，原发灶侵袭前沿也未一致升高。后续优先寻找肝转移情境中的特定状态与条件。','small')
page();h('创新切口：CGNL1与肝转移局部状态的联系')
p('本项目的创新目标，是在已建立的肝环境研究框架内，解析CGNL1由哪类细胞表达、在何种局部条件下改变，以及这种改变如何影响肝内存活或定植。疾病关联、空间情境与功能验证将共同构成新增证据。')
tab([['已有研究基础','本项目可推进的问题'],['肝样转录与增强子重编程；Cell Research 2020 [6]','以器官环境塑造转移状态为生物学依据，进一步定位CGNL1相关细胞与局部条件。'],['CGNL1连接-肌球蛋白耦联；J Cell Biol 2023 [7]','已有蛋白功能可指导表型与救援设计；检验其在结肠癌肝转移中的具体作用。'],['肝细胞Plexin B2-肿瘤semaphorin-KLF4；Nature 2024 [8]','借鉴宿主-肿瘤互作的验证框架，寻找与CGNL1相关的上游条件和功能输出。'],['肝脂肪变性与转移生长方式；Nature 2026 [9]','将肝环境差异作为可检验条件，确定目标状态出现的生物学情境。'],['肝转移生长模式的空间分区；BJC 2026 [10]','结合空间分区和细胞身份，检验CGNL1相关状态与局部结构的联系。']],[234,271])
sub('拟定课题主线')
p('<b>结肠癌肝转移局部微环境中CGNL1的细胞来源、状态关联及功能机制。</b>')
p('优先利用独立单细胞与空间数据确定表达来源，再根据结果选择癌细胞内功能或宿主细胞互作模型。两类模型服务于同一个研究问题：CGNL1相关状态与肝内定植有何联系。')
p('截至本次聚焦检索，未找到直接完整证明“线性CGNL1蛋白介导结肠癌肝定植”的原始研究。CGNL1已有连接力学及其他癌种研究，本课题的创新定位在特定疾病情境、来源与功能条件；circCGNL1为不同分子形式。[7,11-12]','small')
page();h('实施路径清晰，生信结果直接衔接湿实验')
tab([['阶段','核心工作','可交付结果'],['1 / 来源定位','整合独立单细胞与RNA-CNV身份支持，核查恶性上皮、肝细胞及间质来源。','细胞来源图谱与适用实验模型。'],['2 / 空间状态','完成组成拟合、严格边界与外部空间验证，按患者拆分组成和状态贡献。','优先研究的细胞-邻域关系。'],['3 / 条件筛选','结合调控数据、公开干预及文献，筛选与目标状态相关的上游条件。','候选机制、处理条件与检测指标。'],['4 / 功能验证','生信主线收敛后，开展组织定位、至少两个结肠来源模型、两种独立扰动及野生型救援。','与失巢存活、细胞团或定植相关的功能证据。'],['5 / 机制与体内','依据表型开展互作或结构域救援，并区分肝内早期存活与后期扩增。','从局部状态到功能及机制的证据链。']],[81,296,128])
sub('立项与资源安排')
p('建议现在启动来源定位、空间状态和条件筛选三个生信模块；完成后形成统一的湿实验任务书，明确细胞模型、处理条件、主要终点及救援设计。按阶段交付能使实验投入集中在证据最充分的关系上。')
sub('讨论：关键解释与适用范围')
p('严格边缘和组成调整后未达到显著，当前不支持已确认的癌细胞内在机制；也不能据此认定全部来自混合或认定肝环境诱导。Visium为多细胞测量，RNA-CNV为身份支持，均需与独立来源证据结合。')
p('严格结肠队列支持疾病范围；部分空间病例的结肠/直肠信息待核查，发现单细胞来自术前化疗患者。单患者RNA+ATAC用于可检测性证据。约10分为设计目标，最终竞争力取决于新增机制证据、验证完整性与期刊适配。','small')
page();h('分析方法与统计口径')
sub('独立配对bulk')
p('GSE50760使用作者FPKM，剔除28行无法无歧义恢复的Excel日期化基因标签；36份配对组织中至少6份FPKM≥1的14,344基因进入log2(FPKM+1)配对limma-trend。GSE14297使用作者VST及spline标准化值、不再取对数；13对严格结肠样本为主分析，11,511个完整映射基因进入配对limma。每个模型分别计算全测试基因BH及固定126候选家族BH。')
sub('空间区域与边缘排除')
p('采用作者病理注释；spot须UMI≥500、检出基因≥200。按Visium六角阵列100微米中心间距标定物理坐标。主比较为宏转移肿瘤区域近肝小叶中心≤300微米与远端≥600微米，每组至少20个spot。区域计数汇总后计算log2(CPM+1)差，重复切片先等权汇总至患者。')
p('严格边缘验证要求6个相邻spot均有注释，中心及至少5/6邻居为肿瘤区域，近组限定>100且≤300微米；6/6邻居规则作为敏感性。使用相同QC和患者级汇总。组成模型含肝、间质、上皮、内皮标记评分与log UMI；两批次共有EPCAM/KRT8构成统一上皮代理。')
sub('单细胞与多组学')
p('RNA-CNV按患者采用infercnvpy 0.6.1，以独立留出的成纤维和内皮参考细胞校准95%/99%阈值，推断基因排除126候选及肝标记。每部位至少20个细胞后作患者伪bulk。P001 RNA+ATAC按联合QC及上皮标记检出筛选，采用基因边界TSS代理的上游2 kb/下游500 bp窗口与作者峰相交。')
sub('推断、复核与版本')
p('独立患者不足10例仅报告描述性/探索性结果；达到10例报告效应、95%CI、P及适用FDR。单一预设主终点不作FDR。其他候选、来源比较、距离和组成敏感性分别在固定同类家族内校正。所有关键阴性保留；原始数据、患者映射、方案版本、脚本、退出状态和SHA256均归档。')
p('bulk的8个全基因组模型及8个固定候选家族已独立复算；空间主终点从缓存原始计数重建，7个FDR家族复算通过。严格边缘结果采用新版本目录，未覆盖初筛结果。','small')
page();h('主要数据与文献')
refs=[('1','GSE50760：18对CRC原发与肝转移','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760'),('2','GSE14297：18对CRC，本次严格结肠子集13对','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14297'),('3','GSE225857：发现单细胞队列','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225857'),('4','Spatial Multi-Omics Landscape of Colorectal Cancer Macro- and Micrometastases / GSE294385','https://pmc.ncbi.nlm.nih.gov/articles/PMC13390483/'),('5','GSE294559：本次使用P001组织与类器官RNA+ATAC','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294559'),('6','Tissue-specific transcription reprogramming promotes liver metastasis of colorectal cancer','https://pmc.ncbi.nlm.nih.gov/articles/PMC6951341/'),('7','Cingulin and paracingulin tether myosins-2 to junctions to mechanoregulate the plasma membrane','https://pmc.ncbi.nlm.nih.gov/articles/PMC10202830/'),('8','In vivo interaction screening reveals liver-derived constraints to metastasis','https://www.nature.com/articles/s41586-024-07715-3'),('9','Steatosis shapes prognosis-defining liver metastasis heterogeneity in CRC','https://www.nature.com/articles/s41586-026-10686-2'),('10','Spatial transcriptomics differences of histopathological growth patterns in colorectal cancer liver metastases','https://www.nature.com/articles/s41416-026-03567-y'),('11','Lung-residing metastatic and dormant neuroblastoma cells','https://pubmed.ncbi.nlm.nih.gov/21703429/'),('12','circRNA-CGNL1与胰腺癌的研究（环状RNA）','https://pmc.ncbi.nlm.nih.gov/articles/PMC10829403/')]
for n,t,u in refs:p(f'[{n}] <a href="{u}" color="#137F85">{html.escape(t)}</a>','small')
sub('证据文件与交付范围')
p('详细统计及全部敏感性见既有“CGNL1深化分析_证据包”“CGNL1空间初筛_证据包”，严格边缘结果见后台运行run_20260924T052010Z/V01。本报告整合已完成证据，不把尚未完成的来源模型、外部复现或干预分析写成结果。','small')
p('服务器项目：/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1；本地项目：C:/Users/31937/CRC_routes_20260923。','small')
def footer(c,doc):
 c.saveState();c.setStrokeColor(colors.HexColor('#D3E1E5'));c.line(45,40,550,40);c.setFont('yh',8);c.setFillColor(ink);c.drawString(45,26,'结肠癌肝转移 / CGNL1立项可行性 / 客户版  |  2026-09-24');c.drawRightString(550,26,str(doc.page));c.restoreState()
f=O/'结肠癌肝转移_CGNL1立项可行性报告_客户版_20260924.pdf';doc=SimpleDocTemplate(str(f),pagesize=(595.28,841.89),leftMargin=45,rightMargin=45,topMargin=43,bottomMargin=54,title='结肠癌肝转移CGNL1研究可行性报告',author='生信分析');doc.build(story,onFirstPage=footer,onLaterPages=footer);r=PdfReader(str(f));assert len(r.pages)==8,len(r.pages);(O/'可行性报告_可编辑正文.md').write_text('\n\n'.join(md),encoding='utf-8-sig');(O/'REPORT_STATUS.json').write_text(json.dumps({'pages':8,'completed_evidence_cutoff':'2026-09-24','includes_strict_edge_result':True,'future_queue_results_excluded':True,'render_review':'pending'},indent=2));print('8-page feasibility report created')
