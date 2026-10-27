import sys, json, pathlib, hashlib, shutil, zipfile, html
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.enums import TA_LEFT
from pypdf import PdfReader

O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2')
D=O/'delivery'; D.mkdir(exist_ok=True)
DATE='2026-09-24'
def bh(p):
    p=np.asarray(p,float); order=np.argsort(p); n=len(p)
    q=np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1]
    out=np.empty(n);out[order]=np.minimum(q,1);return out
checks=[]
for folder in ['bulk_validation','second_validation']:
    for p in sorted((O/'results'/folder).glob('*_all_genes.tsv.gz')):
        t=pd.read_csv(p,sep='\t'); assert np.isfinite(t['P.Value']).all()
        err=float(np.abs(bh(t['P.Value'])-t['adj.P.Val']).max());assert err<1e-10
        checks.append({'check':'genomewide BH','file':str(p.relative_to(O)),'n':len(t),'max_abs_error':err,'pass':True})
    t=pd.read_csv(O/'results'/folder/'candidate_validation_all_models.tsv',sep='\t')
    for name,x in t.groupby('model'):
        assert len(x)==126 and x.candidate.nunique()==126
        err=float(np.nanmax(np.abs(bh(x['P.Value'].fillna(1))-x.FDR_fixed126)))
        assert err<1e-10
        checks.append({'check':'frozen126 BH','file':folder+'/'+name,'n':126,'max_abs_error':err,'pass':True})
frozen=json.loads((O/'PLAN.json').read_text(encoding='utf-8'))['frozen_genes']
held=json.loads((O/'results/cnv/CNV_PLAN.json').read_text())['heldout_genes']
assert all(g.replace('.','-') in held for g in frozen)
checks.append({'check':'126 candidates excluded from RNA-CNV inference','pass':True})
s=pd.read_csv(O/'results/second_validation/design.tsv',sep='\t')
sp=s[s.site.isin(['PT','LM'])&s.strict_colon]
assert sp.patient.nunique()==13 and len(sp)==26 and (sp.groupby('patient').size()==2).all()
checks.append({'check':'strict colon paired sample design','patients':13,'samples':26,'pass':True})
(D/'统计复核.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')

refs=[
 ('R1','GSE50760：18对原发灶—肝转移灶，FPKM','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50760','本次独立bulk验证；逐患者结肠/直肠部位未明确，按CRC报告。'),
 ('R2','GSE14297：18对CRC，严格结肠子集13对','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14297','作者VST与spline标准化微阵列；本次按解剖部位筛出13对。'),
 ('R3','GSE225857：配对单细胞发现队列','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225857','5名明确结肠癌患者，均有术前化疗；本次增加RNA推断CNV敏感性。'),
 ('R4','GSE294559：同细胞RNA与ATAC','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294559','本次只分析P001肿瘤组织与配对类器官；不是原发—肝转移比较。'),
 ('R5','GSE294385：CRC空间数据','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294385','沿用首轮5名患者的前沿比较；CGNL1无一致升高。'),
 ('R6','Cingulin and paracingulin tether myosins-2 to junctions to mechanoregulate the plasma membrane','https://pubmed.ncbi.nlm.nih.gov/37204781/','原始研究：CGNL1/paracingulin与NM2A/NM2B的连接定位和力学调节已有研究；不能将基础结合关系包装为新发现。'),
 ('R7','Paracingulin recruits CAMSAP3 to tight junctions and regulates microtubule and polarized epithelial cell organization','https://pmc.ncbi.nlm.nih.gov/articles/PMC10184829/','原始研究，DOI 10.1242/jcs.260745：CAMSAP3招募及上皮组织功能已有证据；未等于结肠癌肝转移机制。'),
 ('R8','Cgnl1, an endothelial junction complex protein, regulates GTPase mediated angiogenesis','https://pmc.ncbi.nlm.nih.gov/articles/PMC5852532/','原始研究：内皮/血管生成功能，说明CGNL1不是肿瘤上皮特异基因。'),
 ('R9','Identification of Tumor Suppressive Genes Regulated by miR-31-5p and miR-31-3p in Head and Neck Squamous Cell Carcinoma','https://pmc.ncbi.nlm.nih.gov/articles/PMC8227492/','CGNL1已进入癌症研究；不能称为从未研究过的癌症基因。'),
 ('R10','circRNA-CGNL1与胰腺癌的研究','https://pmc.ncbi.nlm.nih.gov/articles/PMC10829403/','研究对象为环状RNA，不能直接转用于线性CGNL1蛋白的作用机制。'),
 ('R11','Tissue-specific transcription reprogramming promotes liver metastasis of colorectal cancer','https://pmc.ncbi.nlm.nih.gov/articles/PMC6951341/','FOXA2/HNF1A与肝样转录/增强子重编程已有功能研究，泛肝适应叙事创新性不足。'),
 ('R12','Prognostic value of brown adipocyte-related genes in colorectal cancer: a multi-omics and Mendelian randomization study','https://pmc.ncbi.nlm.nih.gov/articles/PMC12819916/','已包含CRC中PCSK5敲低及增殖/迁移/侵袭功能验证；削弱PCSK5同类设计新颖性。'),
 ('R13','THBS2-producing matrix CAFs promote colorectal cancer progression ... CD47-MAPK axis','https://www.sciencedirect.com/science/article/pii/S2211124725003262','CRC的CREB3L1—THBS2/基质CAF方向已有直接研究。'),
 ('R14','CRC中SLC9家族研究','https://pubmed.ncbi.nlm.nih.gov/34759967/','SLC9B2并非完全未涉及CRC；本次更主要的问题是第二结肠队列不支持差异。'),
]
(D/'文献定位与创新性核查.json').write_text(json.dumps({'date':DATE,'scope':'聚焦检索，不等同系统综述；未检出不是证明文献不存在。','queries':['CGNL1 colorectal colon liver metastasis','paracingulin cancer metastasis','JACOP cancer','CGNL1 CAMSAP3 myosin','CREB3L1 THBS2 colorectal','PCSK5 colorectal knockdown','SLC9B2 colorectal'],'sources':[dict(id=a,title=b,url=c,interpretation=d) for a,b,c,d in refs]},ensure_ascii=False,indent=2),encoding='utf-8')

pdfmetrics.registerFont(TTFont('YaHei',r'C:\Windows\Fonts\msyh.ttc',subfontIndex=0))
pdfmetrics.registerFont(TTFont('YaHeiBold',r'C:\Windows\Fonts\msyhbd.ttc',subfontIndex=0))
pdfmetrics.registerFontFamily('YaHei',normal='YaHei',bold='YaHeiBold',italic='YaHei',boldItalic='YaHeiBold')
ink=colors.HexColor('#15384A');teal=colors.HexColor('#147D87'); pale=colors.HexColor('#EDF5F5')
styles={
 'title':ParagraphStyle('title',fontName='YaHeiBold',fontSize=23,leading=34,textColor=ink,spaceAfter=20,wordWrap='CJK'),
 'h1':ParagraphStyle('h1',fontName='YaHeiBold',fontSize=17,leading=25,textColor=ink,spaceAfter=13,wordWrap='CJK'),
 'h2':ParagraphStyle('h2',fontName='YaHeiBold',fontSize=11.8,leading=19,textColor=teal,spaceBefore=11,spaceAfter=6,wordWrap='CJK'),
 'body':ParagraphStyle('body',fontName='YaHei',fontSize=10.1,leading=17.2,spaceAfter=8,wordWrap='CJK'),
 'small':ParagraphStyle('small',fontName='YaHei',fontSize=8.4,leading=13.6,spaceAfter=6,wordWrap='CJK'),
 'cell':ParagraphStyle('cell',fontName='YaHei',fontSize=9,leading=14.7,wordWrap='CJK'),
 'quote':ParagraphStyle('quote',fontName='YaHeiBold',fontSize=12,leading=20,textColor=ink,spaceAfter=10,wordWrap='CJK'),
}
story=[];md=[]
def para(text,style='body'):
    story.append(Paragraph(text,styles[style]));md.append(text.replace('<b>','**').replace('</b>','**').replace('<br/>','\n'))
def heading(text):para(text,'h1');md[-1]='## '+text;md.append('')
def sub(text):para(text,'h2');md[-1]='### '+text
def page():story.append(PageBreak());md.append('\n---\n')
def table(rows,widths):
    data=[[Paragraph(str(v),styles['cell']) for v in row] for row in rows]
    t=Table(data,colWidths=widths,hAlign='LEFT',repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),pale),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),8),('RIGHTPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8),('LINEBELOW',(0,0),(-1,0),.8,teal),('LINEBELOW',(0,1),(-1,-1),.3,colors.HexColor('#DBE4E8'))]))
    story.append(t);story.append(Spacer(1,8))
    mr=['| '+' | '.join(str(x).replace('<br/>','；') for x in r)+' |' for r in rows]
    mr.insert(1,'|'+'|'.join(['---']*len(rows[0]))+'|')
    md.append('\n'+'\n'.join(mr)+'\n')
def figure(name,width=505):
    from PIL import Image as PILImage
    p=O/'results/figures'/name;w,h=PILImage.open(p).size
    story.append(Image(str(p),width=width,height=width*h/w));story.append(Spacer(1,6));md.append(f'![{name}](../results/figures/{name})')

para('结肠癌肝转移<br/>CGNL1方向深化立项建议','title')
para('2026年9月24日  |  生信深化预分析 / 客户讨论版','small')
para('建议：有条件推进CGNL1机制方向，先完成细胞来源关卡，再安排机制湿实验。','quote')
para('这轮分析将宽泛的“结肠癌转移多组学”收窄为一个可检验的问题：<b>CGNL1相关的细胞连接—骨架耦联，是否帮助结肠癌细胞在肝内存活与定植？</b> 跨队列表达关联已有支持；连接—骨架耦联和肝内定植作用仍是机制假说。')
table([
 ['评估维度','目前判断','决定项目上限的条件'],
 ['创新性','中等；存在明确的机制切口','避开已知CGNL1结合关系及泛肝样重编程，证明其在结肠癌肝定植中的阶段特异作用。'],
 ['可行性','表达关联较强；细胞来源部分支持','两个独立bulk队列复现；单细胞CNV敏感性保留方向，仍需独立患者来源证据。'],
 ['验证难度','中高','需要定位、两个独立扰动、救援及肝定植模型；基础表达与迁移实验不足以完成主线。'],
 ],[66,139,300])
sub('已经具备的立项依据')
para('18对CRC中16对、另一个队列13对明确结肠癌中12对，肝转移灶CGNL1表达升高；两队列均通过全基因组FDR。逐一去除患者后，配对平均差始终为正。CNV推断支持的肿瘤细胞子集保留相同总体方向。')
sub('需要先解决的关键问题')
para('正常肝组织CGNL1表达较高，肝背景标记校正使bulk统计证据减弱。现有证据支持“CGNL1与肝转移组织差异相关”，尚不能确认全部差异来自癌细胞，也未证明其促进转移。')
sub('如何理解“目标10分”')
para('可以据此启动以约10分期刊为目标的正式研究设计；不能据此承诺发表。项目的文章潜力取决于能否完成“患者来源证据—特定转移表型—结构域/互作救援—体内阶段验证”的机制链。')

page();heading('两个独立队列确认CGNL1的组织表达关联')
para('主分析均以患者配对，比较肝转移灶与原发灶。第二个队列进一步排除了直肠与直乙交界样本，直接支持结肠癌研究范围。[R1–R2]')
table([
 ['指标','GSE50760 / CRC','GSE14297 / 明确结肠癌'],
 ['独立患者对数','18','13'],
 ['肝转移升高比例','16/18','12/13'],
 ['配对效应（95% CI）','+1.031（0.571–1.491）','+1.251（0.894–1.607）'],
 ['效应单位','log2(FPKM+1)差值','作者标准化表达单位差值'],
 ['P / 全基因组FDR','1.48e-4 / 0.00532','3.12e-6 / 0.000997'],
 ['冻结126候选内FDR','0.00222','0.000173'],
 ['逐患者剔除后的均值范围','+0.877至+1.105','+1.160至+1.367'],
 ],[129,188,188])
figure('lead_independent_replication.png')
para('图1｜连线连接同一患者的原发灶与肝转移灶。两平台的表达尺度不同，不能把效应直接相加，也不能将微阵列差值解释为倍数。13对结肠样本来自第二队列18对CRC之内，未重复计作独立队列。','small')
para('加入间质与上皮表达代理评分后，CGNL1仍为正向，两个队列全基因组FDR分别为0.00245与0.000365。这支持其关联不完全由这两个代理评分解释；代理评分不等于精确的细胞比例。')

page();heading('肿瘤细胞中存在CGNL1信号，但来源尚未完全分离')
para('在5名结肠癌患者的单细胞队列中，增加独立于候选基因的RNA推断CNV筛选后，有4对达到每个部位至少20个细胞的预设门槛。[R3]')
table([
 ['分析层级','CGNL1结果','证据定位'],
 ['CNV支持 / 95%阈值','4对：3升1降；中位差+0.900','探索性；患者数小于10'],
 ['CNV支持 / 99%阈值','4对：4升0降；中位差+0.749','更严格阈值敏感性；非DNA确证'],
 ['RNA+ATAC / P001','组织及类器官上皮样细胞均检出RNA与启动子邻近可及性','单患者；多组学可检测性证据'],
 ],[132,192,181])
figure('lead_source_and_boundary.png')
para('图2｜左、中：RNA推断CNV支持的肿瘤细胞伪bulk，单位为log2(CPM+1)；两种阈值全部展示。右：第二队列正常肝（5例）表达高于原发灶和肝转移灶，作为来源解释的关键对照，仅作描述。','small')
para('P001上皮样细胞中，CGNL1 RNA检出比例为组织11.1%、类器官11.6%；启动子邻近ATAC检出比例为15.6%、9.3%。同细胞同时检出分别为9与77个。组织与类器官来自同一患者，不能当成两个患者或转移验证。[R4]')
para('这些结果说明CGNL1可在肿瘤相关上皮细胞中表达，具备开展调控研究的材料基础；它同时表达于成纤维细胞和内皮细胞，尚不能称为肿瘤上皮特异靶点。','small')

page();heading('创新切口应落在肝定植功能及机制条件')
para('文献已建立CGNL1（paracingulin/JACOP）参与细胞连接、肌球蛋白定位和微管组织的基础功能。创新性应来自这些功能在结肠癌肝转移中的具体作用，而非重新证明已知结合关系。')
table([
 ['已知研究','对本项目的意义'],
 ['CGNL1与NM2A/NM2B结合、参与连接与膜力学调节 [R6]','提供肌球蛋白耦联候选机制；直接结合本身不新。'],
 ['CGNL1招募CAMSAP3、组织微管与上皮极性 [R7]','提供微管/极性候选机制；需要筛选与肝定植表型真正相关的分支。'],
 ['血管内皮及头颈癌相关研究 [R8–R9]','不能称CGNL1“从未在癌症中研究”；需区分癌细胞与血管来源。'],
 ['circCGNL1胰腺癌研究 [R10]','环状RNA与线性CGNL1蛋白不是同一研究对象。'],
 ['CRC肝样转录与增强子重编程 [R11]','泛“肝适应/多组学重编程”主线已有先例，需转向明确的功能问题。'],
 ],[226,279])
sub('建议研究问题')
para('<b>CGNL1是否通过维持连接—骨架耦联，增强结肠癌细胞团在肝环境中的存活与定植？</b> 将NM2与CAMSAP3作为可比较的机制分支，依据功能与救援结果收敛到一条主线；不预先宣布两条轴均成立。')
sub('现有生信尚未证明互作轴被激活')
para('在18对队列内，CGNL1与微管组织程序在肝转移灶中的相关系数为0.507，P=0.0319，但相关性家族FDR=0.141；校正间质/上皮评分后FDR=0.310。连接、肌动蛋白—肌球蛋白程序亦无CGNL1校正后显著结果。<b>当前数据不支持把程序协同或CGNL1—MYH9/MYH10轴写成已验证结论。</b>')
para('聚焦检索截至2026年9月24日，未定位到直接完整证明“线性CGNL1蛋白调控结肠癌肝定植”的原始研究。该结论是本次检索范围内的文献定位，不等同于保证无先行研究。','small')

page();heading('按关卡推进，先补齐生信，再投入机制实验')
para('现有跨队列表达关联足以支持进一步研究。湿实验投入应建立在正式生信来源验证之后，避免让客户用实验弥补尚未处理的组织混杂。')
table([
 ['关卡','核心工作与判定','未通过时的处理'],
 ['A / 正式生信优先','补充独立患者的配对单细胞或显微切割/肿瘤上皮数据；核查治疗、部位、环境RNA与双细胞；以患者伪bulk检验方向，达到10例后报告P/FDR。验证癌细胞内方向及来源。','若信号仅由肝细胞、内皮或基质解释，停止“癌细胞CGNL1”主线，重新评估来源。'],
 ['B / 人组织定位','配对原发灶—肝转移灶及邻近肝，采用RNA原位或蛋白共定位，区分肿瘤上皮、内皮和基质。分析单位为患者。','若缺乏肿瘤细胞内信号，不直接进入大规模动物机制实验。'],
 ['C / 特定功能','至少两种独立扰动及野生型救援；比较细胞团/单细胞的连接稳定、失巢凋亡、3D生长与肝环境存活。记录基础增殖以分离一般毒性。','若只有普遍增殖下降而无特定转移表型，降低机制主张或更换候选。'],
 ['D / 机制判别','在确认互作与定位后，比较野生型和不能结合NM2或CAMSAP3的突变体救援，控制表达与定位；增加互作伙伴扰动/上位性。','若突变体与野生型救援相同，不支持该互作分支，不能以共表达替代因果证据。'],
 ['E / 体内阶段','选择可分离肝内早期滞留/存活与后期病灶扩增的模型和时间点；必要时用原位自发转移补充全程证据。','若效应完全由一般生长差异解释，不宣称肝定植特异机制。'],
 ],[69,316,120])
para('湿实验数量与动物/患者样本量应在模型预实验后，以独立重复的效应和变异估算；不按细胞数、视野数或技术重复充当样本量。肝定植模型检验转移的一部分步骤，不能单独替代完整自发转移过程。','small')
sub('面向约10分目标的结果组合')
para('更有竞争力的组合是：独立患者来源明确 + 疾病阶段相关功能 + 可区分机制的救援 + 体内一致性。仅“差异表达、预后模型、敲低后增殖/迁移下降”不足以支撑本方案的创新主张。')

page();heading('关键边界决定叙事与候选优先级')
sub('肝背景问题保留为主要矛盾，不能用措辞消除')
para('GSE50760加入肝背景代理评分后，CGNL1效应由+1.031降为+0.016，P=0.972、全基因组FDR=0.99975；GSE14297严格结肠子集由+1.251降为+0.795，P=0.0837、全基因组FDR=0.624。正常肝组织表达较高，与这一敏感性结果一致。')
para('第一队列肝背景评分与转移部位高度重叠（患者内R²=0.804，VIF=5.11），使模型较难分离组织部位与肝相关表达。该共线性不能证明校正“过度”，也不能证明全部信号来自污染。现阶段应保留两个事实：组织层面差异可复现，癌细胞独立贡献尚未完全确立。')
sub('空间路线暂不并入同一机制链')
para('原发灶空间数据中，CGNL1前沿相对深部的中位差为+0.045，4名可配对患者中2升2降；相对全部其他区域为3/5升高。当前数据不支持一致的原发侵袭前沿激活。[R5] 因而优先检验肝内定植是假说选择，不能称已证明CGNL1只在定植阶段发挥作用。')
table([
 ['候选','第二队列13对结肠结果','当前优先级原因'],
 ['CGNL1','差值+1.251；P=3.12e-6；全基因组FDR=0.000997','主候选；跨队列复现，机制切口可检验；来源关卡尚需完成。'],
 ['CREB3L1','差值−0.032；P=0.170；全基因组FDR=0.499','暂缓；独立结肠复现不足，CRC基质CAF轴已有研究 [R13]。'],
 ['PCSK5','差值−0.165；P=0.225；全基因组FDR=0.569','暂缓；复现不足，CRC敲低功能已有研究且不支持简单“降低促进转移”叙事 [R12]。'],
 ['SLC9B2','差值+0.041；P=0.235；全基因组FDR=0.580','备选；本轮未通过独立结肠队列，不能以早期小样本信号替代。'],
 ],[70,186,249])
para('发现单细胞患者均接受术前化疗，未治疗患者中的普适性仍需独立评估。RNA推断CNV不是DNA测序；4对CNV子集、单患者多组学和5例正常肝对照均保持探索性/描述性定位。研究主张不外推到药物可用性、临床预测工具或所有结肠癌患者。','small')

page();heading('分析方法与可复核设计')
sub('队列、预处理与主要模型')
para('候选家族固定为上一轮126个基因。GSE50760采用作者FPKM，剔除28行无法无歧义恢复的Excel日期化基因标签；过滤为36个配对组织样本中至少6个FPKM≥1，共14,344基因。log2(FPKM+1)后采用患者配对limma-trend，模型为patient+site。')
para('GSE14297使用作者VST及spline标准化矩阵，不再取对数。依据样本解剖部位，排除4名直乙交界与1名直肠患者，主要分析13对结肠样本；18对CRC作为敏感性。探针按官方符号、唯一同义词或NCBI RefSeq映射，多探针无结局选择地求均值；11,511个完整基因进入配对limma。CGNL1对应单一探针3840670。')
sub('组成评分、校正与敏感性')
para('肝背景评分使用ALB/FGA/FGB/FGG/APOA1/APOA2/TTR/HP/AGT/GC/CRP；间质评分使用COL1A1/COL1A2/DCN/LUM；上皮评分使用EPCAM/KRT8/KRT18/KRT19。以可测标记的样本内平均表达计分，并加入患者配对模型。逐患者剔除分析报告配对平均差的变化范围，不称为逐次校正后显著。')
para('每个模型在全部测试基因内单独BH校正，另在冻结126候选家族内校正；未测/未过滤通过的候选以P=1计入校正、结果显示缺失。不混合差异、相关性和敏感性家族。全基因组FDR指所有通过过滤且进入该模型的基因，并非人类所有已注释基因。')
sub('RNA推断CNV与患者伪bulk')
para('infercnvpy 0.6.1按患者分析，采用成纤维与内皮参考，每类至少40个细胞；随机一半训练、一半阈值校准，种子24。7,109个常染色体基因用于推断，排除全部126候选及肝背景标记；窗口100、步长10。以参考各类型95%分位的最大值为主阈值，99%为敏感性；每部位至少20个细胞。细胞按患者及部位汇总，4对保留，不作推断性检验。')
sub('多组学与机制相关性')
para('P001组织/类器官分别进行RNA与ATAC联合QC，RNA UMI≥500、基因≥200、ATAC检出峰≥200、线粒体≤25%；至少检出两个EPCAM/KRT8/KRT18/KRT19定义上皮样。基于NCBI GRCh38.p14基因边界的TSS代理，取上游2 kb/下游500 bp与作者峰相交，按链方向定义启动子邻近窗口；不同样本峰集不进行直接差异检验。')
para('机制相关性固定4候选×3不含候选的程序×2部位，共24项；原始与组成调整相关分别BH校正。同一18对队列内的候选选择后分析定位为探索性。程序基因、全部系数、P及FDR见证据包。软件、模型、完成标记和独立BH复核均已保存。','small')

page();heading('主要数据与文献索引')
para('链接直达数据或原始研究。文献核查用于限定创新性；因候选名称、环状RNA及蛋白机制易混淆，解释栏完整保存在“文献定位与创新性核查.json”。','small')
for a,b,c,d in refs:
    para(f'<b>[{a}]</b> <a href="{c}" color="#147D87">{html.escape(b)}</a>','small')
sub('交付内容与阅读顺序')
para('本PDF用于客户立项沟通；“完整分析审计.md”包含全部CGNL1模型、关键阴性、替代候选及数据处理记录；“CGNL1深化分析_证据包.zip”包含脚本、统计结果、图件和校验清单。原始大文件及计算缓存保留在sgs10项目目录。','small')
para('项目目录：/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2','small')
para('本阶段完成的是深化预分析及有条件立项建议。正式生信来源验证关卡尚未完成；尚未进行本项目湿实验。','small')

def footer(c,doc):
    c.saveState();c.setStrokeColor(colors.HexColor('#D2E1E7'));c.line(45, forty:=40,550,forty)
    c.setFont('YaHei',8);c.setFillColor(ink);c.drawString(45,26,'CGNL1 / 结肠癌肝转移深化立项建议  ·  2026-09-24');c.drawRightString(550,26,str(doc.page));c.restoreState()
pdf=D/'CGNL1结肠癌肝转移_深化立项建议.pdf'
doc=SimpleDocTemplate(str(pdf),pagesize=(595.28,841.89),rightMargin=45,leftMargin=45,topMargin=43,bottomMargin=54,title='CGNL1结肠癌肝转移深化立项建议',author='生信分析')
doc.build(story,onFirstPage=footer,onLaterPages=footer)
reader=PdfReader(str(pdf));print('PDF pages',len(reader.pages));assert len(reader.pages)==8
assert all('CGNL1' in p.extract_text() for p in reader.pages)

def mdtable(t):
    cols=list(t.columns); out=['| '+' | '.join(cols)+' |','|'+'|'.join(['---']*len(cols))+'|']
    for row in t.itertuples(index=False,name=None):
        vals=[f'{v:.8g}' if isinstance(v,(float,np.floating)) else str(v) for v in row]
        out.append('| '+' | '.join(vals)+' |')
    return '\n'.join(out)
audit=['# CGNL1结肠癌肝转移深化预分析：完整审计\n\n日期：2026-09-24\n\n结论：有条件优先推进CGNL1。组织表达关联跨队列复现，但癌细胞独立贡献尚未完全确立，机制仍为假说。不得将本交付理解为全套正式生信及湿实验已完成。\n']
audit+=['## 客户报告正文\n']+md
audit+=['\n## 全部CGNL1模型（无选择省略）\n\n第一队列logFC实际为log2(FPKM+1)差，第二队列logFC字段实际为VST标准化表达差，不能等同于倍数。']
for f in ['bulk_validation','second_validation']:
    t=pd.read_csv(O/'results'/f/'candidate_validation_all_models.tsv',sep='\t');t=t[t.candidate=='CGNL1']
    audit+=['### '+f,mdtable(t[['model','logFC','CI.L','CI.R','P.Value','adj.P.Val','FDR_fixed126']])]
for name,file,query in [
 ('CNV全部预设子集','cnv/candidate_CNV_summary.tsv',"gene=='CGNL1'"),
 ('CNV可用细胞数','cnv/paired_CNV_cell_coverage.tsv',None),
 ('肿瘤/基质来源','mechanism_candidates/source_audit.tsv',"gene=='CGNL1'"),
 ('多组学QC','multiome/QC.tsv',None),
 ('CGNL1多组学检测','multiome/candidate_RNA_ATAC_evidence.tsv',"gene=='CGNL1'"),
 ('CGNL1机制程序相关','mechanism_candidates/mechanism_program_associations.tsv',"gene=='CGNL1'"),
 ('肝背景与部位共线性','bulk_validation/site_overlap_diagnostics.tsv',None),
 ]:
    t=pd.read_csv(O/'results'/file,sep='\t')
    if query:t=t.query(query)
    audit+=['\n## '+name+'\n',mdtable(t)]
audit+=['''
## 数据处理异常与修复记录

1. GSE50760样本标题中metastasized cancer曾导致分类不唯一，配对断言阻止模型执行；最终以组织特征明确PT/LM，36样本/18对断言通过。
2. FPKM不是原始计数，使用log2(FPKM+1)与limma-trend，未使用负二项或将FPKM送入计数模型。
3. 28行Excel日期化基因标签无法无歧义恢复，剔除而未按同名合并；其他23,477个唯一基因可用，表达过滤后14,344个测试。
4. 126候选中117个经符号标点还原后精确映射，6个经NCBI唯一同义词恢复（ADA2/CECR1、ADIRF/C10orf116、ALKAL1/FAM150A、CAVIN1/PTRF、RUBCNL/KIAA0226L、TMEM131L/KIAA0922）；GPAT3、MGAM2、SNHG18缺失或不唯一。所有126保留在校正家族。
5. GSE14297最初只按RefSeq映射而遗漏直接符号CGNL1，属于映射失败，非生物学阴性；最终联合唯一符号、同义词和NCBI RefSeq，CGNL1为单探针3840670。多重前缀等歧义注释不纳入。
6. GSE14297作者描述QC后18,184探针；下载矩阵去全缺失后为18,185行，其中存在部分缺失；最终按映射且完整基因纳入11,511个。完整矩阵与映射保存，不人为按目标数删除任意探针。
7. R读取strict_colon字段得到True/False字符串，初次模型运行中止；最终显式解析True/TRUE，13名/26样本断言通过。排除患者2/4/14/15（直乙交界）和6（直肠）。该子集依解剖学定义，未依据表达选择。
8. RNA-CNV每患者采用独立参考校准；候选及肝标记不参与推断，减少直接循环判定。s0115肝转移仅10/5个细胞达到95%/99%阈值，未达20细胞门槛，两个阈值均排除该患者；不得将其他4对描述为原始5对全部一致。
9. UCSC/EBI参考下载出现TLS连接失败，最终使用官方NCBI GRCh38.p14 GTF。未关闭TLS验证。infercnvpy依赖安装于项目vendor，未修改全局环境。
10. CGNL1机制相关性为候选选择后在同一18对队列内的探索性分析，不是新的独立验证。全部相关性及阴性均保留。

## 证据分级与禁止外推

- 组织差异复现：18对CRC与13对明确结肠癌，两研究独立；后者18对全集包含13对子集，不能重复计数。
- CNV支持的癌细胞内变化：4对探索性、同一发现队列敏感性；不能称独立单细胞外部验证，不能等同DNA恶性确证。
- RNA+ATAC：单患者组织/PDO，检测性与可及性证据；不支持转移差异、差异增强子、转录因子因果或样本间峰可及性倍数结论。
- 空间前沿：无一致变化，不能包装为主线阳性。偏向定植是研究设计假说，不是已证明的阶段特异性。
- 程序相关：CGNL1无校正后显著的连接/肌球蛋白/微管程序协同，不能据已知蛋白互作强行构建已验证轴。
- 肝背景校正阴性与正常肝高表达是改变主张的重要结果，客户正文保留。共线性不能用于取消该阴性。
- 全基因组BH覆盖各模型全部测试基因；固定126的BH另列，未用筛出的3–4个候选重算FDR以获得显著。候选选择属于探索流程，不宣称整个发现/选择/验证流程具有单一预注册确认性错误率。
- 标记评分是表达代理，不是估计比例；未独立精确校正肝细胞、内皮、间皮的组成，不得写为“完全去除所有细胞混杂”。

## 可复核位置

本地：C:/Users/31937/CRC_routes_20260923/deep_v2

服务器：/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2

输入原始文件位于项目raw/metadata；原始发现单细胞及空间材料位于上一轮CRC_routes_v1与上级02_raw。脚本保存明确输入路径。大体积原始数据、CNV数组及vendor未打包，继续保留服务器；交付包包含报告、全部轻量统计结果、关键表达矩阵、计划、脚本、图件和SHA256。复现完整CNV需要服务器原始文件及依赖环境。

## 统计验收

独立Python复算8个全基因组模型BH、8个冻结126家族BH，最大差异均小于1e-10；检查冻结基因排除和13对严格结肠设计通过。复核不代替生物学来源验证。结果见统计复核.json。
''']
audit+=['## 文献完整定位\n']
for a,b,c,d in refs:audit.append(f'- [{a}] [{b}]({c})：{d}')
(D/'完整分析审计.md').write_text('\n\n'.join(audit),encoding='utf-8-sig')
(D/'交付说明.txt').write_text('建议先阅读PDF。统计与阴性详见完整分析审计，代码和数据详见证据包。\n当前状态：深化预分析完成；正式生信来源关卡待完成；湿实验未开展。\n目标约10分是研究设计目标，不是发表承诺。\n',encoding='utf-8-sig')
print(json.dumps({'pdf':str(pdf),'pages':len(reader.pages),'statistical_checks':len(checks),'audit':str(D/'完整分析审计.md')},ensure_ascii=False))
