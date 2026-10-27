from pathlib import Path
import re, html, json, hashlib
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4
import pymupdf
from PIL import Image, ImageOps, ImageDraw

root=Path(r'C:/Users/31937/CRC_routes_20260923/pure_bioinformatics_plan_20260929')
out=root/'client_pdf_v1.1'; out.mkdir(exist_ok=True)
qa=out/'qa'; qa.mkdir(exist_ok=True)
original=(root/'CGNL1纯生信论文详细计划.md').read_text(encoding='utf-8-sig')
optional=(root/'可选分支_v1.1.md').read_text(encoding='utf-8-sig')
updated=original.replace('版本 v1.0','版本 v1.1')
updated=updated.replace('## 八 不优先投入的内容',optional+'\n\n## 八 可选范围与不适用的外推')
updated=updated.replace('| 常规多算法风险模型和列线图 | 不默认开展 | 现阶段目标是状态与来源，且独立临床验证不足 |','| 机器学习和临床预测 | 可选B01，数据合格后开展 | 与临床基准比较，嵌套验证；外部临床终点需核实 |')
updated=updated.replace('| 从肿瘤细胞系直接推断内皮药敏 | 不开展 | 测量对象和作用来源不匹配 |','| 药物反应与候选优先级 | 可选B02 | 区分真实反应、模型预测与内皮适配，不将肿瘤药敏直接归为内皮药敏 |\n| MR与共定位 | 可选B03 | 先核实强cis工具、祖源、区域数据及适用结局 |')
(root/'CGNL1纯生信论文详细计划_v1.1.md').write_text(updated,encoding='utf-8-sig')

pdfmetrics.registerFont(TTFont('YaHei',r'C:/Windows/Fonts/msyh.ttc',subfontIndex=0))
pdfmetrics.registerFont(TTFont('YaHeiB',r'C:/Windows/Fonts/msyhbd.ttc',subfontIndex=0))
pdfmetrics.registerFontFamily('YaHei',normal='YaHei',bold='YaHeiB',italic='YaHei',boldItalic='YaHeiB')
navy=colors.HexColor('#14374C'); teal=colors.HexColor('#147C82'); gray=colors.HexColor('#546471'); light=colors.HexColor('#EDF4F6')
styles={
 'body':ParagraphStyle('body',fontName='YaHei',fontSize=10.1,leading=16,spaceAfter=7,textColor=navy,wordWrap='CJK'),
 'small':ParagraphStyle('small',fontName='YaHei',fontSize=8.5,leading=13,spaceAfter=5,textColor=gray,wordWrap='CJK'),
 'h1':ParagraphStyle('h1',fontName='YaHeiB',fontSize=21,leading=30,spaceAfter=15,textColor=navy,wordWrap='CJK'),
 'h2':ParagraphStyle('h2',fontName='YaHeiB',fontSize=13,leading=20,spaceBefore=9,spaceAfter=9,textColor=teal,wordWrap='CJK'),
 'h3':ParagraphStyle('h3',fontName='YaHeiB',fontSize=11,leading=17,spaceBefore=6,spaceAfter=5,textColor=navy,wordWrap='CJK'),
 'table':ParagraphStyle('table',fontName='YaHei',fontSize=9,leading=14,textColor=navy,wordWrap='CJK'),
 'th':ParagraphStyle('th',fontName='YaHeiB',fontSize=9,leading=14,textColor=colors.white,wordWrap='CJK'),
}
story=[]; section_titles=[]
def fmt(s):
 s=s.replace('—','-').replace('–','-').replace('‑','-')
 s=html.escape(s)
 s=re.sub(r'\*\*(.+?)\*\*',r'<b>\1</b>',s)
 s=re.sub(r'`([^`]+)`',r'\1',s)
 s=re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',r'<link href="\2" color="#147C82">\1</link>',s)
 return s

def p(s,kind='body'): return Paragraph(fmt(s),styles[kind])
def add(s,kind='body'): story.append(p(s,kind))
def page(title,kicker='研究方案'):
 if story: story.append(PageBreak())
 section_titles.append(title)
 add(kicker.upper(),'small'); add(title,'h1')
def table(headers,rows,widths):
 cells=[[p(x,'th') for x in headers]]+[[p(str(x),'table') for x in row] for row in rows]
 t=Table(cells,colWidths=widths,hAlign='LEFT',repeatRows=1)
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),navy),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),9),('RIGHTPADDING',(0,0),(-1,-1),9),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8),('ROWBACKGROUNDS',(0,1),(-1,-1),[light,colors.white]),('LINEBELOW',(0,0),(-1,0),0.7,teal)]))
 story.append(t);story.append(Spacer(1,10))

def paragraphs(text):
 for block in re.split(r'\n\s*\n',text.strip()):
  block=block.strip()
  if not block:continue
  if block.startswith('### '):add(block[4:],'h2')
  elif block.startswith('**'):add(block)
  else:add(block.replace('\n',' '))

page('CGNL1结直肠癌肝转移\n纯生信研究计划'.replace('\n','<br/>'),'客户沟通版  |  2026年9月29日  |  v1.1')
# Replace the escaped linebreak for the title only
story[-1]=Paragraph('CGNL1结直肠癌肝转移<br/>纯生信研究计划',styles['h1'])
add('研究目标','h2')
add('以CGNL1相关内皮连接状态为主线，通过独立队列、单细胞和空间转录组，明确其细胞来源、局部组织分布及患者层面的关联，形成不依赖新增实验的研究论文。')
add('项目调整','h2')
add('本阶段将投入重点转向跨队列复现、细胞组成与状态分离，以及大样本组织和临床关联。机器学习、药物反应预测与MR纳入可选分支，按数据适配程度启动。')
table(['工作模块','建议投入','主要价值'],[
['来源与内皮状态复现','30%','明确哪类细胞及状态与CGNL1相关'],['空间定位与组成分离','25%','回答状态出现在哪里及组织混合的贡献'],['大样本组织与临床关联','20%','检验群体一致性和临床背景'],['调控与公开干预证据','10%','补充可检验的调控解释'],['稳健性与复现交付','15%','保证证据独立、统计可靠和结果可追溯']], [160,65,266])
add('建议论文方向','h2')
add('结直肠癌肝转移中CGNL1相关内皮连接状态的跨队列与空间转录组研究。')
add('本文件包含22项主体分析、3项可选分支、主图蓝图与阶段安排。投入比例为研究资源建议；预期产出是待检验的目标，不是预先承诺的阳性结果。','small')

page('现有证据与研究重点','01  项目基础')
table(['证据层级','当前结果','下一阶段重点'],[
['配对组织表达','18对CRC中16对升高；独立严格结肠13对中12对升高，两队列均通过原全基因组FDR。','作为研究起点，复用已完成结果。'],
['肝背景解释','肝背景代理评分调整后，两队列CGNL1效应及统计证据减弱。','确定哪些来源贡献组织信号。'],
['内皮状态','已有外部3名合格患者的去身份重叠连接面板均升高；属于探索性结果。','增加独立有效患者，分离深度和亚型因素。'],
['空间定位','严格相对距离分析：发现10人中7人正向，P=0.0819、FDR=0.164；外部3人中2人正向。','统一端点、组成调整和独立空间验证。'],
['大样本资源','两个bulk队列分别有171和167名患者具备CGNL1测量。','开展固定分数、组成及临床关联。'],
['有效内皮人数','表达汇总合格14人，其中外部9人；患者内状态匹配合格6人，其中外部3人。','不同端点分别计数，不能相互替代。']], [86,225,180])
add('计划中的主要突破点','h2')
add('将组织层面的CGNL1升高，推进为来源清楚、患者级可复核、能在独立组织位置中评价的内皮状态问题。内皮表达汇总和患者内检出匹配回答不同问题，两者均保留。')
add('配对bulk及空间数字来自既有项目分析；本计划尚未产生新增验证结果。已探索过的队列继续登记为已探索资料，新增资料按冻结方案验证。','small')

page('数据资源与验证分工','02  数据配置')
table(['资料','用途','资格与使用规则'],[
['GSE50760 / GSE14297','配对组织发现','保留18对CRC与13对严格结肠的区别和原始单位。'],
['GSE225857','单细胞发现','5名既有患者；不再算成全新验证样本。'],
['GSE231559 / GSE285991','已探索外部支持','保留研究来源、亚型和治疗信息。'],
['GSE245552 / GSE178318','来源补充与资格审核','目前分别3人、0人达到既定内皮表达汇总标准。'],
['GSE234804 / multiome RNA','新增验证候选','先核实患者、内皮覆盖及重复；RNA文件到位不等于ATAC合格。'],
['GSE294385','空间发现','项目肝转移子集24张切片、11名患者；每个端点另算有效n。'],
['Niche-DE / GSE217414 / GSE206552 / GSE254364','外部空间候选','先查原始患者重叠和病理距离覆盖，切片数不等于患者数。'],
['GSE159216','主要组织及临床关联','283份表达、171名患者；患者去重，核实随访起点。'],
['GSE73255','组织程序复现','185份肝转移表达、167名患者；当前矩阵无生存字段。'],
['GSE131418 / 正常肝参考','条件扩展与来源解释','审核患者、组织、处理方式；不直接将样本条目数写成n。']], [157,110,224])
add('研究范围按CRC肝转移表述；明确结肠来源的病例作为预设子集。患者在不同病灶、切片及模态中的重复测量不增加独立生物学样本量。','small')

modules=re.findall(r'(?ms)^### (A\d{2}[^\n]*)\n(.*?)(?=^### A\d{2}|^## 五)',original)
assert len(modules)==22,len(modules)
page_titles=['研究入口与组织证据','细胞来源与内皮亚型','内皮状态主检验与患者内复核','独立验证与可迁移状态分数','空间区域与组成解释','空间邻域与外部复现','大样本组织和临床情境','临床关联与预测增量','调控程序与ATAC证据','通讯候选与公开干预','计算扰动与完整证据汇总']
for k in range(0,22,2):
 page(page_titles[k//2],f'03  主体分析  |  A{k+1:02d}-A{k+2:02d}')
 for title,body in modules[k:k+2]:
  title=title.replace(' P0','  核心').replace(' P1','  条件增强').replace(' P2','  可选探索')
  add(title,'h2')
  # Client-language refinements: implementation history remains in internal plan.
  body=body.replace('既定最少spot数','预设最少spot数').replace('不重跑完整靶点筛选','沿用已核查的靶点选择与统计结果')
  paragraphs(body)
  story.append(Spacer(1,8))

branches=re.findall(r'(?ms)^### (B\d{2}[^\n]*)\n(.*?)(?=^### |\Z)',optional)
assert len(branches)==3
for title,body in branches:
 if title.startswith('B03'):
  first,second=body.split('**结局层级：**',1)
  page('MR与共定位的可行性及设计','04  可选分支  |  B03')
  paragraphs(first)
  page('MR的结局解释与预期交付','04  可选分支  |  B03 续')
  paragraphs('**结局层级：**'+second)
  add('投入建议','h2')
  add('先做1至2个工作日的可行性审核；确认合格工具、结局、区域统计与LD资料后，再投入约3至6个工作日开展正式MR和共定位。无合格输入时，以资格审查结论结束本分支。')
 else:
  page('机器学习与临床预测' if title.startswith('B01') else '药物反应与候选药物优先级','04  可选分支  |  '+title[:3])
  paragraphs(body)

page('论文图件与结果组织','05  图件蓝图')
table(['图件','核心问题','预计内容'],[
['Fig.1','为什么研究CGNL1','流程、队列矩阵、配对组织差异及组成调整效应。'],
['Fig.2','信号来自哪些细胞','来源图谱、患者级表达、内皮亚型和比例。'],
['Fig.3','是否有可复现状态','患者级主关联、匹配复核、技术对照及独立验证。'],
['Fig.4','状态出现在哪里','病理区域、细胞组成、残余空间关联及外部患者效应。'],
['Fig.5','是否有群体和临床关联','两个bulk队列、组成调整、临床分层及生存结果。'],
['Fig.6','哪些解释得到支持','调控程序、合格ATAC/干预证据、跨层证据矩阵。'],
['Opt-Fig.1','预测是否有增量','模型比较、外部性能、校准与特征稳定性。'],
['Opt-Fig.2','哪些药物值得跟进','真实或预测反应、跨来源候选一致性和靶点联系。'],
['Opt-Fig.3','是否有遗传支持','MR森林图、区域共定位及适用的敏感性分析。']], [78,136,277])
add('结果呈现原则','h2')
add('主文目标5至6张图，每张围绕一个问题。可选分支有可靠独立验证时可进入主文，其余放补充；ATAC、干预或生存资料不合格时，删去相应主张或面板。所有图件最终以实际结果为准。')
add('每个分析保留完整效应、患者n和统计表。主要阴性与矛盾结果进入正文；全量QC、全部空间切片及敏感性分析进入补充。')

page('统计设计与证据标准','06  统计方案')
for title,body in [
('患者作为独立样本','细胞、spot、ROI、病灶和同患者多模态属于嵌套测量。小于10名有效独立患者的结果定位为探索性；缺少有效重复、全零或模型不适用时只报告描述性结果。达到10例不等于功效充分。'),
('原始P用于候选筛选','有有效统计检验时，以对应目标效应的原始P<0.05作为名义显著筛选标准。FDR未通过的候选可以进入后续分析，但明确称为名义显著或探索性信号。预设主终点无论阳性或阴性均完整呈现。'),
('完整报告多重检验','效应量、方向、患者一致性、P和适用FDR共同解释。差异基因、通路、细胞比例、相关、TF、药物和遗传分析分别定义同类家族；在全部受检成员内校正。单个预设主要终点无需FDR。'),
('分析前冻结关键选择','冻结基因面板、探针规则、模型、距离、患者纳排和临床终点。多因素分析使用目标效应P，不用其他协变量最小P替代。不按显著性更换参数或截点。'),
('验证和训练严格分离','已探索资料如实登记。新增验证冻结完整流程；机器学习按患者分折，预处理和特征选择在训练折内完成，外部资料不用于反复调参。'),
('可选分支分别评价','预测性能以误差、校准和独立验证衡量；MR工具需可靠遗传关联与强度评估，不能只用P<0.05选工具。计算扰动和药物模型统计明确标为模型证据。'),
('保留不一致和不可估计','区分未测量、测得零、QC失败和模型不可估计。方向相反或依赖少数患者的结果完整保留；不能把无数据写成生物学阴性。')]:
 add(title,'h2');add(body)

page('执行安排与阶段交付','07  实施计划')
table(['阶段','预计工作量','主要交付'],[
['输入与方案锁定','2-3工作日','队列登记、有效n、固定面板、模型和空间资格。'],
['来源与内皮状态','5-8工作日','来源/亚型、主关联、患者内匹配和独立验证。'],
['空间与组成','5-8工作日','全部切片资格、组成解释和外部空间效应。'],
['大样本及临床','4-6工作日','跨队列组织关联、临床字段核查及生存分析。'],
['调控与证据整合','3-5工作日','合格调控线索、稳健性及跨层证据矩阵。'],
['图件与稿件交付','4-6工作日','主图、补充材料、统计表、正文和复现材料。']], [143,96,252])
add('主体工作量','h2')
add('预计23至36个工作日，约5至8个工作周。估算以现有服务器和已下载资料为基础；数据审核完成后，独立分析作业可并行安排。受控数据申请或大规模原始数据整理另行估算。')
add('可选分支的增量安排','h2')
add('建议先分别投入1至2个工作日审核三条分支的可行性，再优先选择最多两条。合格输入到位后，机器学习约3至5日，药物候选约3至5日，MR及共定位约3至6日。可选分支单列增量工作量，不挤占主体方案100%的投入。')
add('关键阶段决策','h2')
add('来源与状态通过后推进空间解释；空间证据不足时保留发现性定位或收窄为组成关联；临床与调控增强不决定主体稿件是否必须等待。若核心状态和外部空间均无支持，依据完整结果重新评估题目。')

page('交付与研究定位','08  交付成果')
for title,body in [
('论文与图表','主文图目标5至6张，补充图和完整统计表，配套可编辑图件与source data。每项分析输出方向、效应、有效患者数及适用P/FDR。'),
('正文和方法','形成研究方法、结果及讨论初稿；标题和摘要依据最终证据调整。研究方法客观描述流程，结果先给结论和证据，样本、因果和临床外推边界集中讨论。'),
('数据与复现','提供队列和患者纳排清单、面板与模型参数、脚本、软件版本和运行说明，保留完整QC和敏感性结果。'),
('可选分支','已启动分支交付资格审核、全部模型或药物/遗传结果和验证图件；未通过输入资格的分支交付可行性结论。机器学习拓展A16，不重复计为两套独立证据。'),
('研究定位','论文重点在CRC肝转移情境下的细胞来源、连接状态、空间关系与跨队列一致性。已知内皮连接和肌球蛋白定位功能作为生物学基础；新增价值由本研究的独立验证和情境解释体现。'),
('结论适用范围','有充分证据的结果使用“显示”或“表明”，探索结果清楚标注。公共数据关联、模型药敏和计算扰动不能替代功能实验；一般CRC遗传风险不能等同肝转移因果机制。投稿定位在主要结果完成后按证据强度和期刊范围确定。')]:
 add(title,'h2');add(body)

page('主要依据','09  参考资料')
refs=[
('GSE159216原始队列研究','https://link.springer.com/article/10.1186/s13073-021-00956-1','171名患者的肝转移表达与临床资料，提供分型、多病灶和随访背景。'),
('GSE73255原始研究','https://pmc.ncbi.nlm.nih.gov/articles/PMC5499209/','肝转移表达及CRIS背景；具体字段以实际下载资料为准。'),
('GSE294385官方记录','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294385','CRC原发及转移空间资料；本项目使用核查后的肝转移子集。'),
('Niche-DE原始研究','https://link.springer.com/article/10.1186/s13059-023-03159-6','空间邻域分析及公开病例来源。'),
('CGNL1内皮肌球蛋白定位研究','https://pubmed.ncbi.nlm.nih.gov/40636102/','既有基础功能依据，提示疾病情境与独立验证的重要性。'),
('TRIPOD+AI','https://www.bmj.com/content/385/bmj-2023-078378','预测模型开发、验证和透明报告框架。'),
('oncoPredict原始方法','https://pmc.ncbi.nlm.nih.gov/articles/PMC8574972/','基于表达与药物筛选资料的预测方法。'),
('药物靶点MR方法研究','https://www.nature.com/articles/s41467-020-16969-0','cis遗传工具与组织情境的设计和解释。'),
('cis MR方法研究','https://www.nature.com/articles/s41467-024-50385-y','基因区域相关工具和复杂信号的分析依据。')]
for i,(title,url,desc) in enumerate(refs,1):
 add(f'{i}. [{title}]({url})','h3');add(desc,'small')
add('现有项目数字依据2026年9月24日至25日的分析收口及数据资格核查；本客户版更新于2026年9月29日。网络资料核对用于方案设计，不代替新分析的实际输入审核。','small')

pdf=out/'CGNL1纯生信研究计划_客户版_v1.1.pdf'
def furniture(c,d):
 c.saveState();W,H=A4
 c.setStrokeColor(colors.HexColor('#D6E3E7'));c.line(52,H-42,W-52,H-42)
 c.setFont('YaHei',8);c.setFillColor(gray);c.drawString(52,H-31,'CGNL1  |  结直肠癌肝转移纯生信研究计划')
 c.drawString(52,29,'客户沟通版  ·  2026.09.29  ·  v1.1');c.drawRightString(W-52,29,f'{d.page}')
 c.restoreState()
doc=SimpleDocTemplate(str(pdf),pagesize=A4,rightMargin=52,leftMargin=52,topMargin=57,bottomMargin=49,title='CGNL1结直肠癌肝转移纯生信研究计划',author='研究方案',pageCompression=1)
doc.build(story,onFirstPage=furniture,onLaterPages=furniture)
D=pymupdf.open(pdf); pages=[]; full=''
for i,pg in enumerate(D):
 txt=pg.get_text();full+=txt
 pix=pg.get_pixmap(matrix=pymupdf.Matrix(1.35,1.35),alpha=False);pix.save(str(qa/f'page_{i+1:02d}.png'))
 blocks=[b for b in pg.get_text('blocks') if len(b)>4 and b[4].strip()]
 pages.append({'page':i+1,'characters':len(txt),'text_top':min(b[1] for b in blocks),'text_bottom':max(b[3] for b in blocks)})
# Contact sheets contain every page for inspection, full resolution pages kept alongside.
for start in range(0,len(D),6):
 board=Image.new('RGB',(1000,3*730),'#cbd5df');draw=ImageDraw.Draw(board)
 for offset in range(min(6,len(D)-start)):
  im=Image.open(qa/f'page_{start+offset+1:02d}.png');im.thumbnail((480,690))
  x=10+(offset%2)*500;y=25+(offset//2)*730
  board.paste(im,(x,y));draw.text((x,y-17),f'Page {start+offset+1}',fill='black')
 board.save(qa/f'contact_{start//6+1}.jpg',quality=92)
checks={'pages':len(D),'planned_sections':len(section_titles),'module_A':len(set(re.findall(r'A\d{2}',full))),'module_B':len(set(re.findall(r'B\d{2}',full))),'bad_glyph': '\ufffd' in full,'internal_paths':bool(re.search(r'C:/Users|/mnt/BioDisk|sgs10',full)),'sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'page_stats':pages}
(qa/'content_check.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in checks.items() if k!='page_stats'},ensure_ascii=False))
print(str(pdf))
