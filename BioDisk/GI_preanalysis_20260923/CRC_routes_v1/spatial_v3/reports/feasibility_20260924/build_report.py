from pathlib import Path
import json,html
from reportlab.platypus import SimpleDocTemplate,Paragraph,Table,TableStyle,Spacer,PageBreak,Image
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image as PI
from pypdf import PdfReader
S=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3');O=S/'reports/feasibility_20260924';D=S.parent/'deep_v2'
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
p('结肠癌肝转移多组学研究<br/>可行性报告','title');p('CGNL1相关空间状态与细胞来源  |  已完成证据截至2026年9月24日','small')
p('<b>建议：有条件推进生信深化，保留CGNL1作为重点候选；完成细胞来源验证后，再确定湿实验主线。</b>')
p('已完成分析建立了可复现的组织表达关联，并发现近肝区域的空间富集。严格排除混合边缘后，效应减弱且未达到显著；因此，目前最值得解决的是CGNL1信号的细胞来源及局部响应，而非直接将其定为癌细胞促转移因子。')
tab([['评估维度','当前判断','立项含义'],['创新性','存在具体机制空间','肝适应及肝-肿瘤空间互作已有研究；CGNL1的疾病情境与功能条件需形成新增证据。'],['可行性','数据资源充足，来源问题未解决','两独立bulk队列复现，空间和单细胞提供进一步检验入口。'],['验证难度','中高','需要细胞定位、特定表型、机制救援与体内阶段验证。']],[70,140,295])
sub('已经完成的工作')
p('两个独立配对bulk队列；单细胞RNA推断CNV敏感性；单患者RNA+ATAC检测证据；24张肝转移空间切片；严格边缘排除验证；候选文献定位及统计复核。另有110份文件下载完成、36份人类单细胞样本完成矩阵基础检查。')
sub('尚未计入完成证据的工作')
p('外部单细胞的细胞身份确认、空间组成拟合、独立外部空间复现及公开干预响应验证。下载和矩阵质控完成不等于这些生物学验证已经完成。')
sub('与约10分发表目标的关系')
p('现有证据支持投入下一阶段验证。若能明确来源、发现具体调控关系并完成功能救援，可形成有竞争力的研究设计；目前不足以承诺发表，也不宜仅据现有结果启动全套机制动物实验。')
page();h('跨队列结果支持CGNL1与肝转移组织差异相关')
p('两个公开队列均采用患者配对分析。GSE14297进一步按解剖部位排除直肠和直乙交界病例，提供明确的结肠癌验证。[1-2]')
tab([['指标','GSE50760 / CRC','GSE14297 / 严格结肠子集'],['患者对数 / 升高','18对 / 16对','13对 / 12对'],['平均差与95%CI','+1.031（0.571-1.491）','+1.251（0.894-1.607）'],['效应单位','log2(FPKM+1)差','作者VST标准化表达差'],['P / 全基因组FDR','1.48e-4 / 0.00532','3.12e-6 / 0.000997'],['固定126候选内FDR','0.00222','0.000173'],['逐患者剔除的均值范围','+0.877至+1.105','+1.160至+1.367']],[125,190,190])
fig(D/'results/figures/lead_independent_replication.png')
p('图1｜每条线连接同一患者的原发灶与肝转移灶。两平台尺度不同，效应不可相加；VST差值不能解释为表达倍数。第二队列18对CRC全集包含上述13对，未重复记作独立验证。','small')
p('加入间质与上皮表达代理评分后，两队列关联仍存在，全基因组FDR分别为0.00245和0.000365。组成代理评分并非精确的细胞比例。','small')
page();h('近肝空间富集存在，但严格边缘排除后证据减弱')
p('GSE294385包含11名患者的24张肝转移切片。本轮主分析使用19张非复发切片，以患者为独立单位；比较肿瘤标注区域中距肝小叶spot中心≤300微米与≥600微米的部分。[4]')
tab([['比较','患者一致性','平均差 / 95%CI','统计证据'],['预设空间主分析','9/11升高','+0.541 / 0.120-0.963','P=0.0169'],['排除混合边缘：5/6邻居为肿瘤','7/10升高','+0.231 / -0.193-0.655','P=0.249'],['更严格：6/6邻居为肿瘤','5/8升高','+0.164 / 仅描述','探索性，无推断检验']],[157,86,150,112])
p('前两项分别为预先指定的主要终点，不另作FDR校正；严格结果属于同队列稳健性分析，不是新的独立验证。效应单位均为log2(CPM+1)差。','small')
fig(O/'strict_matched_patients.png')
p('图2｜限于同一批10名患者，原始空间平均差为+0.470，边缘排除后为+0.231。患者8因严格规则下有效spot不足而未纳入；这项共同患者比较表明，点估计减弱并不完全由患者数变化造成。未将两种效应之间的变化宣称为显著差异。','small')
p('预设200/400微米敏感性仍为正向（10/11升高，差+0.529，P=0.00219，距离敏感性家族FDR=0.00876）；400/800微米仅9名有效患者，保持探索性定位。')
page();h('区域来源证据提示肝与基质成分不可忽略')
tab([['已完成分析','结果','支持到的层级'],['肝小叶区域 vs 肿瘤区域','11/11更高；差+1.361；P=0.000350，FDR=0.00105','肝区域表达富集'],['转移基质区域 vs 肿瘤区域','11/11更高；差+1.113；P=0.000109，FDR=0.000434','间质区域表达富集'],['空间主距离，组成评分调整','系数+0.182；P=0.289，敏感性FDR=0.861','当前不支持独立剩余关联显著'],['严格边缘后，再作组成调整','系数+0.112；P=0.484，敏感性FDR=1.000','来源问题仍未解决']],[163,215,127])
p('前两项为区域伪bulk差；后两项为spot层面log2(CPM+1)的调整系数，不能与伪bulk效应直接比较大小。来源比较的FDR在4候选×3来源对比的12项家族内计算。','small')
sub('bulk肝背景敏感性与空间结果方向一致')
p('GSE50760加入肝表达评分后，效应为+0.016，P=0.972、全基因组FDR=0.99975；严格结肠13对为+0.795，P=0.0837、FDR=0.624。肝评分与组织部位高度重叠，因此结果既不能证明全部来自污染，也不能证明肝环境诱导了癌细胞表达。')
sub('癌细胞相关线索仍值得继续核查')
p('发现单细胞队列在排除候选基因参与RNA-CNV推断后，95%参考阈值子集保留4对，CGNL1为3升1降、中位差+0.900；99%阈值为4升0降、中位差+0.749。均属小样本探索性证据，RNA推断CNV不等于DNA确证。[3]')
p('P001组织与类器官上皮样细胞均检出CGNL1 RNA及启动子邻近ATAC可及性，但只来自同一患者，也不是原发-肝转移比较。这支持可检测性，尚不支持差异增强子或转移诱导机制。[5]')
p('连接、肌球蛋白和微管程序相关性分析尚无CGNL1校正后显著结果。原发灶侵袭前沿亦无一致升高。当前不支持把CGNL1-肌球蛋白/CAMSAP3轴写成已由生信验证的结论。','small')
page();h('创新应落在具体细胞来源与功能条件')
tab([['已有研究','对本项目的定位'],['肝样转录与增强子重编程；Cell Research 2020 [6]','FOXA2/HNF1A等机制已有功能研究，且使用过GSE50760；整体“肝适应”现象不新。'],['CGNL1连接-肌球蛋白耦联；J Cell Biol 2023 [7]','CGNL1与NM2A/NM2B的直接结合及连接力学作用已有证据，不能包装为首次发现。'],['肝细胞Plexin B2-肿瘤semaphorin-KLF4；Nature 2024 [8]','已有空间交界分析和体内功能证明，不能只以“加入空间转录组”作为创新。'],['肝脂肪变性与转移生长方式；Nature 2026 [9]','肝环境塑造转移表型已有机制先例；需找到CGNL1相关的特定条件。'],['肝转移生长模式的空间分区；BJC 2026 [10]','肝侧、肿瘤侧和包膜分区已有研究；可借鉴设计，但需新的生物学问题。']],[234,271])
sub('建议保留的研究问题')
p('<b>CGNL1相关的连接状态，是否与某类肝-肿瘤接触邻域相联系；这一联系由哪类细胞贡献，并通过什么条件影响肝内存活或定植？</b>')
p('若主要来自癌细胞，研究局部信号与癌细胞连接/骨架功能；若主要来自宿主或间质，转向该细胞来源及其与癌细胞的相互作用。两条解释由数据区分，不预先固定促转移方向。')
p('截至本次聚焦检索，未找到直接完整证明“线性CGNL1蛋白介导结肠癌肝定植”的原始研究。但CGNL1已有其他癌种及转移相关报道，不能称为从未研究过的癌症基因；circCGNL1研究亦不能直接替代其蛋白机制。[11-12]','small')
page();h('分阶段推进，把湿实验投入放在关键证据之后')
tab([['阶段','下一步要回答的问题','推进条件'],['1 / 来源验证','去除混合边缘后，在独立恶性细胞身份支持下，CGNL1由哪类细胞表达？','保留未确定身份的细胞；不只用上皮标记判癌。'],['2 / 组成与外部空间','细胞组成能解释多少差异；剩余状态是否与局部邻域相关，能否在外部数据复现？','按患者汇总；不同参考模型及合格外部样本支持同一方向。'],['3 / 上游线索','公开干预或调控数据中，哪些肝相关条件能够改变目标状态？','有可解释的处理与生物学重复；不预定已有轴必然成立。'],['4 / 湿实验定位与功能','临床组织定位、至少两个结肠来源模型、两个独立扰动及野生型救援；区分一般生长与细胞团/失巢存活。','来源与特定功能一致，排除普遍毒性。'],['5 / 机制与体内','用互作/结构域救援区分候选机制，并区分肝内早期存活和后期扩增。','形成可复核机制链；原位自发转移按主张需要补充。']],[81,296,128])
sub('项目范围与主要讨论边界')
p('严格结肠的独立bulk队列提供疾病范围支持；空间队列的逐例结肠/直肠信息尚未明确，因此按CRC支持证据报告。发现单细胞患者均接受术前化疗，不能直接外推未治疗人群。')
p('Visium测量多细胞spot；空间邻近不等于直接接触或因果调控。表达评分和探索性组成拟合也不等于测得真实细胞比例。小样本结果保留效应与一致性，不通过细胞数或spot数增加独立n。')
p('下一阶段若无法建立可靠细胞来源及相应状态联系，应更换候选或调整研究主线。仅完成表达差异、预后模型与一般增殖/迁移实验，不足以支撑本项目面向约10分目标的机制创新。')
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
 c.saveState();c.setStrokeColor(colors.HexColor('#D3E1E5'));c.line(45,40,550,40);c.setFont('yh',8);c.setFillColor(ink);c.drawString(45,26,'结肠癌肝转移 / CGNL1研究可行性  |  2026-09-24');c.drawRightString(550,26,str(doc.page));c.restoreState()
f=O/'结肠癌肝转移_CGNL1研究可行性报告_20260924.pdf';doc=SimpleDocTemplate(str(f),pagesize=(595.28,841.89),leftMargin=45,rightMargin=45,topMargin=43,bottomMargin=54,title='结肠癌肝转移CGNL1研究可行性报告',author='生信分析');doc.build(story,onFirstPage=footer,onLaterPages=footer);r=PdfReader(str(f));assert len(r.pages)==8,len(r.pages);(O/'可行性报告_可编辑正文.md').write_text('\n\n'.join(md),encoding='utf-8-sig');(O/'REPORT_STATUS.json').write_text(json.dumps({'pages':8,'completed_evidence_cutoff':'2026-09-24','includes_strict_edge_result':True,'future_queue_results_excluded':True,'render_review':'pending'},indent=2));print('8-page feasibility report created')
