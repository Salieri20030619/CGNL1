from pathlib import Path
import json,datetime,html,base64,hashlib,itertools,zipfile,shutil
import numpy as np,pandas as pd
O=Path.cwd();D=O/'results';R=O/'report';stamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
t=pd.read_csv(D/'CGNL1_endpoint_summary.tsv',sep='\t');sp=pd.read_csv(D/'target_specificity_summary.tsv',sep='\t');ov=pd.read_csv(D/'overlap_summary.tsv',sep='\t');osp=pd.read_csv(D/'overlap_specificity.tsv',sep='\t');ge=pd.read_csv(D/'geodesic_summary.tsv',sep='\t');gsp=pd.read_csv(D/'geodesic_specificity.tsv',sep='\t');qa=json.loads((O/'NUMERICAL_QA.json').read_text());lit=json.loads((O/'literature/LITERATURE_REGISTER.json').read_text())
def check_summaries(prefix,spfile):
 p=pd.read_csv(D/(prefix+'_patient_effects.tsv'),sep='\t');ss=pd.read_csv(D/(prefix+'_summary.tsv'),sep='\t');sps=pd.read_csv(D/spfile,sep='\t');err=0.
 for _,r in ss.iterrows():
  x=p[(p.endpoint==r.endpoint)&(p.gene==r.gene)].effect.to_numpy();dist=np.array([0.])
  for v in x:dist=np.r_[dist+v,dist-v]
  pp=float(np.mean(abs(dist)>=abs(x.sum())-1e-11));err=max(err,abs(pp-r.P),abs(x.mean()-r['mean']))
 for _,r in sps.iterrows():
  a=p[p.endpoint==r.endpoint].pivot(index='patient',columns='gene',values='effect_SD');x=(a.CGNL1-a.drop(columns='CGNL1').median(axis=1)).to_numpy();dist=np.array([0.])
  for v in x:dist=np.r_[dist+v,dist-v]
  pp=float(np.mean(abs(dist)>=abs(x.sum())-1e-11));err=max(err,abs(pp-r.P),abs(x.mean()-r['mean']))
 return err
additional_qa={'overlap_summary_error':check_summaries('overlap','overlap_specificity.tsv'),'geodesic_summary_error':check_summaries('geodesic','geodesic_specificity.tsv')};assert max(additional_qa.values())<1e-8
(O/'FINAL_SUMMARY_QA.json').write_text(json.dumps(additional_qa,indent=2),encoding='utf-8')
decision={'utc':stamp,'candidate_question_zh':'在结直肠癌肝转移的肿瘤—肝界面，CGNL1与壁细胞邻域的联系是否区别于普通上皮分化标志，并超出局部组织组成？','mainline_status':'not_qualified_as_established_primary_story','candidate_status':'retain_exploratory_relative_spatial_coupling','recommendation':'Do not expand the current dataset into a definitive CGNL1 interface mechanism manuscript. Retain one bounded spatial candidate; advancement requires qualified independent spatial material and robust prespecified replication.','feasibility':'regional mixed-spot modeling feasible; epithelial cell-intrinsic attribution and independent interface replication not qualified','novelty':'conditional potential for target-relative interface coupling; generic zonation, stromal barriers, HGP and NM2 function already studied','reason_not_FDR_only':'main absolute endpoints not nominal; conditional/weighted/geometry robustness and independent input limitations, while nominal specificity signal retained','user_review_pending':True,'this_validation_complete':True,'whole_CGNL1_project_complete':False,'report':'report/CGNL1_spatial_mainline_review.html'}
(O/'MAINLINE_DECISION.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2),encoding='utf-8')
def table(df,cols=None):
 if cols:df=df[cols]
 return '<div class="tablewrap">'+df.to_html(index=False,border=0,escape=True,na_rep='不可评估',float_format=lambda x:f'{x:.6g}')+'</div>'
def img(name,caption):
 data=base64.b64encode((O/'figures'/name).read_bytes()).decode();return f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption)}"><figcaption>{caption}</figcaption></figure>'
def details(title,body):return f'<details><summary>{title}</summary>{body}</details>'
primary=t[t.variant=='primary'][['endpoint','n','mean','mean_SD','positive','negative','P','FDR']].copy();primary['解释']=['调整后的近肝减内部差异' if x=='E1' else '近肝与内部的壁细胞邻域关联差异' for x in primary.endpoint]
specrows=[]
for name,df in [('主模型',sp[sp.variant=='primary']),('邻域组成调整',sp[sp.variant=='neighbor_adjusted']),('100 μm邻域',sp[sp.variant=='radius100']),('300 μm邻域',sp[sp.variant=='radius300']),('250 μm空间块',sp[sp.variant=='blocks250']),('组成平衡权重',osp),('组织路径距离',gsp)]:
 r=df[df.endpoint=='E2'].iloc[0];specrows.append({'分析':name,'患者n':int(r.n),'相对效应':r['mean'],'正向患者':int(r.positive),'P':r.P,'FDR':r.FDR})
specdf=pd.DataFrame(specrows)
elig=pd.read_csv(D/'model_eligibility.tsv',sep='\t');support=pd.read_csv(D/'target_measurement_support.tsv',sep='\t');balance=pd.read_csv(D/'overlap_qualification.tsv',sep='\t');geom=pd.read_csv(D/'geometry_diagnostic.tsv',sep='\t');patient=pd.read_csv(D/'patient_effects.tsv',sep='\t')
source_rows=[['GSE294385现有非复发切片','实测目标、作者病理和坐标可用','19切片/11患者；主模型14切片/10患者','可作区域层面探索；未确证每个spot中的上皮来源'],['GSE217414','4患者有原始矩阵和坐标','同口径独立病理界面掩膜尚未落实','不能用已有一般邻域相关冒充本轮E1/E2外部复现'],['GSE206552 / GSE254364','已有4份/2份空间样本；本轮复核官方登记','本轮未完成患者独立性及逐spot病理界面映射','不能把样本数直接加成新增患者验证'],['源研究Visium HD等','最新概念记录解析到Zenodo 21431802','返回restricted且files为空；单独latest接口也发生TLS失败','不能开展直接高分辨率来源核验；访问限制不是生物阴性'],['2026 HGP ISS','正式论文、作者代码及公开项目可访问','直接.h5ad响应实际为HTML；官方归档入口有效但本轮未完整取得','矩阵/目标覆盖未通过资格；不插补CGNL1以充当实测'],['M242T','邻域技术流程已完成','单患者面板不测CGNL1','仅技术背景，不作目标验证']]
source_table=table(pd.DataFrame(source_rows,columns=['材料','已有条件','本轮资格','结论']))
literature_table=table(pd.DataFrame([{'研究':r['title'],'年份':r['year'],'已经回答':r['prior_result'],'本项目需要增加':r['increment_requirement']} for r in lit['records']]))
body=f'''<header><div class="eyebrow">CGNL1 · 空间主线核定 · 内部审核版 v1</div><h1>界面—壁细胞邻域联系保留候选资格，<br>尚未达到正式主线的证据要求</h1><p class="lead">本轮从“空间转录组能增加什么科学认识”出发，完成原始空间计数重提取、固定问题检验、组成平衡与几何敏感性，并核查相关文献及新材料资格。没有将原来RNA关联的复现失败直接等同为空间路线失败，也没有因新的局部阳性而宣布主线成功。</p><p class="meta">记录时间：{stamp}（UTC；纽约当地2026-10-15）。研究目标：提出值得研究且可验证的问题；本报告供用户审核，不是客户稿或论文定稿。</p></header>
<section class="decision"><h2>核定意见：保留一个问题，暂不推进为已成立的主线</h2><p><strong>建议保留的问题：</strong>在结直肠癌肝转移的肿瘤—肝界面，CGNL1与壁细胞邻域的联系是否区别于普通上皮分化标志，并超出局部组织组成？</p><p>这是有明确空间含义、可被否定的候选问题。研究价值在于区分“一般上皮状态/组织混合”与“界面中的特殊空间关系”。但本轮结果还不能把它写成稳定、独立或细胞自主的CGNL1机制。</p><p><strong>当前建议：</strong>停止在同一批材料上继续扩展故事或挑选有利阈值。保留该候选及全部结果；只有获得合格独立空间材料，能直接解决来源与界面复现问题时，才进入下一阶段。当前没有核定一条已获充分数据支持的正式论文主线。</p></section>
<section><h2>可行性、创新性与实证支持分别判断</h2>{table(pd.DataFrame([['科学问题','成立','位置与邻域关系不同于普通表达相关；有明确的竞争解释和失败标准。'],['区域分析可行性','通过','实测CGNL1、作者病理与坐标支持同患者的区域建模。'],['上皮细胞内来源','未通过','混合spot加代理评分不能确证CGNL1来自哪个细胞。'],['潜在创新性','有条件','靶点相对普通标志物的界面特异空间耦合具有候选增量；一般界面/基质/HGP叙事已被研究。'],['主要绝对效应','未通过本轮名义检验','E1 P=0.8867；E2 P=0.08789。'],['相对普通标志物的线索','名义显著，保留','主模型E2相对效应P=0.02539，FDR=0.05078；不是因FDR未过而剔除。'],['稳健性和独立复现','未通过主线要求','空间块方向不一致；平衡权重/组织路径分析未确证；独立同口径材料资格未落实。']],columns=['维度','判断','依据']))}</section>
<section><h2>从19张切片重建测量，14张切片进入10患者主分析</h2><p>从GSE294385现有非复发肝转移原始矩阵重新读取19张切片，得到11位患者的163,691个QC spot。目标原始计数及UMI与既有经核查记录一致。仅使用作者Layer2/Layer3一致的宏转移肿瘤区域，在同一连通肿瘤区域内比较近肝50–200 μm与内部500–1500 μm。每侧至少20个spot，最终主模型为14张切片、10位患者；P8未达到连通区域内的双侧覆盖要求。</p><p>先按切片估计效应，同患者切片等权平均，再以患者作推断。空间块、spot和细胞都没有替代独立患者数。几何资格筛选没有使用CGNL1表达强弱。</p></section>
<section><h2>主要绝对效应没有支持稳定界面富集或明确的壁细胞耦合增强</h2>{table(primary)}<p>E1是控制深度、上皮身份、成熟分化、成纤维、内皮、肝细胞、壁细胞和缺氧代理及肿瘤连通区域后的近肝减内部差异。平均效应为−0.00373 ln(1+CP10k)，6/10患者正向，P=0.8867。不能据此定义CGNL1为普遍的界面高表达标志。</p><p>E2直接检验“近肝与内部的CGNL1—邻域壁细胞关联是否不同”，而不是分别看两组P值。主尺度200 μm下平均交互效应+0.02903，7/10正向，P=0.08789、FDR=0.17578。该方向值得完整保留，但未达原始P&lt;0.05。</p>{img('01_patient_evidence.png','图1：每个点为一位患者。E1单位为ln(1+CP10k)；E2单位为邻域壁细胞分数增加一个切片内SD时，近肝与内部斜率之差对应的ln(1+CP10k)。右图为CGNL1与六个普通上皮标志的中位标准化效应之差。')}</section>
<section><h2>相对普通上皮标志物的空间差异出现名义信号，但未形成稳定独立证据</h2><p>与CDH1、KRT20、CLDN3、CLDN4、VIL1、CDX2的中位标准化效应比较，CGNL1的E2差异为+0.1495，8/10患者正向，P=0.02539、两终点家族FDR=0.05078。这是预设比较的名义显著结果，可进入候选清单；不是“无信号”。该结果表示相对所选标志物的空间关联不同，不等于CGNL1特异的分子作用。</p>{table(specdf)}<p>300 μm分支仍为正向并通过其敏感性家族FDR，但不能把主尺度从200 μm改成300 μm以宣布成功。邻域组成调整、平衡权重和组织路径分析的相对效应仍偏正，但原始P均未低于0.05；空间块分析则没有保持清楚的E2差异。各分支样本支持也不同。</p></section>
<section><h2>组成平衡和空间尺度检查说明，局部阳性不能替换主结果</h2><p>初次近邻匹配的所有11个足够配对切片均未达到协变量平衡标准，因此没有把匹配结果称为“已经控制组成”。在记录这一结果后，采用不使用目标表达的逻辑回归重叠权重，得到11切片、8患者的可比支持，各合格切片测量协变量均值的最大标准化差异低于6×10⁻⁹。有效样本量不足的切片保留为失败资格，没有放宽阈值。</p>{table(ov[ov.gene=='CGNL1'],['endpoint','n','mean','positive','negative','P','FDR'])}<p>平衡后E1平均+0.00898、6/8正向，E2平均+0.05345、6/8正向。小样本效应与方向保留，但当前不支持将其升级为稳定独立联系。平衡的是所选协变量的一阶矩，不是已控制所有混杂。</p><p>250 μm空间块的E1为7/7正向、P=0.015625、FDR=0.03125。这个阳性完整保留。但同样的切片/患者集合在原spot模型中平均效应仅+0.01606，空间块模型为+0.03866；支持人群变化和聚合尺度均影响结果。E2在该分支为3/7正向、平均−0.01683。不能只选择空间块E1阳性替代原主终点。</p>{img('02_sensitivity.png','图2：不同定义和支持人群下的效应及描述性t区间。权重分析为事后资格/稳健性检查；空间块n=7。图不包含后续组织路径分支，其完整结果列于下表。')}</section>
<section><h2>沿组织网格计算距离后，界面结论仍未获得确认</h2><p>最近肝实质spot距离是解剖接近程度的代理，不能直接当作精确组织界面。检查发现，肝实质标签层级一致性仅改变不到1%的主带spot，但允许路径跨越组织空隙的几何差异在个别切片中改变约31.5%的分带。因此另行记录方案，采用作者组织网格上的相邻110 μm最短路径，保留不可到达区域为不可评估。</p>{table(ge[ge.gene=='CGNL1'],['endpoint','n','mean','positive','negative','P','FDR'])}<p>该检查保留9患者。E1为+0.01249，E2为+0.02838；相对普通标志物E2仍偏正（7/9，P=0.05859、FDR=0.11719），但未达到名义阈值。几何敏感性没有把该方向确认为稳定主线。</p></section>
<section><h2>创新性需要超越已有的界面分区、基质屏障与HGP研究</h2><p>原研究已对同一空间来源的微/宏转移状态和邻域结构进行了深入研究；其他研究也已描述纤维包膜分区、不同生长模式和内皮/免疫区域差异。CGNL1的连接处肌球蛋白组织作用也不是新发现。仅增加一个基因名称、若干热图或另一种空间算法，不足以构成新的主线。</p>{literature_table}<p>本轮候选的潜在增量在于：以明确病理肿瘤区域、患者内对照和普通上皮标志物比较，检验CGNL1在界面中的相对空间关系。检索到的几篇核心空间论文正文未检出CGNL1，不代表附件中没有相关结果，也不证明全球首创。本轮为聚焦创新性核查，不是系统综述；创新性最终仍需要稳定的新证据支撑。</p></section>
<section><h2>独立材料的关键缺口是测量与病理资格，而不是算法数量</h2>{source_table}<p>本轮下载并保存了论文、官方登记、作者代码目录、公开项目和接口响应。HTML查看器响应未当作HDF5输入，部分下载归档未当作完整数据。Zenodo受限记录没有绕过访问控制；当前可用文件清单为空。没有向作者或其他外部对象发送消息。</p></section>
<section><h2>建议的最终主线决策及停止条件</h2><p><strong>候选方向保留，正式主线暂不通过。</strong>不建议继续以“CGNL1在界面升高”“CGNL1驱动肝定植”或“CGNL1普遍调控骨架转录”组织论文。当前最清楚的候选是相对普通标志物的界面—壁细胞邻域联系，但它仍属于单来源探索。</p><p>值得投入的下一步，应当能够带来一个新的决定性信息：独立患者的同口径空间复现，或在实测CGNL1的高分辨率材料中定位到明确细胞来源，并与独立病理界面对应。单纯在已有切片增加模型、分数、距离阈值或通讯网络不会自动填补这一缺口。</p><p>如果合格独立材料到位，可冻结同一个E2相对比较问题及主要协变量，再验证方向、效应、患者一致性和空间尺度。原始P&lt;0.05用于名义筛选，FDR据实报告；不以FDR作为唯一门槛。如果材料资格不成立，或者固定验证仍不支持稳健的空间关系，应关闭这一主张并保存负结果。当前没有证据要求更换靶点，也没有理由强行宣布阳性主线。</p></section>
<section><h2>Methods：固定问题、患者层面推断及完整记录</h2><p>原始输入、样本映射和旧结果在BASELINE_AND_AUTHORIZATION.json及baseline/SOURCE_MANIFEST.json中留存。ANALYSIS_SPECIFICATION.json在新效应计算前冻结，固定CGNL1、六个比较标志、两终点、主距离带、邻域半径、组成面板及推断家族。IMPLEMENTATION_FREEZE.json补充执行细节，仍在新效应计算前保存。</p><p>表达为ln(1+每spot基因计数/全基因UMI×10,000)，零表达保留。空间邻域排除中心本身；孤立spot的邻域测量记缺失。每个模型在切片内控制宏转移连通区域截距；常数协变量剔除后要求设计满秩及至少10残差自由度。切片系数在患者内等权平均，患者间等权总结。</p><p>主要P值为患者效应的双侧精确符号翻转检验；另保留符号检验、描述性t区间及逐患者排除均值范围。少于10患者按探索性呈现，不以spot为独立重复。两项CGNL1终点在同一分支内构成BH家族，普通标志物和靶点相对比较另立家族；预设敏感性分支分别记录。后续组成权重与组织路径检查明确标为事后稳健性分析，不改写初始冻结方案。</p><p>初始匹配、后续重叠权重和模型资格均不使用目标或比较基因的结局值来选择可比样本。失败资格、优化器精度提示、下载/解析失败及修复记录均保留。数值核查使用独立QR求解、残差化公式、患者汇总及手工BH重算。</p></section>
<section><h2>Discussion：研究价值与尚未跨越的边界</h2><p>已有材料可以研究区域层面空间关系，但不能由混合spot确证CGNL1的上皮细胞内来源。壁细胞与其他谱系分数是RNA代理，不等于实测比例或细胞膜接触；病例没有经过本轮病理专家重新标注HGP。最近肝组织及组织路径距离均是操作性定义。</p><p>所有新统计仍使用既有已探索患者；多个方法、半径和图件不会增加独立样本。病例、治疗、测量平台及取样位置的差异限制外推。符号翻转推断依赖患者独立及零假设下效应对称性；少量患者和不同资格门槛限制精度。协变量代理的平衡不证明无未测混杂，调整也不能作为因果中介分析。</p><p>名义相对信号应保留，但当前数据没有把候选问题推进为确认性机制。缺乏同口径外部验证和高分辨率来源资格是实质性缺口；它们不是CGNL1无功能的证据。报告核定的是当前投入和主张的边界，不保证发表可行性或期刊层级。</p></section>
<section><h2>数值验收与可复核附件</h2><p>独立重拟合196项主模型切片×基因系数、核对994项患者汇总和105项主/敏感性统计；最大数值差约3.6×10⁻¹⁵，P值复算一致。重叠权重主效应与交互效应独立重算通过；组织路径分支残差化校验和最终汇总复算通过。旧快照文件哈希和初始方案哈希均保持不变。</p><p>完整运行材料位于服务器：<code>{html.escape(str(O))}</code>。原始输入不复制到本地。下载目录的审核包包含报告、图件、全部汇总/患者效应/资格表、方案、QA和来源登记；全部原始响应、逐spot测量、权重与运行日志保留在服务器。</p>
{details('展开：CGNL1全部预设分支',table(t))}
{details('展开：全部靶点与普通标志物',table(pd.read_csv(D/'all_endpoint_summary.tsv',sep='\t')))}
{details('展开：主模型逐患者效应',table(patient[(patient.variant=='primary')&(patient.gene=='CGNL1')]))}
{details('展开：目标原始计数与检出支持',table(support))}
{details('展开：组成平衡的有效样本量与优化诊断',table(balance))}
{details('展开：几何资格诊断',table(geom))}
{details('展开：所有14张主模型切片的空间图',img('03_all_section_maps.png','图3：展示全部主模型切片，未按效应大小挑选。左图为候选距离带，实际模型另应用连通区域覆盖限制；右图为实测CGNL1，统一使用全体QC spot的99百分位颜色上限。图不替代H&E判读。'))}
{details('展开：数值QA', '<pre>'+html.escape(json.dumps(qa,ensure_ascii=False,indent=2))+'</pre>')}
</section><section><h2>核心原始来源</h2><ol>{''.join('<li><a href="'+html.escape(x['url'])+'">'+html.escape(x['title'])+'</a> ('+str(x['year'])+') — '+html.escape(x['verified'])+'</li>' for x in lit['records'])}</ol><p class="meta">本文为当前授权范围内的研究问题验证报告；未发送外部接收方。用户审核状态：待审核。</p></section>'''
css='''body{font-family:"Microsoft YaHei","Noto Sans CJK SC",Arial,sans-serif;background:#eef2f4;color:#203541;margin:0;line-height:1.8}main{max-width:1120px;margin:32px auto;background:white;padding:44px 52px;box-shadow:0 6px 26px #15384915}h1{font-size:30px;line-height:1.45;color:#123b4a;margin:15px 0}h2{font-size:21px;color:#176879;margin-top:0}p{margin:12px 0}.eyebrow{font-size:13px;letter-spacing:1px;color:#447985}.lead{font-size:17px}.meta{color:#5a6c75;font-size:13px}.decision{background:#fff8e9;border-left:5px solid #c48c30;padding:24px}section{margin-top:32px;border-top:1px solid #dde5e8;padding-top:26px}.tablewrap{overflow-x:auto;margin:18px 0}table{border-collapse:collapse;font-size:13px;min-width:90%}th,td{padding:9px 11px;border-bottom:1px solid #dce5e8;text-align:left;vertical-align:top}th{background:#edf5f6;white-space:nowrap}tr:nth-child(even){background:#f8fafb}img{max-width:100%;height:auto}figure{margin:24px 0}figcaption{font-size:13px;color:#536a73}details{border:1px solid #dbe4e8;padding:12px 16px;margin:14px 0}summary{cursor:pointer;color:#176879;font-weight:600}code{overflow-wrap:anywhere;font-size:12px;background:#eff3f4}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}a{color:#146f90}@media(max-width:700px){main{margin:0;padding:22px}h1{font-size:25px}}@media print{body{background:white}main{box-shadow:none;margin:0;padding:0}section,figure{break-inside:avoid}}'''
doc='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CGNL1 空间主线核定报告</title><style>'+css+'</style></head><body><main>'+body+'</main></body></html>'
(R/'CGNL1_spatial_mainline_review.html').write_text(doc,encoding='utf-8')
abstract=f'''CGNL1 空间主线核定报告：审核摘要
时间：{stamp}

核定结论：保留一个有明确空间含义的候选问题，但当前尚不能核定为已获充分数据支持的正式论文主线。

候选问题：在结直肠癌肝转移的肿瘤—肝界面，CGNL1与壁细胞邻域的联系是否区别于普通上皮分化标志，并超出局部组织组成？

本轮实际完成：重读19张非复发切片原始矩阵（11患者、163,691个QC spot）；冻结方案后，14切片、10患者进入主要比较。进行了200μm主邻域、100/300μm敏感性、距离带改变、普通标志对照、250μm空间块、组成匹配/重叠权重、组织网格路径距离及独立数值验证；保存全部失败和阴性。

主要结果：
1. 调整后的近肝减内部CGNL1效应−0.003728 ln(1+CP10k)，6/10正向，P=0.886719，FDR=0.886719。
2. 界面特异壁细胞邻域关联差异+0.029032，7/10正向，P=0.087891，FDR=0.175781。
3. 相对六个普通上皮标志的中位标准化效应差+0.149535，8/10正向，P=0.025391，FDR=0.050781。名义显著，保留为候选；未因FDR未通过而剔除。
4. 组成平衡保留8患者，相对空间差异仍正（+0.249966；6/8），P=0.054688、FDR=0.109375；组织路径距离保留9患者，相对差异+0.157412、7/9，P=0.058594、FDR=0.117188。均为探索性，未构成确认性复现。
5. 250μm空间块的界面差异7/7正向，P=0.015625、FDR=0.03125；但空间块E2方向不一致，且支持人群和聚合方式改变，不能替换原主终点。

可行性：区域层面问题可计算；混合spot的上皮细胞内来源和独立同口径界面复现尚未通过资格。
创新性：一般界面分区、基质屏障、HGP、内皮免疫区域差异及CGNL1/NM2基础功能已有研究。潜在增量只在明确靶点对照的相对空间关系，尚不能用“正文未检出CGNL1”宣布首创。

材料核查：源研究Zenodo概念记录当前解析到21431802，restricted且files为空；ISS公开查看器和官方归档入口存在，但直接.h5ad响应为HTML，归档只取得不完整前段，未当作数据使用。其他已有空间来源尚无同口径的合格病理界面映射，不能冒充本轮外部验证。

建议：不继续在同一批切片上扩大模型和阈值筛选，也不把当前结果写成CGNL1界面机制。保留这一空间候选；只有合格独立材料能解决来源与界面复现时再推进。原始数据、方法、失败、日志及完整结果均保存在服务器；本轮验证已完成，报告待用户审核，整个项目未宣布完成。

完整报告：CGNL1_spatial_mainline_review.html
服务器目录：{O}
'''
(R/'CGNL1_空间主线核定_审核摘要.txt').write_text(abstract,encoding='utf-8')
print(json.dumps({'report':str(R/'CGNL1_spatial_mainline_review.html'),'bytes':len(doc.encode('utf-8')),'decision':decision['mainline_status'],'additional_qa':additional_qa},ensure_ascii=False))
