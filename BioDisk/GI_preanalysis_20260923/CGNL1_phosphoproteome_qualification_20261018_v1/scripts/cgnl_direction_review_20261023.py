from pathlib import Path
import json,csv,re,hashlib,collections,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
O=P/'audit/CGNL1_DIRECTION_REVIEW_20261023_v1';O.mkdir(exist_ok=True)
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,d):
 f=O/name;f.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert read(f)==d
paths={
 'P01':('胰腺导管腺癌','PDAC_phosphosite_357/ANALYSIS_COMPLETE_357.json','PDAC_phosphosite_357/INDEPENDENT_QA_360.json'),
 'H01':('肝细胞癌','H01_582/RESULT.json','H01_582/INDEPENDENT_QA583.json'),
 'H02':('早期HBV相关肝细胞癌','H02_recovered_590/RESULT.json','H02_recovered_590/INDEPENDENT_QA591.json'),
 'H03':('肝细胞癌','HCC_exploratory_572/RESULT572.json','HCC_exploratory_572/INDEPENDENT_QA574.json'),
 'B03':('肝内胆管癌','B03_recovered_590/RESULT.json','B03_recovered_590/INDEPENDENT_QA591.json'),
 'E01':('食管癌','EC_phosphosite_validation_388/ANALYSIS_COMPLETE_388.json','EC_phosphosite_validation_388/INDEPENDENT_QA_390.json'),
 'G01':('早发胃癌','GASTRIC_validation_415/ANALYSIS_COMPLETE_415.json','GASTRIC_validation_415/INDEPENDENT_QA_419.json'),
 'G02':('胃癌','GC2026_validation_418/ANALYSIS_COMPLETE_418.json','GC2026_validation_418/INDEPENDENT_QA_420.json')}
limits={
 'P01':'同一44例发现队列；CDAP替代处理不是独立患者验证；源观测值/处理来源边界保留。',
 'H01':'S283和S486对额外中心化方向敏感；HCC跨研究患者重叠未排除。',
 'H02':'53对合格但位点完整数不同；缺失/MBR和定位证据边界；HCC重叠未排除。',
 'H03':'GRR插补及观测值掩码未解决；逐位点定位信息受限；HCC重叠未排除。',
 'B03':'单一iCCA配对队列；S283统计强度受处理影响；跨源尺度不可直接合并。',
 'E01':'S296-S299区域合并特征不能拆为独立位点重复证据；跨癌种而非PDAC内验证。',
 'G01':'技术通道已按患者聚合；当前合格位点验证未见原始P<0.05。',
 'G02':'当前合格位点验证未见原始P<0.05；S708仅9对描述。'}
baselines={};records=[];originals={};qas={}
def pick(x,*keys):
 for k in keys:
  if k in x:return x[k]
 return None
for cohort,(disease,rel,qa) in paths.items():
 f=P/'results'/rel;baselines[str(f)]=sha(f);d=read(f);originals[cohort]=d
 qf=P/'results'/qa;baselines[str(qf)]=sha(qf);qas[cohort]=read(qf)
 for i,x in enumerate(d.get('contrasts',d.get('results',[]))):
  site=x.get('site')
  if site and '|' in site:site=site.split('|')[2]
  if site is None:site='/'.join(x['feature_labels'])
  n=pick(x,'n','n_pairs','n_patients');pv=pick(x,'P','p');qv=pick(x,'FDR','fdr','FDR_within_GC','FDR_within_GC2026')
  record=dict(cohort=cohort,disease=disease,site=site,n=n,
    adjusted_mean=pick(x,'mean_adjusted_delta','mean'),site_mean=pick(x,'mean_unadjusted_site_delta','mean_raw_site_delta','mean_site','mean_delta_site'),
    protein_mean=pick(x,'mean_total_protein_delta_same_pairs','mean_total_protein_delta','mean_protein_delta','mean_protein','mean_delta_protein'),
    up=pick(x,'up','positive'),down=pick(x,'down','negative'),P=pv,FDR=qv,
    alternative_centering_mean=pick(x,'uncentered_mean','sensitivity_mean'),alternative_centering_P=pick(x,'uncentered_P','sensitivity_p'),
    status='不可评估' if n==0 else ('按原冻结规格仅描述' if pv is None else ('原始P达到0.05' if pv<0.05 else '当前有效检验未达原始P阈值')),
    original_source=str(f),source_row=i,measurement='源分析定义下的肿瘤-对照位点变化减总CGNL1蛋白变化；不同来源不直接合并效应',limits=limits[cohort])
  records.append(record)
assert len(records)==81 and all(x['n']>=0 for x in records)
save('ALL_SITE_EVIDENCE_LEDGER.json',records)
save('SOURCE_RESULTS_SNAPSHOT.json',originals);save('NUMERICAL_QA_RECORDS.json',qas)
with (O/'ALL_SITE_EVIDENCE_LEDGER.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
common={}
for cohort,rel,idkey in [('P01','PDAC_phosphosite_357/patient_site_values.tsv','case'),('H01','H01_582/PATIENT_VALUES.json','case'),('H02','H02_recovered_590/PATIENT_VALUES.json','case'),('H03','HCC_exploratory_572/patient_site_values572.json','patient'),('B03','B03_recovered_590/PATIENT_VALUES.json','case')]:
 f=P/'results'/rel;baselines[str(f)]=sha(f);ss=collections.defaultdict(set)
 if f.suffix=='.tsv':
  with f.open(encoding='utf-8') as h:a=list(csv.DictReader(h,delimiter='\t'))
 else:a=read(f)
 for x in a:ss[x['site'].split('|')[2] if '|' in x['site'] else x['site']].add(x[idkey])
 common[cohort]={'+'.join(pair):len(set.intersection(*(ss[s] for s in pair))) for pair in [('S198','S297'),('S283','S297'),('S283','S708'),('S283','S486'),('S198','S297','S708')]}
save('COMMON_PATIENT_FEASIBILITY.json',common)
registry=read(P/'metadata/PROJECT_CHECKPOINT_595.json')['registry'];save('26_SOURCE_REGISTRY_SNAPSHOT.json',registry)

candidates=[
 dict(id='C01',title='PDAC中CGNL1多位点相对变化是否反映功能状态重排',priority='优先进行可行性分析',
  hypothesis='疾病组织中CGNL1位点不按总蛋白一致变化，而呈现可重复的位点组合。',
  evidence=['P01 S198 n35 Δ-0.263659 P0.00173010 FDR0.00692038','P01 S297 n35 Δ+0.606412 P0.000283359 FDR0.00212288','P01 S708 n21 Δ+0.628323 P0.000353813 FDR0.00212288','P01 S283 n44 Δ+0.276201 P0.00758637 FDR0.0227591'],
  counterevidence=['仍只有一个PDAC独立发现队列','CDAP重处理支持S283/S708但不是独立生物学验证','胃癌/食管癌不能支持现阶段泛消化道统一结论','相邻S296-S299定位需谨慎'],
  immediate_test='先固定S198-S297及S283-S708两个候选位点对，用共同患者的Δ(siteA-siteB)检验；共有28及21对。三位点共同患者仅14对，不建立高维评分。该新对比尚未计算。',
  rationale='同患者两位点差值可代数抵消共同总蛋白分母和同一组学全局中心化项；不消除肽段特异响应、定位或缺失偏差，不是占有率。',
  go='定位合格、共同患者对比稳定且获得独立PDAC材料支持',stop='差异由少量异常/缺失或定位不清驱动，或仅能复现总蛋白分母变化',
  resources='现有患者级数据可完成探索，独立PDAC位点+总蛋白配对材料仍待取得；机制需实验',novelty='组合所对应的功能状态是可检验问题；差异表或事后组合本身不构成机制创新'),
 dict(id='C02',title='肝内胆管癌CGNL1 S486相对总蛋白降低',priority='优先进行可行性分析',
  hypothesis='iCCA中S486变化超出总CGNL1变化，并与定义明确的细胞状态相关。',
  evidence=['B03 n60，42降18升，Δ-0.319315，P0.000120399，FDR0.000240798','不额外中心化Δ-0.396850，P0.00373336，同方向且仍达原始P阈值','B02现有143例肿瘤S486与总蛋白完整值，可用于后续受限背景问题'],
  counterevidence=['B03只有单一配对发现队列','H01 S486方向随中心化改变，不能当稳定跨病种确认','P01 S486 n26 P0.920455；E01 n25 P0.140915，不能称泛癌通用','B02仅肿瘤，不能验证肿瘤-正常差异'],
  immediate_test='先检查B03的S486定位、位点/蛋白分量、缺失和个体影响；与S283共同患者60对，可检验区别于整体CGNL1变化的位点模式。之后确定具体细胞表型再使用B02，不能直接再筛所有临床终点。',
  go='同癌种独立配对或细胞分辨材料确认位点变化',stop='源信号主要由细胞组成/缺失或定位不清解释，且无独立确认',
  resources='现有B03可做来源/稳健性核查；独立配对材料待补；B02受限原始谱图不是当前必需下载项',novelty='疾病背景下的位点功能可形成研究问题；仅再报告S486差异不保证创新'),
 dict(id='C03',title='HCC中CGNL1位点模式为何在队列间不同',priority='条件性候选，先判定测量可比性',
  hypothesis='观察到的队列差异由测量处理或疾病背景造成；需要先区分二者。',
  evidence=['H02 S283 n38 Δ-1.483808 P1.89505e-5 FDR3.79010e-5；S297 n21 Δ-1.645554 P0.000161171 FDR0.000161171','H01 S297 n34 Δ+0.800823 P3.74042e-7 FDR7.48085e-7，未额外中心化仍同方向','H03 S283 n118 Δ-1.295436 P3.71096e-9 FDR1.85548e-8，但插补边界未解'],
  counterevidence=['H01 S283方向翻转；H03 GRR插补与定位缺口','HCC患者重叠尚未排除','H03 S297主分析P0.134890，不是可靠同向重复','不能把H02早期选择直接解释为分期机制'],
  immediate_test='先统一估计量、配对及定位规则，再在共同患者比较S283-S297，H01/H02/H03分别34/17/27对。只有差异在可比处理中保留且有临床分层资料，才检验预先定义的背景交互；不比较各队列P值代替交互检验。',
  go='排除主要技术来源后存在可重复的背景依赖位点模式',stop='差异随测量口径统一消失或来源无法核对',
  resources='现有数据可先做技术归因；临床背景和独立性需外部信息',novelty='若能解释位点状态转换才有疾病意义；技术处理差异本身不是新机制'),
 dict(id='C04',title='非肿瘤炎症/屏障损伤中的CGNL1磷酸化',priority='储备方向，先过材料资格门槛',
  hypothesis='炎症/应激时的CGNL1位点状态是否不同于单纯蛋白丢失，并影响屏障修复。',
  evidence=['2026 VEO-IBD会议摘要报告CGNL1变异/丢失及肠道表型','CGNL1与微管/连接组织的直接蛋白功能研究提供可操作读出'],
  counterevidence=['现有项目配对磷酸化证据集中于肿瘤，不能外推炎症疾病','会议摘要不是本项目独立原始磷酸化数据','CGNL1缺失-屏障的一般命题已有披露，不是新的创新点'],
  immediate_test='先确认是否能获得同样本CGNL1位点、总蛋白和独立生物学重复的炎症/对照或刺激/恢复数据；不预先指定S283，不用RNA替代磷酸化。',
  go='测量资格通过，并能提出区别于已知CGNL1缺失效应的位点调控问题',stop='只有RNA/总蛋白或没有有效重复',
  resources='当前不能承诺有合格数据，需要定向来源检查；机制实验另需资源',novelty='位点依赖的可逆调控而非简单重复蛋白缺失结论')]
save('CANDIDATE_DIRECTIONS.json',candidates)
lit=[dict(id='L1',url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11513637/',type='primary_biochemical_study',finding='CGNL1微管结合的250-420区域及365-377片段证据；未证明S283或S297的位点依赖功能。'),
 dict(id='L2',url='https://pmc.ncbi.nlm.nih.gov/articles/PMC10184829/',type='primary_cell_and_mouse_study',finding='CGNL1-CAMSAP3与微管/上皮极性；不等于消化癌症磷酸位点机制。'),
 dict(id='L3',url='https://academic.oup.com/ibdjournal/article/32/Supplement_1/izag006.073/8436212',type='conference_abstract',finding='CGNL1丢失/VEO-IBD先前披露；不是磷酸化位点验证。')]
save('LITERATURE_DECISION_REGISTER.json',dict(accessed_date='2026-10-23',records=lit,exhaustive_novelty_review=False,
 exclusions=['circCGNL1不能替代线性CGNL1','CGN的AMPK磷酸化机制不能移植为CGNL1已有证据','仅motif/数据库记录不能指定上游激酶或证明功能']))
summary='''CGNL1磷酸化研究回顾与候选方向（2026-10-23）

建议：暂不确定唯一主线，保留三个现有数据候选和一个非肿瘤储备方向。优先评估PDAC多位点模式与iCCA S486；HCC队列差异作为条件性候选；炎症/屏障方向先确认材料。各候选不是并行开展所有临床检验的授权。

一、回到原始问题
原始目标是CGNL1磷酸化与不同消化系统疾病的关系，S283是后续候选。早期规格为所有符合测量条件的CGNL1位点，而非仅S283。当前实际资料主要是胰腺、肝胆、胃、食管和结直肠肿瘤；不能称覆盖所有消化系统疾病。
本轮回顾8个具有配对位点结果的队列/来源，共81条既有位点或区域记录，另保留26个登记来源的资格快照。81条包含低覆盖和不可评估，不是81个阳性检验；8个来源也不等于8个已证明独立的验证队列。
未重新拟合模型或计算新P值。原始结果、原始P、原FDR、效应、覆盖和敏感性完整保留。不同预处理的效应不直接合并；原文用缺失/小样本处理得到空P时不补造P。当前默认原始P<0.05识别名义信号，但本文保留历史冻结规格，不偷偷重算。

二、候选C01：PDAC的CGNL1多位点模式
问题：CGNL1是否存在疾病相关、不能由总蛋白统一变化解释的位点组合？
支持：P01 S198降低（n35，Δ=-0.263659，P=0.00173010，FDR=0.00692038）；S297相对总蛋白升高（n35，Δ=+0.606412，P=0.000283359，FDR=0.00212288）；S708相对总蛋白升高（n21，Δ=+0.628323，P=0.000353813，FDR=0.00212288）。S283也有信号（n44，Δ=+0.276201，P=0.00758637，FDR=0.0227591），但并非唯一或应天然优先的位点。
可行性：S198/S297共同患者28对，S283/S708共同患者21对，可预先固定两组同患者位点对比。三位点S198/S297/S708共同患者仅14对，不适宜直接建立复杂评分。
下一步：先做共同患者的位点间变化对比和分子/分母分解。共同总蛋白和同组学全局中心化项可在差值中抵消；定位、肽段响应和缺失偏差仍需核对。此新对比尚未执行。
反证与边界：PDAC目前仍是一个独立发现队列。CDAP重处理支持S283/S708但不能算独立验证。S296-S299邻近位点的定位必须严格区分。两套胃癌和食管结果不支持直接升级为泛消化道统一标志。探索价值在于可检验的位点功能状态，而不是多列几个显著P值。
判断：优先候选；计算可行，最终仍需独立PDAC配对材料及位点功能验证。

三、候选C02：肝内胆管癌S486相对总蛋白降低
问题：iCCA中S486是否存在超出总CGNL1变化的稳定改变？
支持：B03有60对完整患者，42降18升，Δ=-0.319315，P=0.000120399，FDR=0.000240798；不额外中心化时Δ=-0.396850，P=0.00373336，方向和原始显著性均保留。它比B03 S283对现有中心化敏感性更稳定。
可行性：B03可直接进行位点、总蛋白、缺失及个体影响核查；S283/S486在同60例中都有值。B02有143例肿瘤S486及总蛋白完整值，但只能作为后续疾病背景材料，不能验证肿瘤-正常变化。
下一步：先确认S486定位和源处理，排除总蛋白/细胞组成及少数个体解释；通过后再提出一个有生物学依据的具体表型，不把S283血管侵犯阴性后再试所有终点当成研究策略。
反证：P01 S486 n26、P=0.920455；食管E01 n25、P=0.140915，未支持跨病种普遍信号。H01 S486方向受中心化影响，不能作稳定确认。不同队列一显著一不显著不证明癌种特异性。
判断：优先候选；当前单队列疾病差异较清楚，独立iCCA配对确认仍缺。

四、候选C03：HCC位点模式的队列差异
问题：不同HCC队列的CGNL1位点变化，是测量处理不同，还是在可比测量下真实存在背景依赖？
支持：H02 S283 n38、Δ=-1.483808、P=1.89505e-5、FDR=3.79010e-5；S297 n21、Δ=-1.645554、P=0.000161171、FDR=0.000161171。H01 S297却上升（n34、Δ=+0.800823、P=3.74042e-7、FDR=7.48085e-7），且未额外中心化仍同向。H03 S283降低，但其GRR插补来源尚有限制；H03 S297现行主分析P=0.134890，不能写成显著同向复现。
可行性：S283/S297共同患者H01/H02/H03分别34/17/27对。可先开展固定口径的位点间对比，再判断是否值得引入预设疾病背景。
优先排除：H01 S283中心化翻转、H03插补与定位信息不足、不同入选条件，以及HCC患者重叠未排除。不能先把H02早期队列差异解释成分期机制。
判断：条件性候选。先做技术归因，只有测量可比后仍保留差异，才考虑疾病背景交互。不能用队列间P值差异替代交互检验。

五、候选C04：非肿瘤炎症/屏障损伤，作为储备
问题：炎症或应激过程中CGNL1磷酸化状态是否可逆变化，并区别于单纯蛋白缺失？
支持：已有CGNL1与连接、微管/极性的直接蛋白功能证据；2026 VEO-IBD会议摘要已报告CGNL1变异/丢失与肠道表型。
目前不足：本项目尚未确认适用的非肿瘤患者/模型配对位点-总蛋白数据；摘要不是独立审计的磷酸化数据，也不是S283功能证明。“CGNL1丢失导致屏障问题”的一般命题已有披露。
下一步：仅做材料资格检查，要求同样本位点与总蛋白、明确刺激/对照或恢复设计及生物学重复。获得合格材料后再选位点；不以RNA充当磷酸化读出。
判断：储备，不立即投入大规模关联分析。

六、共同的功能验证出口与创新性判断
2024原始研究支持CGNL1 250-420区域与微管结合，删除365-377影响结合。S283/S297位于该区域，给出可以检验的位置依据，但不证明这两个位点已经调控该功能；S486/S708不能套用同一区域解释。CGNL1-CAMSAP3/连接定位的研究提供微管组织和上皮极性读出。
如果某个疾病-位点模式经独立材料确认，下一步可检验“位点状态是否改变CGNL1定位/结合及细胞功能”。需要确证位点、匹配总蛋白表达和适当救援；磷酸模拟突变不能单独替代实际位点磷酸化证据。不要把CGN的AMPK机制、circCGNL1研究或motif预测直接移植给线性CGNL1。
创新性目前是待验证的空间，不是已获证明的新机制，也不以检索未命中断言从未有人研究。

七、暂不作为主线的方向
不强推S283泛消化道统一改变：胃癌63对P=0.0647711、30对P=0.807830，食管31对P=0.946101，已有关键阴性。保留完整效应/P/FDR，阴性不等于等效或没有任何疾病功能。
不继续将B02 S283-血管侵犯作为优先问题：四情景P=0.487–0.580，结果一致不支持该关联；它不能否定其他位点或其他疾病。
不先铺开泛预后、药敏、激酶网络或机器学习模型。必须由合格位点、明确定义的疾病问题和可排除的竞争解释决定方法。

建议顺序：先整理完整疾病×位点证据和同患者覆盖（本轮已完成），再有限推进C01/C02；C03先做可比性诊断，C04先查材料。此处给出候选及判断依据，不启动新检验、不确定唯一主线。

主要文献：
https://pmc.ncbi.nlm.nih.gov/articles/PMC11513637/
https://pmc.ncbi.nlm.nih.gov/articles/PMC10184829/
https://academic.oup.com/ibdjournal/article/32/Supplement_1/izag006.073/8436212
全部数值以服务器既有结果为准；ALL_SITE_EVIDENCE_LEDGER.json/tsv保留81条结果、来源、效应、P、FDR及边界。
'''
(O/'CGNL1_DIRECTION_REVIEW_20261023.txt').write_text(summary,encoding='utf-8')
for k,v in baselines.items():assert sha(Path(k))==v
save('REVIEW_ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),review_scope='8 paired sources;81 stored site/region records;26-source registry;targeted primary literature',
 source_result_files_unchanged=len(baselines),source_hashes=baselines,new_inferential_tests=0,new_models=0,
 candidates=4,unique_mainline_selected=False,whole_project_complete=False,
 report_sha256=sha(O/'CGNL1_DIRECTION_REVIEW_20261023.txt'),artifacts=[dict(path=str(f),sha256=sha(f)) for f in O.iterdir() if f.is_file() and f.name!='REVIEW_ACCEPTANCE.json']))
print(json.dumps(dict(report=str(O/'CGNL1_DIRECTION_REVIEW_20261023.txt'),sha256=sha(O/'CGNL1_DIRECTION_REVIEW_20261023.txt'),rows=len(records),common_patients=common),ensure_ascii=False))
