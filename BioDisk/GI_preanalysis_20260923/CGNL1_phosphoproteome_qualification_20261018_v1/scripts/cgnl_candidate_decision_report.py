from pathlib import Path
import json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==z
w1=read(O/'wave1/RESULTS.json');w2=read(O/'wave2/RESULTS.json');purity=read(O/'wave3_source_purity/RESULTS.json')
evidence=[]
for wave in ['wave1','wave2','wave3_source_purity']:
 a=read(O/wave/'ACCEPTANCE.json')
 assert a.get('QA_pass',a.get('R_QA_pass'))
 for r in a['artifacts']:assert sha(Path(r['path']))==r['sha256'];evidence.append(r)
for n in ['source_methods/HARMONIZED_METHODS_AUDIT.json','functional_sources/PKAKO_CGNL1_QUALIFICATION.json','functional_sources/MER_CGNL1_EVIDENCE_QUALIFICATION.json','C04_sources/PAK_BIOLOGICAL_SAMPLE_MAPPING.json','C04_sources/PAK_RAW_PILOT_SPEC.json','C04_sources/PAK_RAW_PILOT_RETRY.json','source_requalification/COAD_REQUALIFICATION.json']:
 f=O/n;assert f.exists();read(f);evidence.append(dict(path=str(f),sha256=sha(f)))
assert len(w1)==7 and len(w2)==12 and sum(x['P'] is not None for x in w2)==7
assert next(x for x in w1 if x['test']=='C03_H02')['P']>.05
text='''CGNL1磷酸化与消化系统疾病：候选比较和后续决策依据
版本：decision_v1；以服务器当前结果、源材料和验收记录为依据。

我们现在能回答什么
当前数据支持部分消化系统肿瘤中CGNL1不同位点的相对变化并不一致。最有价值的后续问题是：这些位点模式是否在独立样本中重复，以及是否对应细胞连接/屏障相关调控。当前不能回答哪一个位点导致疾病，也不能建立通用的“CGNL1磷酸化升高或降低”结论。
本次探索目的不是复现原论文全部结论。原论文主要提供数据和处理背景；CGNL1位点组合是本项目探索性问题。方案是在已知部分位点结果后提出的，冻结发生在新对比计算前，不能称为盲态验证。

四条候选如何选择
C01 胰腺癌S198−S297：保留为条件性共同主候选。28对患者的相对位点改变具有一致方向和留一稳健性；两种DNA纯度代理未解释该变化。优势是多位点问题明确，并有S198区域Mer/Gas6与S297区域PKA扰动的异器官功能先验。缺口是同疾病独立验证、精确导出阶段/插补来源和逐观测定位。不能用同队列CDAP再处理充当外部验证。
C02 肝胆疾病S486相对模式：保留为条件性共同主候选。iCCA有60对患者、源肽段与数据成分可复核；H01中S486−S283同向支持。更适合研究肝胆疾病中的共同或不同调控，而非直接主张iCCA特异。缺口是逐观测定位、兼容测量下的疾病差异检验，以及S486功能证据。
C03 HCC S283−S297：保留为异质性比较分支。H01/H03有统计支持，H02未通过预设检验；H03插补及三来源患者重叠未消除。不得合并宣称三次独立验证。PKA文献增加了位点分别调控的可检验性，不消除疾病数据的这些问题。
C04 非肿瘤炎症/屏障：保留为材料资格分支。TCT小鼠研究已完成8例身份、5对照/3炎症及TMT通道映射，但Cgnl1检出仍未知。不能把下载失败写成阴性。其科学价值是连接原始疾病范围，但目前还不能与前三条按效应量或P值比较。
当前不选择唯一主线：C01与C02各有不同优势，没有证据支持只因C01的P值更小就胜出。若新增资源优先补材料，先确认C01精确位点测量与C02定位可复核性；若能开展实验，先做共同的靶向测量可行性小试，再决定投入哪一条。

主要结果（直接读取已验收JSON；效应为各源数据锁定尺度，不能跨源比较绝对大小）
位点对终点为同一患者：(A肿瘤−A正常)−(B肿瘤−B正常)。负值表示A的肿瘤-正常变化相对B更低，不等于A绝对去磷酸化或占有率下降。C02_A另为原有蛋白调整终点的复核。
检验采用冻结的双侧精确符号秩方法；BH分别在C01两项、C02两项、C03三项内计算。样本数是独立患者，不是通道或肽段数。
'''
for r in w1:
 text+=f"{r['test']} / {r['cohort']} / {r['contrast']}：n={r['n']}，均值={r['mean']:.6f}，向下/向上={r['down']}/{r['up']}，P={r['P']:.12g}，FDR={r['FDR']:.12g}。\n"
text+='''
外部/跨疾病覆盖检查（12个预先锁定场景，7项可检验；其余保留不可评估或描述性）
下列结果不是同疾病独立复现的统一集合；复用同一来源不增加独立性。P值差异不能证明癌种特异性。
'''
for r in w2:
 value='不可评估' if r['mean'] is None else f"均值={r['mean']:.6f}"
 stat='未检验' if r['P'] is None else f"P={r['P']:.12g}，FDR={r['FDR']:.12g}"
 text+=f"{r['test']}：n={r['n']}，{value}，{stat}，{r['status']}。\n"
text+='''
必须保留的阴性和敏感性
PDAC S283−S708在完整21对中P=0.320457，不支持预设差异。共同14例在不同处理下出现信号，是子集/处理敏感性，不能替代完整样本结果。
H02 S283−S297的P=0.108887；均值bootstrap区间与符号秩检验针对不同统计量，不能挑选区间来改称检验阳性。
C02在H02/H03只有4例/3例，均值方向相反且仅描述性；保留这些结果，不能宣称跨所有肝癌来源一致。
PDAC的WES/WGS纯度代理模型P分别0.970245/0.742481；这不等于排除组织组成，尤其没有正常组织的同等组成测量。
原B02四种合理批次方案均阴性，不继续扩展临床终点搜寻。Fudan、Britton及COAD再资格检查未提供当前预设的精确双位点患者数据，不能拆开多位点肽段替代。

测量与功能证据
统一重处理的补充方法规定PTMProphet定位概率≥0.75；当前BCM位点ID末尾1表示主要异构体而非概率。当前矩阵保留NA，不足以单独证明每个非缺失值均直接观测。CDAP没有S198/S297也不能反推统一流程有误。
Mer/Gas6作者数据涉及小鼠脑内皮Cgnl1 S199/T205组合特征；全长比对对应人S198/S204区域，不能分解成孤立S198效应。
PKA双催化亚基敲除源表涉及3对小鼠肾上皮克隆：S298（人S297）log2变化−1.73，S284（人S283）+0.09，总蛋白+0.08。支持位点分别调控的候选解释，不是直接激酶作用或消化道疾病机制证明。两条文献均不支持宣称这些区域的上游调控是首次发现。

下一步所需材料与验收条件
G01/C01测量：取得精确矩阵导出说明、观测掩码和S198/S297定位证据；确认蛋白异构体/肽段/样本连接。验收看来源追溯和定位，不看能否获得更小P。
G02/独立重复：独立PDAC患者需同时测到S198和S297，预先锁定相同终点；不能把同队列重处理或其他癌种当替代。肝胆比较需要兼容测量、明确患者独立性与直接疾病对比，不能根据分别显著/不显著判定特异性。
G03/功能：在相关消化系统模型中先建立S198、S283、S297、S486的靶向检测与总CGNL1测量；测试急性上游扰动后，才进入位点突变/救援和细胞连接/屏障表型。PKA或Mer干预需区分间接调控，磷酸模拟突变不能单独证明真实磷酸化机制。当前未开展湿实验，也未假定已有可用位点抗体。
G04/炎症材料：需要可读取的处理矩阵或RAW。单文件试验已冻结，转换工具安装/版本通过；但RAW连续两轮HTTP429，未转换。先通过221702430字节预期大小、哈希/读取、MS2-MS3关联与报告离子检查，再决定全磷酸化分馏搜库。单分馏未检出不能作为队列阴性。总蛋白与磷酸化不能混用TMT通道。来源中0.75分离特异性不能替代定位概率。
G05/HCC：逐来源核对患者重叠、H03原始观测/插补掩码和定位，再考虑合并；当前只做分来源比较。

完成状态
已完成可用数据的三轮分析/QA、关键来源和功能先验资格核查、四候选比较以及具体缺口与验收条件。尚未完成C04原始文件试验及C01精确观测来源闭合，未完成任何疾病特异的位点功能验证。因此不把本报告或计算QA称为整个项目完成；Goal暂保持active，整体完成审计未通过。
数据、代码、冻结方案、源文件和验收凭证位于服务器continuous_exploration_20261023_v1；本报告为决策交付，不替换原结果。
'''
f=O/'CGNL1_CANDIDATE_DECISION_v1.txt';f.write_text(text,encoding='utf-8');assert f.read_text()==text
requirements=[
 dict(id='R01',requirement='frozen exploratory contrasts and appropriate tests',status='proven',evidence=['wave1/SPECIFICATION_BEFORE_RESULTS.json','wave2/SPEC_BEFORE_RESULTS.json','wave3_source_purity/SPECIFICATION_BEFORE_RESULTS.json']),
 dict(id='R02',requirement='independent numerical QA and preserved null/sensitivity',status='proven',evidence=['wave1/ACCEPTANCE.json','wave2/ACCEPTANCE.json','wave3_source_purity/ACCEPTANCE.json']),
 dict(id='R03',requirement='C04 qualification before inference',status='qualification_progress_proven_analysis_unachieved',evidence=['C04_sources/PAK_BIOLOGICAL_SAMPLE_MAPPING.json','C04_sources/PAK_RAW_PILOT_SPEC.json']),
 dict(id='R04',requirement='candidate comparison by reliability replication disease function',status='proven_with_explicit_gaps',evidence=[f.name,'CANDIDATE_COMPARISON_v1.json']),
 dict(id='R05',requirement='C01 exact observation provenance',status='unachieved',evidence=['source_methods/HARMONIZED_METHODS_AUDIT.json']),
 dict(id='R06',requirement='site-specific GI mechanism or independent biological validation',status='unachieved_explicit_experimental_gap',evidence=['functional_sources/PKAKO_CGNL1_QUALIFICATION.json','functional_sources/MER_CGNL1_EVIDENCE_QUALIFICATION.json'])]
save(O/'GOAL_REQUIREMENT_AUDIT.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),requirements=requirements,completion_proven=False,artifacts=evidence,report=str(f),report_sha256=sha(f),new_tests=0))
s=read(O/'STATUS.json');s['latest_decision_report']=str(f);s['completion_audit']='not_proven_C04pilot_andC01observation_provenance_open';save(O/'STATUS.json',s)
print(json.dumps(dict(report=str(f),sha256=sha(f),requirements=len(requirements),evidence_files=len(evidence),completion_proven=False)))
