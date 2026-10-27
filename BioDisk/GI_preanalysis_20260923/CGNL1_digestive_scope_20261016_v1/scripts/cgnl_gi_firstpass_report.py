from pathlib import Path
from datetime import datetime,timezone
import json,hashlib
O=Path.cwd();now=datetime.now(timezone.utc).isoformat()
records=[
 {'id':'GI01','disease':'gastric cancer','question':'peritoneal dissemination and epithelial state transition','molecule':'linear CGNL1 transcript','source':'https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.939969/full','location':'Figure 2E; cohort results','evidence':'31 gastric cancer patients;12 P+/Cy+,19 P0/Cy0;CGNL1 reported fold change4.82,P0.0081;figure visually verified;FDR not shown in this figure','interpretation':'clinical association, not CGNL1 perturbation or established epithelial source','priority':'first direct-evidence lead, not selected final mainline'},
 {'id':'GI02','disease':'gastric cancer','question':'peritoneal dissemination program','molecule':'linear gene-expression signature member','source':'https://www.nature.com/articles/srep22371','location':'Figure1a gene list; Results','evidence':'CGNL1 appears in published graphical gene list; experimental focus is DDR2','interpretation':'supporting screening lead, do not transfer DDR2 functional validation to CGNL1'},
 {'id':'GI03','disease':'gastric cancer','question':'ascites survival, plasticity and therapy response','molecule':'CGNL1 relevance untested in this source','source':'https://www.nature.com/articles/s41467-023-36310-9','evidence':'single-cell ascites resource,35 participants acrosswith/withoutGCPM;not35 qualified malignant paired donors','interpretation':'candidate material; malignant coverage, pairing and independence must be qualified; published paligenosis result is not novel CGNL1 evidence'},
 {'id':'GI04','disease':'stomach and gastroesophageal junction','question':'metaplasia to malignant progression','molecule':'CGNL1 relevance untested','source':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288235','evidence':'public spatial resource for GEJ/metaplasia/cancer; EAC biology must be separated from ESCC','interpretation':'problem-led data opportunity, not a target-positive result'},
 {'id':'GI05','disease':'esophageal squamous cancer, cross-disease reference HNSCC','question':'junctional regulation during malignant progression','molecule':'linear CGNL1','source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC8227492/','evidence':'miR31-3p direct CGNL1 regulation demonstrated in HNSCC;esophageal expression discussion is not same-system functional validation','interpretation':'cross-squamous rationale only;ESCC-specific target qualification required'},
 {'id':'GI06','disease':'pancreatic cancer','question':'apoptosis','molecule':'hsa_circ_0035435 / circCGNL1','source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC10829403/','evidence':'circRNA functional study; distinct from linear CGNL1 protein','interpretation':'separate molecular project if chosen;cannotborrowthismechanismforlinearCGNL1'},
 {'id':'GI07','disease':'HBV-related liver cancer','question':'growth and AKT signaling','molecule':'hsa_circ_0035431 / circRNA_10156','source':'https://medsci.org/v17p1351.htm','evidence':'different circular isoform from pancreatic hsa_circ0035435','interpretation':'not evidence of linear CGNL1 function;not same circular molecule across diseases'},
 {'id':'GI08','disease':'pancreatic and other GI tumors','question':'FGFR2 alteration biology','molecule':'CGNL1 as a fusion partner','source':'https://aacrjournals.org/mct/article/24/10_Supplement/B048/766223','evidence':'conference abstract listsCGNL1amongFGFR2fusionpartners;noindividualdiseaseassignmentinexcerpt','interpretation':'fusion biology separate from wild-type expression;not prior evidenceforordinaryCGNL1role'},
 {'id':'GI09','disease':'gastric cancer','question':'prognostic model','molecule':'linear CGNL1 as multigene signature component','source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC8806649/','evidence':'previous multigene model includesCGNL1','interpretation':'generic CGNL1/prognosis report would have limited novelty;not proof of specific mechanism'},
 {'id':'GI10','disease':'small bowel/appendiceal/anal/gallbladder','question':'not yet prioritized','molecule':'linear CGNL1','source':'focused search log','evidence':'no convincing direct disease-specific functional basis retrieved in first pass','interpretation':'search limitation, not biological absence;noinventedpositiveornegativeconclusion'}
]
(O/'metadata/FIRST_PASS_EVIDENCE_REGISTER.json').write_text(json.dumps({'utc':now,'records':records,'scope':'focused search, not systematic review; no cross-cancer expression or survival scan performed'},ensure_ascii=False,indent=2),encoding='utf-8')
text=f'''CGNL1 消化系统肿瘤范围重置：首轮检索记录
记录时间：{now}
状态：问题发现与文献初筛，不是已经通过验证的新主线。

判断
之前过早同时固定了CGNL1、结直肠癌肝转移和空间深化，容易把“如何解释局部关联”替代成研究主线。新的工作应先比较疾病问题的重要性、靶点的直接证据、创新性和可取得的材料，再选择分析方法。
广泛检索用于选择最终聚焦的问题，不是把所有癌种合成一篇泛癌表达/预后筛选，也不是寻找最小P值。旧空间结果保留为已验证边界，不继续为维持该候选扩展同一批数据。

范围
消化道：食管（鳞癌/腺癌分开）、胃、小肠、阑尾、结直肠、肛管。
相关消化系统：胰腺、肝细胞癌、胆管/胆囊癌，另列比较。原发肝癌不与结直肠癌肝转移混同。
暂以线性CGNL1/paracingulin为线索；环状RNA、融合基因另列，不能互相替代。

首轮得到的实质线索
1. 胃癌腹膜播散具有优先核查价值。
   2022年原始研究31患者（12例P+/Cy+、19例P0/Cy0），Figure2E列CGNL1 Fold Change=4.82、P=0.0081。本轮已直接查看原图确认；图未列对应FDR，只作为名义显著关联线索。不能据此确认肿瘤细胞来源或因果作用。
   2016年DDR2研究的播散表达图中也出现CGNL1，但直接功能实验研究的是DDR2，不可移植成CGNL1机制。已有报道同时意味着“发现CGNL1升高”本身不新。
   值得比较的问题包括播散过程中的上皮状态转换、腹水环境中存活与治疗残留、重新黏附/生长的关系。这里只提出问题域，未确认CGNL1参与这些过程。

2. 胃/胃食管结合部癌前演进可以作为另一类问题域。
   已有肠化生、早癌及结合部谱系/空间资料。可比较“适应性上皮重塑与恶性演进的分界”，但目前没有完成CGNL1在这些人群中的资格验证，不因数据多就指定为主线。已有肠化生和谱系研究丰富，创新需要具体新增认识。

3. 食管鳞癌应与食管腺癌分开。
   HNSCC原始研究提供miR-31-3p直接调控CGNL1的证据，可作跨鳞状上皮背景；不能当作食管癌功能已被证明。食管表达名单或跨癌数据库观察的证据等级更低。

4. 胰腺癌和肝癌文献必须先辨别分子。
   胰腺癌2024年功能论文研究hsa_circ_0035435（circCGNL1）；HBV相关肝癌研究的是hsa_circ_0035431（circRNA_10156）。两者既不是线性CGNL1蛋白，也不是同一个环状转录本。把这些论文并成一个CGNL1机制会造成错误。
   融合伙伴记录同样不能替代野生型CGNL1表达/功能的证据。

5. 小肠、阑尾、肛管、胆囊等少见病种暂不优先。
   本轮聚焦检索未找到足够直接的线性CGNL1研究起点。这不是系统综述，也不能写成CGNL1在这些病种无作用。

数据资格提醒
胃癌可用资源并不等于有足够腹膜转移患者。GSE163558官方记录为6患者、10样本，其中只有1份明确腹膜转移样本，不能用全体细胞数支撑腹膜播散推断。2023年胃癌腹水单细胞研究涵盖35名不同情境参与者，也不能直接当作35对恶性细胞样本。当前只核对来源入口，没有重算靶点效应。

下一阶段的决策框架
先在各病种分别填写：重要临床/生物学问题；CGNL1已有直接与间接证据；竞争解释；别人已经回答什么；真正独立且可取得的材料；一个可被否定的主要终点。
据此比较候选，再选择少数方向做靶点资格验证。癌种、过程、细胞来源和方法不同时过早锁死。腹膜播散当前是较具体的新线索，尚未被选为最终主线。

来源
[1] 胃癌腹膜受累队列及Figure2E：https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.939969/full
[2] 胃癌播散表达与DDR2研究：https://www.nature.com/articles/srep22371
[3] 胃癌腹水单细胞研究：https://www.nature.com/articles/s41467-023-36310-9
[4] GSE163558：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163558
[5] 结合部空间数据：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288235
[6] HNSCC直接调控研究：https://pmc.ncbi.nlm.nih.gov/articles/PMC8227492/
[7] 胰腺circCGNL1：https://pmc.ncbi.nlm.nih.gov/articles/PMC10829403/
[8] 肝癌circRNA_10156：https://medsci.org/v17p1351.htm
[9] 胃癌既有预后签名：https://pmc.ncbi.nlm.nih.gov/articles/PMC8806649/

完整范围、下载回执和分子分类登记保存在：{O}
本轮没有宣称某一新癌种已通过CGNL1机制或预后验证。
'''
(O/'report/CGNL1_消化系统范围重置与初筛.txt').write_text(text,encoding='utf-8')
(O/'STATUS.json').write_text(json.dumps({'utc':now,'stage':'broad_literature_first_pass_complete','new_mainline_selected':False,'priority_lead':'gastric peritoneal dissemination, qualification pending','scope_reset_recorded':True,'raw_or_target_statistics_generated':False,'report':'report/CGNL1_消化系统范围重置与初筛.txt'},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'report':str(O/'report/CGNL1_消化系统范围重置与初筛.txt'),'records':len(records)},ensure_ascii=False))
