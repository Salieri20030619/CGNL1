from pathlib import Path
import json,datetime,hashlib,collections
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
O=P/'continuous_exploration_20261023_v1';W=O/'wave3_source_purity';F=O/'functional_sources'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert read(p)==x
a=read(W/'ACCEPTANCE.json');assert a['R_QA_pass']
for r in a['artifacts']:assert sha(Path(r['path']))==r['sha256']
alignment=read(F/'HUMAN_MOUSE_ALIGNMENT.json');maps=[x['mapping'] for x in alignment['records']]
assert all(m==maps[0] for m in maps)
assert maps[0]['199']['human_position']==198 and maps[0]['205']['human_position']==204
ns={};exec((P/'audit/CGNL1_FOLLOWUP_20261023_v1/followup.py').read_text().split('baseline=read')[0],ns)
source=dict(ns['rows'](F/'Mer_Source_Data.xlsx','Fig. 1',{3,4,9}));assert source[9]['A9']=='Cgnl1 S199/T205'
values={source[4][c+'4']:source[9][c+'9'] for c in ['B','C','D','E']}
evidence=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),doi='10.1038/s42003-023-05287-y',
 source_workbook=str(F/'Mer_Source_Data.xlsx'),source_sha256=sha(F/'Mer_Source_Data.xlsx'),sheet='Fig. 1',row=9,original_feature='Cgnl1 S199/T205',
 source_unit=source[3]['A3'],values=values,source_unit_log_base_not_verified=True,
 descriptive_WT_Gas6_difference=values['WT+G6']-values['WT'],descriptive_KO_Gas6_difference=values['MerKO+G6']-values['MerKO'],
 mapping=maps[0],alignments_checked=len(maps),gap_settings=sorted(set(x['gap_open'] for x in alignment['records'])),
 supports='Published receptor/ligand perturbation changes a mouse Cgnl1 combined phosphosite-region signal; supports regulatability of this region.',
 does_not_establish=['isolated humanS198 response','humanS297 orS486 regulation','Mer directly phosphorylates these serine/threonine residues','CGNL1 phosphosite necessity/sufficiency','digestive disease mechanism','biological replicate count from three technical injections'],
 total_proteome_boundary='Methods describe total-proteome samples from full-medium confluent cultures versus serum-starved +/-15minGas6 phospho samples;do not assume matched same-condition totalprotein correction.',
 statistical_policy='Four source-summary values only;no newP,confidenceinterval or foldchange conversion without verified scale and biological replicate data.',
 novelty='Generic upstream regulation of this CGNL1 region already has published evidence. Any newGI-specific,site-specific functional claim requires additional evidence.')
save(F/'MER_CGNL1_EVIDENCE_QUALIFICATION.json',evidence)
r=read(W/'RESULTS.json');audit=read(W/'SITE_MEASUREMENT_AUDIT.json')
text='''CGNL1持续探究第三轮：来源定位与纯度竞争解释（2026-10-23）

Goal保持active。本轮完成两个预先限定的纯度模型、候选位点来源审核，以及新找到的Mer/Gas6功能文献源数据资格检查。没有确定唯一主线。

一、PDAC S198/S297的来源证据需要单独补强
S198、S297等五个候选位点的15aa窗口与人CGNL1参考序列的坐标一致。这个检查证明坐标对应，不能替代谱图定位概率。
已下载的另一套PDAC CDAP磷酸位点导出只含S23、S112、S283、S486、S514、S708共6个CGNL1位点，没有S198/S297。其S283/S708替代处理核验不能用于证明S198/S297也已跨处理复现。这是处理/覆盖差异，不据此宣布harmonized的S198/S297错误。
C01仍保留：同患者S198-S297发现效应较强，但需原始定位与处理来源证据，不能只凭P值排在所有候选之前。
B03的S486/S283两条源肽段及其坐标均能直接核对，来源可追溯性较好；导出表仍不提供逐观测定位概率。H01的S297映射到富含丝氨酸的肽段，也应保留这一定位证据边界。

二、纯度代理检查
28例C01_A患者均与现有表型文件精确匹配，WES_purity和WGS_purity均无缺失。只检验这两个预先指定字段；没有展开所有细胞评分或通路寻找显著性。
模型为位点差值~DNA推算肿瘤纯度，HC3稳健标准误，效应单位为纯度增加10个百分点：
WES：斜率+0.004287，95%CI[-0.229700,+0.238275]，P=0.970245，FDR=0.970245，描述性Spearman rho=-0.043574。
WGS：斜率+0.046131，95%CI[-0.239403,+0.331666]，P=0.742481，FDR=0.970245，描述性Spearman rho=+0.009582。
Python与R独立HC3计算一致。当前两个代理变量没有解释位点模式，但样本量小、区间较宽，正常组织组成未纳入，且未证明纯度与磷酸化来自同一物理组织块。不能由这些阴性结果宣称排除了细胞组成或证明了细胞自主机制。

三、新功能线索：Mer/Gas6与小鼠Cgnl1位点区域
来源为Regulation of brain endothelial cell physiology by the TAM receptor tyrosine kinase Mer，Communications Biology 2023，DOI10.1038/s42003-023-05287-y。
原文使用WT/Mertk缺失小鼠脑内皮细胞，并给予或不给予Gas6刺激15分钟。源数据工作簿Fig.1第9行为合并特征“Cgnl1 S199/T205”，四个条件各一个汇总值：WT13.8767088，WT+G613.256?（以以下精确源值为准）。
'''
# Write exact values programmatically to avoid transcription errors.
text=text[:text.index('，四个条件各一个汇总值：')]+'，四个条件各一个汇总值：\n'
for k,v in values.items():text+=f'{k}: {v:.12f}\n'
text+='''源表单位为“log converted total intensity”，对数底未核实，未将数值换算成倍数，也未补造P值。
人鼠全长BLOSUM62比对在三种gap-open设置、24条受检最优比对中一致：小鼠S199对应人S198；小鼠T205对应人S204（氨基酸由T变为S），不是按编号直接平移。小鼠S284/S298分别对应人S283/S297。
因此，原文支持“这一位点区域可以响应上游受体/配体扰动”，但不能拆解成孤立S199或人S198效应。它没有证明Mer直接磷酸化CGNL1，没有CGNL1位点突变救援，也不是消化系统疾病验证。文中三次进样不能自动当作三份生物学重复。其总蛋白培养条件与急性磷酸化实验不同，不能直接混用为同条件分母。
这条文献增加了可检验的功能入口，同时表明“该区域存在上游调控”不是全新发现。创新仍需要疾病特异、位点特异的功能证据。

四、候选比较的更新
C01：发现信号保留，纯度代理未解释，但单个位点定位/处理独立性不足，仍为条件性候选。
C02：源肽段可核查，S486相对S283模式在B03和H01均有支持；不能据此宣称iCCA专属。
C03：HCC相对位点模式仍保留，H02阴性、H03插补和跨队列独立性问题未消失。
C04：上一轮相关材料访问缺口仍在，未把其标为生物学阴性。本轮没有重复请求同一失败附件。

接续工作：优先补齐C01位点定位/处理证据，核对新功能研究是否提供可独立复核的观测层数据，并在新的合法材料入口出现时继续C04资格检查。新的结果不得替代前两轮关键阴性。
'''
f=O/'CGNL1_CONTINUOUS_PROGRESS_20261023_v2.txt';f.write_text(text,encoding='utf-8');assert '13.256?' not in f.read_text()
status=read(O/'STATUS.json');status.update(utc=evidence['utc'],wave3='source_and_purity_checks_complete_QA_pass',
 C01='discovery signal retained;WES/WGS purity association null;S198/S297 absent from alternativeCDAP export;localization/processing evidence still needed',
 functional_lead='Mer/Gas6 mouseCgnl1 combinedS199/T205 region;humanS198/S204 mapping;notGI or isolatedsite validation',next_task='C01 source-level localization and functional observation availability;C04 new-source access')
save(O/'STATUS.json',status)
save(W/'DELIVERY_ACCEPTANCE.json',dict(utc=evidence['utc'],report_sha256=sha(f),model_QA_pass=True,functional_source_sha256=sha(F/'Mer_Source_Data.xlsx'),
 alignment_sha256=sha(F/'HUMAN_MOUSE_ALIGNMENT.json'),new_primary_tests=2,new_functional_tests=0,goal_status='active',whole_project_complete=False))
print(json.dumps(dict(report=str(f),sha256=sha(f),mapping=maps[0],values=values),ensure_ascii=False))
