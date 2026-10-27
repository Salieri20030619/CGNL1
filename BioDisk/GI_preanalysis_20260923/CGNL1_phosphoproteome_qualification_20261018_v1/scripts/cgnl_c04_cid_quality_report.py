from pathlib import Path
import json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';W=O/'C04_search_pilot';H=W/'run_v2_all_spectra';C=W/'run_v3_CID';N=O/'C04_raw_pilot/noise_MS3'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==z
q=read(C/'INDEPENDENT_QA.json');assert q['pass_QA'];r=read(C/'IDENTIFICATION_SUMMARY.json');assert read(C/'SEARCH_EXIT.json')['exit']==0
candidate=next(x for x in read(C/'CGNL1_CANDIDATE_PSMS.json') if x['scan']=='26695');assert candidate['q_global']<=.01 and candidate['q_phospho']<=.01
noise=read(N/'REPORTER_NOISE_RESULTS.json');isolation=read(H/'ISOTOPE_ENVELOPE_PROXY_RESULTS.json');quality=read(H/'SPECTRAL_QUALITY_RESULTS.json');scores=read(H/'LOCALIZATION_SCORE_FIELD_AUDIT.json')
text=f'''CGNL1炎症队列互补CID搜索与报告离子质量核查（2026-10-24）
主要结论：Cgnl1候选磷酸肽在同一前体的HCD和CID谱中均有鉴定支持，MS3报告离子可读；具体磷酸化位置仍不一致，不能进入位点特异的炎症效应分析。

CID搜索与独立复核
在新CID结果评估前冻结全鼠蛋白组搜索，数据库、目标—诱饵竞争、TMT/磷酸化修饰和4线程保持不变，仅改为CID、片段bin宽1.0005及offset0.4等低分辨率设置。独立输出目录run_v3_CID保留完整记录。运行10分59秒，正常退出。
输出22243条顶级PSM，13683目标、8560诱饵；2545目标PSM通过全局q≤0.01，710条磷酸化目标PSM同时通过全局与磷酸化分组q≤0.01。独立R与Python q值最大差约5.55e-16。PSM q不代表位点FDR或生物学差异显著性。

Cgnl1候选的保留与未通过项
CID scan26695对应此前通过HCD门槛的同一前体，均鉴定肽段SPTAPSPQAYSETK。CID全局q=0.00082305、磷酸化组q=0.00240964；HCD对应值为0.00068918和0.00194932。
HCD倾向肽段第6位修饰，位点评分0；CID倾向第1位，位点评分5.24。尽管CID匹配碎片更多，不能据此忽略两者的位置不一致。总体AScorePro分数与位点评分不是同一指标，不能用较高总体分数替代定位。
第二个电荷态的CID scan26747全局q=0.01724745、磷酸化组q=0.03766105，未通过；另外两个低分候选也完整保留。不能把四个CID候选都称为可靠检出，更不能把同一前体的HCD/CID当两个独立生物学重复。
全局定位字段检查显示，通过HCD门槛的770条磷酸肽有772个位点评分，其中484个为正、386个≥13。因而候选的零分不是字段普遍失效；这里仅描述分数分布，没有用13替代预先验证的位点概率门槛。

MS3报告离子与噪声诊断
两个关联MS3（26697、26749）的8个预设TMT报告离子均为正且各有唯一匹配峰，最大质量偏差分别约2.72与3.18ppm。另存含厂商噪声数组的MS3导出，确认噪声m/z与强度峰逐项一致。
8通道合计I/noise分别为621.747与585.645；(I−baseline)/noise分别为621.954与585.859。它们是基于导出噪声的诊断量，尚未证明与原作者的SEQUEST summed S/N算法完全相同，不能直接声称复现了原S/N>386门槛。没有计算组间效应。

共隔离诊断及必要修正
第一前体的实际隔离中心为668.33905，报告单同位素m/z为668.0051、电荷3，差值对应约一个同位素间距。只统计单同位素峰会得到不适用的零占比；原结果和修正原因均保留。按已知电荷纳入同位素包络后，在真实隔离窗口内的强度占比代理为0.78684；第二前体为0.72489。
这是基于相邻MS1质心峰的共隔离诊断，不是原作者隔离特异性算法的精确复现，也不覆盖SPS干扰。没有根据这些代理临时改变阈值或选择患者。

下一步决策
单文件读取、全库鉴定、独立数值QA及谱间连接已完成；材料具备扩大检查的技术基础。下一步先固定完整磷酸化分馏及pY文件清单，再按同一流程检查其他分馏能否提供更清晰的Cgnl1位点，而不选择性挑选有利文件。
当前不做S257/S262等命名位点的炎症对比，不将鼠位点直接平移成人位点。即便后续定位改善，仍需总蛋白、通道归一化、同位素杂质/干扰及5对照3炎症的生物学样本资格核查，才能界定可回答的疾病问题。
原C01/C02/C03结果及关键阴性不变。Goal保持active，整体机制与主线验证未完成。
证据保存在C04_search_pilot/run_v2_all_spectra、run_v3_CID及C04_raw_pilot/noise_MS3。
'''
f=O/'CGNL1_C04_CID_QUALITY_20261024.txt';f.write_text(text,encoding='utf-8');assert f.read_text()==text
save(C/'ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),CIDsearch_complete=True,independent_R_QA_pass=True,paired_HCD_CID_peptide_agreement=True,site_localization_pass=False,reporter_peaks_available=True,quantification_biological_gate_pass=False,full_fraction_material_qualification_technically_feasible=True,report_sha256=sha(f),new_disease_tests=0,artifacts=[dict(path=str(x),sha256=sha(x)) for x in [C/'SPEC_BEFORE_RESULTS.json',C/'IDENTIFICATION_SUMMARY.json',C/'INDEPENDENT_QA.json',C/'CGNL1_CANDIDATE_PSMS.json',H/'SPECTRAL_QUALITY_RESULTS.json',N/'REPORTER_NOISE_RESULTS.json',H/'ISOTOPE_ENVELOPE_PROXY_RESULTS.json']]))
s=read(O/'STATUS.json');s.update(C04='PairedHCD/CID identifyCgnl1phosphopeptide;sitesdiscordant;8reporterswithnoise recovered;isolationproxy checked;no biologicalinference',latest_C04_quality_report=str(f),next_task='Freezecompletephosphofraction+pYmateriallist and expandsourcequalification,notselectivepositivefilepicking');save(O/'STATUS.json',s)
print(json.dumps(dict(report=str(f),sha256=sha(f))))
