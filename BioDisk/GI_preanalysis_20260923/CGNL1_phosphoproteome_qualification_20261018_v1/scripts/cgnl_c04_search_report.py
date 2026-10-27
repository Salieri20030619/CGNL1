from pathlib import Path
import json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';W=O/'C04_search_pilot';V=W/'run_v2_all_spectra'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==z
r=read(V/'IDENTIFICATION_SUMMARY.json');q=read(V/'INDEPENDENT_QA.json');links=read(V/'CGNL1_SCAN_LINKS.json');assert read(V/'SEARCH_EXIT.json')['exit']==0 and q['independent_R_max_error']<1e-12 and q['batch_sensitive_PSMs']==0
assert all(len(x['CID_candidates'])==1 and len(x['MS3_linked'])==1 and x['HCD']['refs'] for x in links)
text=f'''CGNL1炎症队列单分馏搜库试验结果（2026-10-24）
主要结论：在MSV000081198的phospho_fraction_1.raw中得到一条通过预设PSM门槛的Cgnl1磷酸肽候选，但磷酸化位点尚未定位，不能据此进行位点特异的炎症效应推断。

方案与来源
使用Comet2026.03rev0（fa08489），官方二进制SHA256核对通过。数据库为当前UniProt小鼠参考蛋白组及异构体63328条，加113条由GPM官方accession清单取得的污染蛋白。此前用字符计数得到63338，已在搜索前按FASTA标题行更正；未改变序列。2个旧污染蛋白accession P04188/P04745未解析，明确保留此限制。
采用全库内部拼接目标—诱饵竞争，4工作线程，HCD MS2、前体±20ppm、片段0.02Da、胰酶全酶切、2漏切；TMT标记N端及K、C烷基化为固定修饰，M氧化与STY磷酸化为可变修饰，总可变修饰最多3个。CID与MS3原谱保留。此为新引擎/新数据库探索性重处理，不是原论文2014鼠ORF数据库及SEQUEST的逐参数复现。
所有搜索/后处理规则均在结果评估前固定。首次小批次运行因重复遍历全库改为一次载入；仅改变批大小，旧运行日志和部分输出保留。完整运行11分27秒、约4.8GB内存、正常退出。与旧运行可比较的1424条完整PSM，肽段和搜索分数一致。

鉴定层面结果
载入23055条可搜索HCD谱，输出{r['total_top_PSMs']}条顶级PSM（唯一scan-charge及scan均为18969），其中目标11259、诱饵7710。按e-value分组处理相同分数，以(D+1)/T估计并单调化q值；全体目标PSM中2326条达到q≤0.01。
对磷酸化PSM另作同样的分组估计，并同时要求全局q≤0.01和磷酸化组q≤0.01，共770条目标PSM通过；其中没有污染蛋白PSM。通过者的同位素偏差校正后绝对质量误差中位数约5.389ppm。PSM级q不等于蛋白或位点级FDR。
独立R与Python的18969行q值最大差约5.0e-16，数值QA通过。这不能代替独立生物学验证。

Cgnl1候选与关键未通过项
两条候选均为肽段SPTAPSPQAYSETK，映射到Cgnl1的4个蛋白记录/异构体，不是4次独立发现。
HCD scan26694：全局q=0.00068918，磷酸化组q=0.00194932，通过本轮预设PSM门槛。
HCD scan26745：全局q=0.00879121，但磷酸化组q=0.01328502，未通过双门槛，完整保留。
两条候选给出不同残基位置的修饰建议，AScorePro位点评分均为0。因此不能因为PSM通过或总体AScorePro分数较高，就指定位点，不能直接映射为人S261等位点并纳入原主线。该Comet定位评分目前属探索性功能，本轮不将其换算为PTMProphet定位概率。

与互补谱和MS3的连接
按同一MS1引用、相同电荷、前体m/z相差≤10ppm且保留时间相差≤2秒查找唯一CID，不单靠扫描号相邻。scan26694唯一连接CID26695及MS3 26697；scan26745唯一连接CID26747及MS3 26749。
这说明后续可以核查互补碎片和报告离子，但本轮未比较5对照与3炎症动物的效应，也未将技术扫描当作生物学重复。单个分馏结果不能代表整个队列。

后续验收与方向判断
下一步先冻结CID全库搜索/互补定位核查及MS3定量质量规则；明确位点或明确无法定位后，再决定能否做位点特异分析。不得为获得阳性而降低定位标准，不直接扩展到临床或机制结论。
C04已由“无有效材料”推进到“原始数据可读取、候选磷酸肽有PSM支持、定位未解决”。这值得继续资格核查，但尚不能与C01/C02/C03的配对效应量直接比较。原三条方向及其关键阴性不变。Goal保持active。

证据：C04_search_pilot中保存数据库/工具/参数/冻结方案；run_v2_all_spectra中保存全部PSM、q值、两条Cgnl1候选、独立R复算和CID/MS3连接。未覆盖原分析。
'''
f=O/'CGNL1_C04_SEARCH_PILOT_20261024.txt';f.write_text(text,encoding='utf-8');assert f.read_text()==text
save(V/'SEARCH_ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),search_exit_verified=True,independent_R_QA_pass=True,batchchange_score_QA_pass=True,CGNL1_candidate_retained=True,site_localization_pass=False,biological_inference_not_performed=True,new_disease_tests=0,report_sha256=sha(f),artifacts=[dict(path=str(V/n),sha256=sha(V/n)) for n in ['SEARCH_EXIT.json','IDENTIFICATION_SUMMARY.json','INDEPENDENT_QA.json','CGNL1_CANDIDATE_PSMS.json','CGNL1_SCAN_LINKS.json','ALL_PSM_QVALUES.json']]))
s=read(O/'STATUS.json');s.update(C04='HCDfullproteomepilotcomplete;770phosphoPSMs passdualPSMq;Cgnl1onepasscandidate,site score0;CID/MS3linked;biology nottested',latest_C04_search_report=str(f),next_task='FreezeCIDcomplementarysearch/localization andMS3quantquality gates');save(O/'STATUS.json',s)
g=read(O/'GOAL_PLAN.json')
for t in g['tasks']:
 if t['id']=='T04':t['progress']='RAWacquisitionandHCDwholeproteomesearch completed;Cgnl1candidateidentified atPSMlevel,site-localization unresolved'
save(O/'GOAL_PLAN.json',g)
print(json.dumps(dict(report=str(f),sha256=sha(f))))
