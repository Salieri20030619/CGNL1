from pathlib import Path
import json,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'H02_localization_reaudit_20261024/evidence_link'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==z
r=read(D/'RESULTS.json');q=read(D/'INDEPENDENT_PARSER_QA.json');assert q['independent_parser_max_error']==0 and q['additional_direct_pass_records']==0 and r['all_four_supported_n']==1
text='''CGNL1 H02底层鉴定与原患者测量连接结果（2026-10-24）
主要结论：底层evidence表补强了S283定位依据，但不足以把原H02的S283/S297对比升级为逐样本均直接定位的独立确认。

核查规则与范围
在查看新连接结果前冻结：仅核对原C03_H02的17对患者、S283和S297共68个组织-位点测量。按Experiment、位点表Evidence IDs和evidence表Phospho(STY)site IDs精确连接；概率序列须与未修饰肽段一致，肽段须唯一映射到CGNL1参考序列。
主核查要求Type为MSMS或MULTI-MSMS、MS/MS Count>0、MS/MS IDs非空，且目标残基定位概率≥0.75。MULTI-MATCH不作为独立MS/MS定位证据，MULTI-SECPEP单列。本轮不计算疾病效应P值、不改变原病例纳入。

结果
双向链接得到248条相关证据。127条存在可解析概率序列，121条概率字符串为空；后者全部为MULTI-MATCH。独立逐字符解析与主解析的目标概率完全一致，未发现非空字符串解析失败。
S283：34/34个组织测量至少有一条满足主核查标准的支持。
S297：6/34个组织测量至少有一条满足主核查标准的支持。
仅L012这1对患者的四个测量全部满足该标准。因此当前不能构造一个有有效生物学重复、四测量均直接定位的H02确认性对比；也没有对这1对患者补造P值。
另外73条记录在evidence侧列出S297候选site ID，但不在S297位点表的Evidence IDs中。这些不是被自动判为源数据错误；我们单独核查了其概率。非空目标概率最高为0.486，没有增加任何符合≥0.75的直接支持记录。因此在本次材料范围内，放宽为单向候选链接不改变6/34的结论。

证据边界
缺乏当前标准下的直接支持，不等于该磷酸化状态不存在，也不等于匹配传递强度全部错误。原H02采用全局定位≥0.75和正强度，其17例原P=0.108887的阴性结果继续保留。此前逐实验Localization prob字段虽为空，但现在可以通过底层记录补充一部分实际定位信息。
支持记录中未发现Reverse或Potential contaminant标记。其鉴定PEP范围为1.6892e-26至0.051252；本次没有临时增加PEP阈值，不能将“定位概率≥0.75”写成“每条鉴定PEP<0.01”或完整鉴定FDR验收。该核查也没有消除S297多磷酸化形式混合问题。

对方向选择的影响
C03应继续作为分来源相对模式与异质性问题，不能称H02已经独立确认S283/S297机制。H01与H03也应分别保留其定位、插补和患者独立性边界。C02的S486因源肽段较明确仍适合优先补定位与独立重复，但尚未完成机制验证。
C01的同序列形式比较阴性继续保留；C04单文件读取已通过，下一步推进预先冻结的全鼠蛋白组搜库与定位资格。Goal保持active。
来源证据、68个测量连接、248条相关证据、73条单向候选核查、独立解析QA均保存在H02_localization_reaudit_20261024/evidence_link。原数据和统计未覆盖。
'''
f=O/'CGNL1_H02_EVIDENCE_LINK_20261024.txt';f.write_text(text,encoding='utf-8');assert f.read_text()==text
save(D/'ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),parser_QA_pass=True,source_link_sensitivity_no_additional_pass=True,original_H02_null_retained=True,new_effect_tests=0,report_sha256=sha(f),artifacts=[dict(path=str(D/n),sha256=sha(D/n)) for n in ['SPEC_BEFORE_COUNTS.json','RESULTS.json','MEASUREMENT_LINKS.json','EVIDENCE_DETAILS.json','LINK_DISAGREEMENTS.json','INDEPENDENT_PARSER_QA.json','ONEWAY_CANDIDATE_CHECK.json']]))
s=read(O/'STATUS.json');s.update(H02_evidence_link='S28334/34,S2976/34 directMSMS+localprob>=.75;all4only1patient;original17patientnullunchanged',latest_H02_evidence_report=str(f),next_task='FreezeC04wholemouseproteomesearch/localization and execute singlefilepilot');save(O/'STATUS.json',s)
print(json.dumps(dict(report=str(f),sha256=sha(f))))
