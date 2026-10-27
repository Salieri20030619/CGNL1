from pathlib import Path
import json,hashlib,datetime,collections
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';H=O/'H02_localization_reaudit_20261024';C=O/'C04_raw_pilot'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==z
hr=read(H/'H02_target_evidence_receipt.json');dr=read(H/'H02_evidence_download_receipt.json');records=read(H/'H02_CGNL1_evidence_rows.json');assert hr['crc_verified'] and dr['crc_verified'] and sha(H/'H02_CGNL1_evidence_rows.json')==hr['sha256_target'] and len(records)==hr['target_rows']>0
q=read(C/'ACQUISITION_QA.json');assert read(C/'ACQUISITION_ACCEPTANCE.json')['acquisition_pilot_pass']
text=f'''CGNL1缺失材料获取与单文件读取验收（2026-10-24）
本轮取得了H02底层鉴定记录和C04炎症队列的有效RAW文件，原有两项访问缺口取得实质进展。尚未据此新增疾病效应或机制结论。

H02 evidence.txt已完整取得
从PXD008373官方Phospho_MaxQuant_results.zip读取ZIP目录，定位txt/evidence.txt，仅下载该成员的698235104个压缩字节，没有重下17879997918字节的全归档。中途一次区段截断后保留实际有效前缀，按断点续传；最终解压字节数3413116245及成员CRC32=3212690271均与目录一致。
压缩成员SHA256：{dr['sha256_compressed']}。
随后流式扫描全部{hr['total_evidence_rows']}条鉴定记录，提取{hr['target_rows']}条包含Q0VF96/CGNL1标识的记录，再次通过完整成员CRC。相关记录的鉴定Type分布：{json.dumps(hr['types'],ensure_ascii=False)}。
这提供了Modified sequence、Phospho(STY) Probabilities、Type、Experiment、MS/MS IDs等底层字段。下一步须与原17对患者及位点ID精确连接，区分直接鉴定、匹配传递、定位不足或无对应记录；目前不把这些记录数量当生物学样本数，也不把“已取得表”写成“全部定位通过”。

C04公开FTPS路线打通
普通匿名FTP返回需要TLS；Python FTPS未通过证书链校验。Windows原生curl在正常用户环境使用系统证书库、保留证书校验后成功访问目录并下载预先锁定的phospho_fraction_1.raw。未关闭TLS校验，未重复HTTP429接口。
RAW为221702430字节，与公开目录一致；本地与服务器SHA256均为{q['raw_sha256']}。
ThermoRawFileParser固定2.0.0-dev版本转换正常退出，生成{q['mzML_bytes']}字节mzML；日志保留3个空谱警告（74195、74532、75054）。此开发版本的使用已记录，不冒充正式稳定版验证。
共75653条谱：MS1=5495、MS2=46774、MS3=23384。MS2中HCD与CID各23387条，符合原文双MS2设计。所有256655个导出MS3前体引用均能解析；22760条MS3在8个预设TMT报告离子位置中的至少一处有正峰，8个位置均有覆盖。
这是文件读取、采集结构和报告离子可用性检查。未搜库、未做位点定位、未证实Cgnl1检出；也没有由技术扫描数构造生物学P值。MS3报告离子存在不能替代鉴定及干扰/QC检查。
QA写入后曾因JSON把数值键转为字符串触发类型断言，已修正记录格式并核对已生成结果和输出哈希，没有因此重转RAW或改动科学数值。

下一阶段
H02：先按患者/实验、证据ID、序列和位点概率核对68个原主分析测量，保留17例原阴性和多磷酸化形式限制。C04：先冻结全鼠蛋白组target-decoy搜库、位点定位和PSM-MS3定量规则，完成单文件搜库资格，再决定是否扩至全部磷酸化分馏。单个分馏未见Cgnl1不能视为整个队列阴性。
原C01/C02/C03结果不变。Goal保持active，整体项目未完成。
材料、代码与验收在服务器continuous_exploration_20261023_v1下H02_localization_reaudit_20261024、C04_sources、C04_raw_pilot目录。
'''
f=O/'CGNL1_MATERIAL_ACQUISITION_20261024.txt';f.write_text(text,encoding='utf-8');assert f.read_text()==text
save(O/'MATERIAL_ACQUISITION_ACCEPTANCE_20261024.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),H02_fullmemberCRC_pass=True,H02_target_rows=len(records),C04_acquisition_pilot_pass=True,new_disease_tests=0,report_sha256=sha(f),whole_project_complete=False))
s=read(O/'STATUS.json');s.update(H02_source_evidence=f'Full evidence member CRCpass;{len(records)}CGNL1records acquired;sample/site linkage pending',C04='RAW FTPSacquired;conversion/MS2-MS3/reportermz passed;fullproteomesearchqualification pending',latest_material_report=str(f),next_task='H02 exactevidence linkage;freezeC04 wholemouseproteomesearch andlocalization');save(O/'STATUS.json',s)
g=read(O/'GOAL_PLAN.json')
for t in g['tasks']:
 if t['id']=='T04':t['progress']='FTPSraw access solved;singlefileacquisitionQApass;CGNL1identificationandlocalizationpending'
save(O/'GOAL_PLAN.json',g)
spec=read(O/'C04_sources/PAK_RAW_PILOT_SPEC.json');spec.update(status='single_file_acquisition_QA_pass_search_pending',current_blockers=['Wholemouseproteome identification/localization andPSM-MS3 quantification notyetqualified']);save(O/'C04_sources/PAK_RAW_PILOT_SPEC.json',spec)
print(json.dumps(dict(report=str(f),sha256=sha(f),H02_records=len(records))))
