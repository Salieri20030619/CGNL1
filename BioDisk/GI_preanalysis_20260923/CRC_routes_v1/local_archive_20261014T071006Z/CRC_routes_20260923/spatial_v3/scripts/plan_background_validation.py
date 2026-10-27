from pathlib import Path
import json
O=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3');B=O/'background_validation';B.mkdir(exist_ok=True)
plan='''# CGNL1空间初筛后的后台验证计划

日期：2026-09-24。状态：规划完成，后续验证计算尚未启动。当前初筛已完成并独立复核。所有后续结果写入background_validation的新版本目录，不覆盖初筛结果。

## 验证目标

区分三种解释：①肝/基质细胞混合使肿瘤标注区域CGNL1升高；②癌细胞内CGNL1状态与局部肝环境相联系；③信号主要来自宿主/基质细胞，需要改为研究相应细胞来源。验证依据是来源、患者一致性和外部证据，而不是要求“加入肝评分后必须显著”。

当前依据：11名患者9名近肝区域升高，主效应+0.541、P=0.0169；肝小叶和转移基质均11/11高于肿瘤区；spot等权及评分调整结果未达到显著。优先解决来源与混合spot问题，不继续堆叠同类bulk相关性。

## 任务顺序与交付

| 编号 | 任务 | 输入与执行规则 | 完成标志 |
|---|---|---|---|
| V01 / 优先 | 去除混合边缘后的空间稳健性 | 复用已下载24张切片和spot_target_coordinates.tsv.gz。主分析仅非复发样本。要求中心与至少5/6个100微米邻居均为作者肿瘤区，且6邻居均有注释；近肝组改为距肝小叶中心>100且≤300微米，远组≥600微米，每组≥20合格spot。严格6/6邻居作为敏感性。 | 逐片纳入/排除理由、有效患者数、患者级效应、原始及组成调整结果、空间掩膜图。结果定位为同一队列内部稳健性，不称独立验证。 |
| V02 | 独立单细胞来源参考 | 先审核GSE231559、GSE285991及正常肝参考。区分肝细胞、恶性上皮、内皮、成纤维/周细胞和免疫细胞。CGNL1及另外3候选不得参与细胞身份、恶性判定或解卷积标记选择。按患者生成表达和检出比例。 | 患者/样本/治疗/部位/物种映射、参考质量表、各细胞类型CGNL1来源图。 |
| V03 | 分離组成与细胞状态 | 在V02参考合格后，以共有基因做空间解卷积及组成预期表达估计。用独立参考中的CGNL1表达预测各spot混合表达，再评估残余空间梯度；与细胞类型特异模型交叉核对。肝评分不直接视为细胞比例。 | 每患者组成估计、预测/实测差异、肿瘤相关残余梯度、参考留一患者敏感性。模型不稳定或缺少参考类型时明确标记不可解释。 |
| V04 | 外部空间复现 | 优先Niche-DE公开患者4/5原始Visium；核对临床来源、治疗及患者重叠后再纳入。GSE217414仅1名明确乙状结肠患者，其余3名按CRC补充。Nature2024人空间为1名患者2切片，仅作探索性。 | 同口径来源、边缘排除与近/远比较；病例重复审查；结肠与CRC证据分别汇总。不能把不同患者的原发与转移拼成配对。 |
| V05 / 条件执行 | 公开干预数据中的响应 | 只有来源和空间关联仍值得推进时，再审查Nature2024配套bulk/scRNA/multiome干预数据及其他可用肝信号实验。核实CGNL1/Cgnl1同源映射、物种、处理、时间和独立重复。已知Plexin B2轴只作为现成干预条件，不预定为新的机制结论。 | 有重复的目标响应和候选上游排序；单次或无生物学重复只报告方向，不做因果外推。 |

V01可最先独立执行；V02与V04的元数据审核、公开输入获取可以同时准备；V03依赖V02；V05依赖前述结果判定。这里的同时执行指后台计算任务，不涉及调用子代理。

## 数据资格审查中的已知重点

- GSE231559包含原发、肝转移、邻近正常等组织。L8T1/L8T2属于同一患者标识的重复样本，必须核对并合并处理；C1T与L1T的名称不能证明同一患者配对。
- GSE285991是含人和小鼠、多种实验类型的SuperSeries。先定位人P01_LM至P10_LM对应的子系列和临床资料，禁止把小鼠Kupffer细胞实验或Stereo-seq混入人恶性上皮参考。治疗分层必须保留。
- GSE294385的11名空间患者尚无逐例结肠/直肠明确映射，本轮只作为CRC支持证据；不能将其改写成11例纯结肠癌。
- Niche-DE的5例来自不同来源，先检查与其他引用队列重叠，不能直接与现有11例相加。
- GSE335552及BJC2026 ISS已确认未覆盖CGNL1，排除直接目标验证；不通过插补将未测量变成“测得”。
- 作者受限Visium HD数据不作为首批依赖，暂不发出数据申请或联系作者。

## 统计和判定规则

以独立患者为统计单位。同患者切片、复发取样、ROI和spot不增加n。小于10名有效患者仅报告方向、效应量和一致性；n≥10时报告效应、区间、P及适用的家族FDR。每阶段在结果查看前保存具体终点、纳入规则、候选家族和分析版本。

V01以CGNL1的5/6邻居纯度规则为预设主要稳健性检查；6/6规则为敏感性。若同时报告其他3候选，在固定3候选家族内校正。多来源比较、相关性、外部验证和敏感性各自定义同类家族，不能混池，也不能筛出阳性后缩小校正家族。

减少混合边缘后病例数下降本身不是阴性；记录效应和置信范围。如果参考模型解释了大部分区域差异且无可靠癌细胞内残余，不把它强行写成癌细胞调控；可以转为宿主/间质来源定位。若严格区域、独立恶性细胞参考和外部数据支持同一联系，再优先设计针对该联系的湿实验。

空间邻近、解卷积残余、配体受体推断和横断面相关不等于因果。公开刺激/扰动后的响应可增强机制线索，但不替代CGNL1本身的功能扰动与救援。发表目标约10分仍取决于最终机制证据。

## 后台运行方式与完成验收

- 服务器：sgs10；根目录：/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/spatial_v3/background_validation。
- 每项任务使用独立run目录，保存输入清单/SHA256、参数JSON、软件版本、标准输出/错误日志、PID、退出码和COMPLETE.json。用锁避免重复启动同一任务。
- 下载仅对公开及已授权地址进行；使用可续传临时文件，完成后校验大小、压缩文件结构及SHA256，再改名。网络失败有界重试，不无限循环或静默遗漏样本。
- 同时至多2项计算任务；每项建议4个CPU、先按16 GB内存预算安排，实际按数据体量调整。资源值为调度预算，不是已测运行耗时。
- 调度器采用依赖完成事件或原生进程等待；无中间判断需求时不反复在对话层查询。确需轮询时由一个持久监控程序每15分钟检查进程。进程消失不视为成功，必须核查退出码、完成标记和结果完整性。
- 发现输入/注释错配、物种错误、患者重复或模型秩缺失时该任务失败并保留日志；独立任务继续。证据不足使用insufficient_data状态，不把它写成生物学阴性。
- 预计首批交付：V01稳健性表和掩膜图、V02/V04数据资格表、下一阶段是否可执行的明确判定。此文档为执行规划；尚未启动V01–V05。

## 数据入口

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231559
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285991
- https://doi.org/10.6084/m9.figshare.24803571
- https://doi.org/10.6084/m9.figshare.24803616
- https://doi.org/10.5281/zenodo.7737590
'''
(B/'后台验证计划.md').write_text(plan,encoding='utf-8-sig')
jobs=[{'id':'V01','name':'strict_spatial_masks','depends_on':[],'input_status':'available','status':'planned'},{'id':'V02','name':'independent_cell_reference','depends_on':[],'input_status':'metadata_available; raw_data_and_eligibility_pending','status':'planned'},{'id':'V03','name':'composition_and_state_separation','depends_on':['V02'],'status':'planned'},{'id':'V04','name':'external_spatial_replication','depends_on':[],'input_status':'public_download_links_verified; clinical_eligibility_pending','status':'planned'},{'id':'V05','name':'public_perturbation_response','depends_on':['V01','V02','V03','V04'],'condition':'earlier evidence warrants mechanism followup','status':'planned'}]
(B/'tasks.json').write_text(json.dumps({'status':'planned_not_launched','max_concurrent_compute_jobs':2,'jobs':jobs},ensure_ascii=False,indent=2),encoding='utf-8')
print('Background validation plan saved; not launched')
