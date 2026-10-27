from pathlib import Path
import json,shutil,hashlib
O=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3');D=Path(r'C:\Users\31937\Downloads\CGNL1立项建议_20260924');D.mkdir(exist_ok=True)
s='''# CGNL1肝转移空间研究：数据可行性与创新定位

日期：2026-09-24。状态：已完成文献、患者映射、数据可获取性及部分平台CGNL1覆盖审查；尚未开展肝转移交界空间差异或邻域关联分析。本轮无推断性统计结论。

## 判断

空间转录组能够帮助分离组织来源与局部癌细胞状态，有实质分析价值。创新性应落在CGNL1的具体细胞来源、邻域条件和功能机制。加入一种技术或重新观察肝适应现象，本身不足以形成新机制。

空间转录组在2016年已有Science原始技术论文（https://pubmed.ncbi.nlm.nih.gov/27365449/）。2020年Cell Research的肝转移重编程研究主要采用转录组、ChIP-seq和功能实验，未使用本次设想的空间转录组定位。2024年Nature已有肝转移交界空间分析与肝细胞来源Plexin B2—肿瘤semaphorin—KLF4功能研究，必须作为近邻先行工作核对：https://www.nature.com/articles/s41586-024-07715-3 。其人空间数据是1名患者的2张非连续切片，不能记为2名患者。

2026年British Journal of Cancer已针对肝转移生长模式划分肝侧、肿瘤侧及包膜两侧空间区域，因此“分析肝—肿瘤交界区”同样不是空白领域：https://www.nature.com/articles/s41416-026-03567-y 。

## 已核实的数据

| 数据 | 样本与现有材料 | CGNL1及用途判断 |
|---|---|---|
| GSE294385 | 既有作者样本映射中，常规Visium肝转移24张切片来自11名患者；5名同时有原发CRC切片；包含宏/微转移和部分复发病灶。已具有全队列作者注释和患者映射；本地此前仅下载并分析8张原发切片。 | 优先主队列。需下载肝转移矩阵并核查每例CGNL1检出、交界是否完整、病理区域及治疗信息。11人不是24人；有效比较患者数待QC。逐患者结肠/直肠部位尚未明确，当前按CRC队列，不能直接声称11例纯结肠癌。 |
| GSE217414 | 4名患者、4张肝转移Visium；本轮已下载4个表达H5。患者1乙状结肠，2高位直肠，3低位直肠，4直乙交界；全部术前治疗。 | 四张均测得CGNL1；检出spot比例26.43%、36.43%、12.48%、19.51%。只代表检测可行性，不代表癌细胞定位或交界富集。严格结肠子集只有1人，仅能补充描述。空间坐标/组织图可公开下载，本轮未下载。 |
| Niche-DE / Genome Biology 2024 | 研究整合5名患者；3例来自OEP001756，另外2例原始Visium分别开放在Figshare 24803571、24803616，文件约147/367 MB。已下载API文件清单。 | 外部候选；先核对结肠/直肠来源、治疗、与其他队列患者重复和CGNL1检出。未下载全矩阵，不提前宣称外部验证成功。 |
| Nature 2024 mosaic liver | 人空间1名/2切片；Zenodo公开Visium.zip约303 MB，另有单细胞、多组学和bulk干预数据。本轮已核实API及下载入口。 | 可用作交界定位与上游机制的探索参考。人空间不能作多患者复现；干预数据的物种、处理和生物学重复另行审查，不能把小鼠结果直接当人结肠癌验证。 |
| GSE335552 / Xenium | 3张肿瘤和1张匹配邻近肝；本轮实际读取全部4个cell_feature_matrix.h5。 | 四个矩阵特征均无CGNL1，不能直接验证CGNL1。未测量不等于阴性，不通过插补把CGNL1制造成测量结果。 |
| BJC 2026 / ISS | 19人、26 ROI的设计；QC后Panel 1为19人25 ROI，Panel 2为17人22 ROI。已下载并检索官方两套探针表。 | 两套探针表均无CGNL1，不能直接验证该基因；可参考病理区域和肝—肿瘤界面设计。患者为推断单位，即便原文部分分析另作处理，本项目仍不把同患者ROI当独立患者。 |

GSE294385作者研究为Spatial Multi-Omics Landscape of Colorectal Cancer Macro- and Micrometastases（https://pmc.ncbi.nlm.nih.gov/articles/PMC13390483/）。常规Visium公开；Visium HD等资料所在作者Zenodo最新记录为restricted，未发起申请，不将受限数据写成已经可获取。主队列使用公开部分即可启动可行性分析。

## 建议的分析主线

1. 确定细胞来源：结合H&E、作者病理分区、独立单细胞参考与肿瘤标记，区分CGNL1来自癌细胞、肝实质、血管内皮或成纤维细胞。Visium spot混合细胞，不能将一个spot直接称为一个癌细胞。不能只用EPCAM判定恶性。
2. 在同一患者内部比较：先在不查看CGNL1的条件下确定交界，再比较肿瘤侧近肝区域与肿瘤内部；同时分析肝侧近端/远端。按物理坐标标定距离，主距离带与敏感性范围在结果分析前固定，不能选显著阈值。显式控制spot组成和测序深度，并评估邻近转录本扩散。
3. 研究状态条件：在病理确证的宏转移与微转移、包膜与直接接触区域中考察关联。有相应区域及足够患者时才比较；截面病灶大小不等于纵向定植进程，不能把微转移直接当作所有宏转移的前驱。
4. 筛选上游机制：从相邻肝细胞/基质信号到癌细胞受体及下游状态寻找一致链条。CGNL1是胞内连接支架，不能硬套为分泌配体或受体。空间接近、配体受体推断和横断面相关均不是因果；公开干预数据可提供更直接的生信支持。
5. 患者级复核：切片内聚合区域效应，再按患者汇总；重复切片、复发切片和ROI不增加独立n。小于10人仅描述性/探索性；达到10名有效患者再报告效应、P与适当同类假设家族FDR。严格结肠与CRC支持证据分开报告。

## 立项价值与判定条件

优先检验的具体问题是：CGNL1相关的癌细胞连接状态，是否集中于某类肝—肿瘤接触邻域，并在独立患者中复现，进而指向一个可扰动和救援的上游调控过程。

如果只有肝/基质区域表达高，研究定位应改为该来源及其局部作用；若癌细胞内无可复现空间关系，不把bulk阳性强接到空间机制。若患者内关系、外部复现与公开扰动结果支持同一链条，再围绕该链条安排客户湿实验。这是提高文章竞争力的依据，不能将“加了空转”作为约10分发表保证。

## 原始入口与审查文件

- GSE294385：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294385
- GSE217414：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217414
- Niche-DE：https://link.springer.com/article/10.1186/s13059-023-03159-6
- Figshare：https://doi.org/10.6084/m9.figshare.24803571 ，https://doi.org/10.6084/m9.figshare.24803616
- Nature 2024数据：https://doi.org/10.5281/zenodo.7737590
- GSE335552：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE335552

本地材料：C:/Users/31937/CRC_routes_20260923/spatial_v3。服务器：/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/spatial_v3。

本轮完成标志：8个H5均成功解析；ISS两套探针表完成精确基因名/Ensembl检查；4份外部API下载入口审查通过。CGNL1表达覆盖详见results/CGNL1_platform_coverage.tsv。本轮尚未计算交界梯度和患者级空间效应。
'''
p=O/'空间方向可行性与创新定位.md';p.write_text(s,encoding='utf-8-sig');shutil.copy2(p,D/p.name)
manifest=[]
for p in (O/'raw').glob('*.h5'):
 manifest.append({'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()})
(O/'results/INPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
(O/'results/AUDIT_COMPLETE.json').write_text(json.dumps({'status':'complete','scope':'dataset eligibility, literature and CGNL1 measurement coverage only','spatial_association_analysis':'not_started','H5_files':8,'ISS_panels':2,'main_cohort_liver_sections':24,'main_cohort_unique_patients':11},indent=2),encoding='utf-8')
print('Saved',D/'空间方向可行性与创新定位.md')
