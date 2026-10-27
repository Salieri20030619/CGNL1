from pathlib import Path,PureWindowsPath
import json,csv,hashlib,datetime,zipfile,xml.etree.ElementTree as E
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'material_search_restart_20261024'
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert json.loads(f.read_text())==z
verified=[]
for label in ['BCM','UMich']:
 for r in json.loads((D/(label+'_extraction_receipts.json')).read_text()):
  f=D/label/PureWindowsPath(r['local_file']).name;assert sha(f)==r['sha256'];verified.append(dict(path=str(f),sha256=sha(f),bytes=f.stat().st_size,zip_member_crc_verified=r['zip_member_crc_verified']))
matched=[]
for f in (D/'BCM').glob('PDAC_phospho*'):
 old=P/'raw'/f.name;assert sha(f)==sha(old);matched.append(dict(current=str(old),official=str(f),sha256=sha(f),byte_identical=True))
f=next((D/'UMich').glob('*.tsv'));targets=[]
with f.open() as h:
 r=csv.DictReader(h,delimiter='\t');headers=r.fieldnames
 for row in r:
  if row['Gene']=='ENSG00000128849.10':
   targets.append(dict(index=row['Index'],site=row['Index'].rsplit('_',1)[1],peptide=row['Peptide'],MaxPepProb=row['MaxPepProb'],ReferenceIntensity=row['ReferenceIntensity'],missing=sum(row[x]=='NA' for x in headers[5:]),sample_columns=len(headers)-5))
s198=next(x for x in targets if x['site']=='S198');s297=next(x for x in targets if x['site']=='S297');assert s198['peptide']=='PSNsQPTsPSLEDPAK'
human=''.join(x.strip() for x in (P/'raw/CGNL1_Q0VF96_380.fasta').read_text().splitlines() if not x.startswith('>'))
peptide=s198['peptide'];start=human.find(peptide.upper());assert start>=0
positions=[start+i+1 for i,x in enumerate(peptide) if x.islower()];assert positions==[198,202]
doc=next((D/'BCM').glob('*.docx'))
with zipfile.ZipFile(doc) as z:
 xml=E.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
 paragraphs=[''.join(x.itertext()) for para in xml.findall('.//w:p',ns) for x in []]
 paragraphs=[''.join(t.text or '' for t in para.findall('.//w:t',ns)) for para in xml.findall('.//w:p',ns)]
(D/'BCM_README_extracted.txt').write_text('\n'.join(paragraphs),encoding='utf-8')
audit=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files_verified=verified,BCM_current_matrix_match=matched,UMich_CGNL1_rows=targets,S198_annotated_phosphoresidues=positions,
 source_qualification_update='Exact current PDAC files match official BCM archive byte-for-byte;upstreamUMichsite-report recovered. Per-observation localization and aggregation lineage still incomplete.',
 new_interpretation_issue='S198 source peptide includes both S198 and S202 modification annotations. A site-level label does not establish isolated S198 occupancy or functional effect.',
 MaxPepProb_boundary='Do not equate maximum peptide probability with per-site/per-sample localization probability.',
 version_discrepancy='BCM README listsMSFragger3.4/Philosopher4.0.1;UMichREADME lists3.2/3.3.12 andTMTIntegrator2.1.5/2.4.0. Matchnumerical observations before claiming exact software lineage.',
 C04='Newpublisherroute403;massiveHTTPSdirectorySSLfailure;no newC04matrix or RAW obtained',new_tests=0,baseline_results_unchanged=True)
save(D/'MATERIAL_QUALIFICATION.json',audit)
text='''CGNL1材料重新搜寻结果（2026-10-24）
本轮取得新的官方源材料，C01来源缺口已部分补齐，C04尚未解除。

已下载并验收
通过PDC官方公开代码确认文件目录API与官方下载接口；使用HTTP Range按需读取ZIP成员，避免下载无关泛癌数据。6个目标文件均通过ZIP成员CRC与本地/服务器SHA256校验：
1. BCM版PDAC正常组织磷酸化矩阵，22,119,041字节。
2. BCM版PDAC肿瘤磷酸化矩阵，48,706,052字节。
3. BCM proteomics README（DOCX）。
4. BCM data-freeze文件说明（XLSX，已归档，未用于新统计）。
5. UMich README。
6. UMich PDAC单位点汇总表，143,365,858字节。
其中两份BCM矩阵与项目原有输入逐字节SHA256一致。现在可以确认当前PDAC输入来自官方BCM统一版，不再仅依据文件名或一般方法推断。官方目录明确将BCM版与Sinai预处理版分列；但不能仅据此宣布所有上游处理或逐观测插补状态已经闭合。

新发现及其解释
UMich表已找到S198和S297源行。S198肽段为PSNsQPTsPSLEDPAK，按人CGNL1参考序列定位，两个小写标记分别在S198和S202。S297行肽段为RSsSSSTTPTSANSLYR;SsSSSTTPTSANSLYR。
这说明当前S198标签至少关联一个双磷酸化肽段。现有结果仍是锁定的位点汇总信号，但不能直接解释为孤立S198占有率或功能变化。需要继续核对单位点汇总对多磷酸化肽段的处理，以及哪些谱图/样本实际贡献了该行。
MaxPepProb不是逐位点/逐样本定位概率，不能因为显示1就称100%位点定位。
两个包的README还列出不同软件版本：BCM为MSFragger3.4、Philosopher4.0.1；UMich为MSFragger3.2、Philosopher3.3.12及两个TMT-Integrator版本。这一版本差异已保留，下一步需通过患者/样本列与数值对应核对处理关系，不能只凭上下游描述宣称完全同一版本。

C04重新搜索
出版社原文入口仍403；新尝试的MassIVE HTTPS目录出现SSL连接失败。本轮没有取得有效C04矩阵/RAW，没有将错误页当数据，也未重复原限流RAW请求。没有对外联系作者。

下一步
优先利用新取得的UMich表完成S198/S297样本列、缺失模式和BCM数值的对应；核对多磷酸化肽段汇总规则及定位字段，必要时再按需提取肽段/多位点表。继续寻找C04可用公开表格入口。现有模型不变，不因新材料自动重选终点；Goal已由用户恢复active。
证据目录：continuous_exploration_20261023_v1/material_search_restart_20261024。
'''
report=D/'CGNL1_MATERIAL_SEARCH_RESTART_20261024.txt';report.write_text(text,encoding='utf-8');assert report.read_text()==text
save(D/'ACCEPTANCE.json',dict(utc=audit['utc'],verified_files=len(verified),BCM_matrix_identity=True,source_site_rows=len(targets),report_sha256=sha(report),new_tests=0,whole_project_complete=False))
for f in [O/'STATUS.json',O/'GOAL_PLAN.json',P/'metadata/PROJECT_CHECKPOINT_CONTINUOUS_20261023.json']:
 s=json.loads(f.read_text());s.update(goal_status='active',material_search_restart=str(D),latest_material_report=str(report),whole_project_complete=False);save(f,s)
print(json.dumps(dict(report=str(report),sha256=sha(report),files=len(verified),CGNL1_rows=len(targets),S198_phosphoresidues=positions,BCM_identity=True)))
