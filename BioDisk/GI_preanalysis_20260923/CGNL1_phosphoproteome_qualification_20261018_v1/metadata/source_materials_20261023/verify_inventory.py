from pathlib import Path
import json,hashlib,datetime,zipfile,collections
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
M=P/'metadata/source_materials_20261023'
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_source_materials_20261023')
B=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576')
V=P/'audit/CGNL1_SOURCE_CONSISTENCY_20261023_v2'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
inventory={}
def add(path,group,status,expected=None,verification=None):
 p=Path(path);h=sha(p)
 if expected:assert h==expected,(str(p),'hash mismatch')
 inventory[str(p)]=dict(path=str(p),group=group,status=status,bytes=p.stat().st_size,sha256=h,verification=verification or 'SHA256 recorded')
for r in read(M/'CDAP_PINNED_SOURCE_RECEIPTS.json'):
 assert r['exit']==0 and r['git_blob_pass'];add(r['path'],'CDAP_2019_workflow','new_download',r['sha256'],'pinned commit and Git blob SHA1 verified')
tree=read(M/'CDAP_2019_tree.json')['tree'];f=D/'CDAP_20190522/lib/scripts/dltools.sh';b=f.read_bytes()
assert hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==next(x['sha'] for x in tree if x['path']=='lib/scripts/dltools.sh')
add(f,'CDAP_2019_workflow','new_download',verification='Git blob verified; downloaded for provenance only, not executed')
for r in read(M/'H01_QC_TSV_RECEIPTS.json'):
 assert r['exit']==0 and r['md5_pass'] and r['size_pass'];add(r['file'],'H01_QC','new_download',r['sha256'],'official size and MD5 verified')
for exp in ['OEX00011296','OEX00011347']:
 f=M/(exp+'_verified_scope.json');r=read(f)['data'];assert r['total']==len(r['list'])==528
 assert all(x['expNo']==exp and x['security']=='Restricted' and x['accessible']==False for x in r['list'])
 add(f,'B02_raw_access_inventory','new_metadata_download',verification='correct experiment ID,528/528 rows, all Restricted')
for name in ['CDAP_2019_commits.json','CDAP_2019_tree.json']:
 add(M/name,'CDAP_version_metadata','new_metadata_download')
manifest=read(P/'metadata/source_search_593/B02_NODE_file_manifest.json')['data']['list']
mapping={'OED00658855':'B02_groups_NODE.xlsx','OED00658856':'B02_phospho_NODE.xlsx','OED00658860':'B02_protein_NODE.xlsx'}
for r in manifest:
 f=B/'source_search_593'/mapping[r['datNo']];b=f.read_bytes();assert len(b)==r['fileSize'] and hashlib.md5(b).hexdigest()==r['md5']
 with zipfile.ZipFile(f) as z:assert z.testzip() is None
 add(f,'B02_processed_and_channels','existing_verified',verification='publisher MD5, size and ZIP CRC rechecked')
f=B/'source_search_593/B02_TableS1.xlsx';expected=read(P/'audit/CGNL1_FOLLOWUP_20261023_v1/SOURCE_HASHES_BEFORE.json')[str(f)];add(f,'B02_clinical','existing_verified',expected)
for r in read(P/'metadata/resume_581/MEMBERS_COMPLETE585.json')['files']:
 add(r['file'],'H02_MaxQuant_phospho','existing_verified',r['sha256'],'prior member CRC plus SHA256 rechecked')
for r in read(P/'metadata/H02_GLOBAL_EXTRACTION589.json')['records']:
 add(r['saved'],'H02_MaxQuant_global','existing_verified',r['sha256'],'prior archive extraction validation plus SHA256 rechecked')
for r in read(P/'results/H01_582/RESULT.json')['source_files']:
 add(B/'H01'/r['file'],'H01_quantification_mapping','existing_verified',r['sha256'])
h03=read(V/'HCC_COHORT_CHRONOLOGY_CONSISTENCY.json')['H03'];add(h03['file'],'H03_clinical_dates','existing_verified',h03['sha256'])
prior=read(V/'CONSISTENCY_ACCEPTANCE.json')['artifacts'];wanted=['B02_supplement.pdf','NODE_samples.json','NODE_proteomic_experiments.json','NODE_getRelatedAnalysis.json','Zhou_Liver_Cancer_Proteome.qcmetrics.html','Zhou_Liver_Cancer_Phosphoproteome.qcmetrics.html']
for r in prior:
 if Path(r['path']).name in wanted:add(r['path'],'source_methods_metadata','existing_verified',r['sha256'])
f=P/'audit/CGNL1_FOLLOWUP_20261023_v1/sources/B02_Dong2022.pdf'
old=read(P/'audit/CGNL1_FOLLOWUP_20261023_v1/FOLLOWUP_ACCEPTANCE.json')['artifacts']
add(f,'B02_methods','existing_verified',next(x['sha256'] for x in old if x['path']==str(f)))
restricted=[]
for exp in ['OEX00011296','OEX00011347']:
 a=read(M/(exp+'_verified_scope.json'))['data']['list'];restricted.append(dict(experiment=exp,files=len(a),bytes=sum(x['fileSize'] for x in a),status='Restricted',downloaded=False))
missing=[dict(item='B02 uncollapsed channel-level quantitative matrix and final repeat aggregation rule',status='not_found_in_current_public_MS_analysis_materials',required_action='author/data-manager supplied processed table or processing code'),
 dict(item='B02 original spectra',status='access_restricted',required_action='repository-approved access required before download; no application submitted'),
 dict(item='H01/H02/H03 cross-study donor mapping',status='not_publicly_available_in_checked_materials',required_action='de-identified mapping or overlap statement from data manager'),
 dict(item='H01 exact executed workflow version and parameters',status='not_yet_linked',required_action='run-specific provenance; downloaded 2019 configuration is a comparison candidate, not an execution certificate')]
counts=dict(collections.Counter(x['status'] for x in inventory.values()))
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),download_scope='materials needed for A01/A04/A06; B02 modeling remains paused',
 public_materials_acquisition_complete_for_identified_sources=True,all_required_materials_available=False,
 new_material_root=str(D),counts=counts,inventory=list(inventory.values()),remaining=missing,restricted_raw=restricted,
 incomplete_archive=dict(path=str(D/'CDAP_galaxy_tools_20190522.zip.partial'),status='preserved_partial_not_accepted',replacement='25 required historical workflow/configuration files downloaded separately and Git-blob verified'),
 invalid_scope_queries=['OEX00011296_public_files.json','OEX00011347_public_files.json','NODE_getPublicDataList_scoped.json','NODE_getDataList_scoped.json'],
 invalid_scope_note='These exploratory API queries did not reliably scope the requested experiment and are not used as absence evidence.',
 scripts_executed=False,clinical_models_run=0,whole_project_complete=False)
f=M/'MATERIAL_DOWNLOAD_MANIFEST.json';f.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==result
report=f'''CGNL所需材料下载清单（2026-10-23）

本轮新增27份可用材料：25份固定commit的2019年CDAP报告工作流/版本配置，2份H01官方批次QC明细TSV。另保存B02两实验完整受限文件目录及CDAP版本/文件树记录。
新增材料目录：{D}
统一清单：{f}

已有材料已复核并纳入清单：B02通道表、蛋白组和磷酸化最终表及临床表；H01定量矩阵、样本映射和官方QC报告；H02两组学MaxQuant参数、结果和临床映射；H03临床日期表；B02主文、补充图与仓库样本元数据。大型文件保留原路径，未重复下载。

验收：CDAP逐文件核对Git blob哈希；H01新增TSV核对官方大小与MD5；B02三份仓库工作簿再次核对官方MD5、大小和ZIP CRC；既有文件核对SHA256。统一清单共{len(inventory)}个可用文件。没有运行所下载的外部脚本，没有启动B02模型。

仍未取得：
1. B02未合并通道级定量及作者重复合并/选择规则，当前公开MS分析材料没有找到。
2. B02原始谱图：蛋白组528份、磷酸化组528份，合计1,589,335,200,757字节，均Restricted且当前不可访问。需仓库批准访问，未代发申请。
3. HCC跨研究患者映射/重叠声明，公开材料尚未提供。
4. H01实际运行的精确版本与参数绑定。已下载CPTAC-GALAXY-1.3.35相关配置用于比对，不等于已证明H01使用此版本。

完整GitHub压缩包超时部分保留为.zip.partial，不计入成功材料；所需25个小文件已单独下载并验收。下载任务在已识别公开来源范围内完成，但不能声称全部证据缺口已经补齐。
'''
rf=M/'CGNL_MATERIAL_DOWNLOAD_RECEIPT_20261023.txt';rf.write_text(report,encoding='utf-8');assert '\ufffd' not in rf.read_text(encoding='utf-8')
accept=dict(utc=result['utc'],manifest_sha256=sha(f),report_sha256=sha(rf),files=len(inventory),counts=counts,all_inventory_hashes_pass=True,all_required_materials_available=False,B02_modeling_paused=True)
(M/'DOWNLOAD_ACCEPTANCE.json').write_text(json.dumps(accept,indent=2)+'\n',encoding='utf-8')
print(json.dumps(accept))
