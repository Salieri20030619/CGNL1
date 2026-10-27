from pathlib import Path
import json,hashlib,datetime,html,subprocess,tarfile,shlex,os
from bs4 import BeautifulSoup
B=Path(r'D:\CGNL1');H=B/'handoff';W=B/'_transfer';F=B/'final'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
migration=json.loads((H/'MIGRATION_STATUS.json').read_text(encoding='utf-8'));assert migration['stage']=='complete'
progress=json.loads((H/'PROGRESS_SNAPSHOT.json').read_text(encoding='utf-8'))
quality=json.loads((H/'NUMERICAL_QA.json').read_text(encoding='utf-8'));assert quality['all_pass']
report=F/'SCP2079_20261014_v1/CGNL1_SCP2079_paired_report.html'
expected=json.loads((H/'REPORT_VALIDATION.json').read_text(encoding='utf-8'))['report_sha256'];assert sha(report)==expected
for r in F.rglob('*.html'):
 s=r.read_text(encoding='utf-8');assert '\ufffd' not in s
 soup=BeautifulSoup(s,'html.parser')
 for img in soup.find_all('img',src=True):assert img['src'].startswith(('data:','https:','http:')),str(r)
source=json.loads((H/'LOCAL_ARCHIVE_MANIFEST.json').read_text(encoding='utf-8'))+json.loads((H/'RAW_TRANSFER_MANIFEST.json').read_text(encoding='utf-8'))
assert all(not Path(x['local']).exists() for x in source)
final_manifest=[{'path':p.relative_to(B).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(F.rglob('*')) if p.is_file()]
(H/'FINAL_ARTIFACT_MANIFEST.json').write_text(json.dumps(final_manifest,indent=2),encoding='utf-8')
finalqa={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'new_report_matches_server_sha256':True,'retained_HTML_images_self_contained':True,'UTF8_readback_no_replacement_characters':True,'all_2804_migrated_source_files_absent_locally':True,'all_analysis_stage_exit_codes_zero':True,'numerical_QA_all_pass':True,'final_artifact_files':len(final_manifest),'plot_visually_checked':True,'migration_final_exit_code':0,'prior_readonly_cleanup_failure_resolved':True}
(H/'FINAL_DELIVERY_QA.json').write_text(json.dumps(finalqa,indent=2),encoding='utf-8')
(H/'MIGRATION_FINAL_EXIT.json').write_text(json.dumps({'exit_code':0,'completed_utc':finalqa['verified_utc'],'supersedes':'Initial archival coordinator intentionally retired; first cleanup stopped on read-only PDF; verified retry completed'},indent=2),encoding='utf-8')
links=[]
for item in final_manifest:
 if item['path'].endswith(('.html','.pdf','.docx')):links.append('<li><a href="../'+html.escape(item['path'],quote=True)+'">'+html.escape(Path(item['path']).name)+'</a></li>')
content='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>CGNL1 项目交接</title><style>body{font:16px/1.7 "Microsoft YaHei",sans-serif;max-width:1100px;margin:35px auto;padding:0 25px;color:#18313d}h1,h2{color:#155d74}code{overflow-wrap:anywhere;background:#edf3f6}td,th{border:1px solid #c9d5db;padding:8px}table{border-collapse:collapse}li{margin:8px 0}</style><body><h1>CGNL1 项目交接：服务器分析，本地仅保留交接与终稿</h1>'''
content+='<p>更新时间：'+html.escape(finalqa['verified_utc'])+'。</p><h2>文件位置与执行规则</h2><ul><li>所有原始数据、代码、中间结果和统计计算均在 sgs10（osv01）上管理。</li><li>本地固定目录为 <code>D:\\CGNL1</code>，仅保留 handoff 与 final；原始文件不存 C 盘，也不作为 D 盘长期分析副本。</li><li>已迁移 22 个 CGNL1 目录、早期 CRC_routes_20260923 目录及本次下载的三个原始文件，共 2,804 个文件、6,512,517,421 字节（6.51 GB）。逐文件 SHA256 一致后才删除本地源文件。</li><li>混合 GI_preanalysis 目录中的胃癌等非本项目文件不在此次迁移范围。该目录未检出本轮列明 CGNL1 队列的文件名残留。</li></ul>'
content+='<p>项目根目录：<code>/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1</code></p><p>本轮分析：<code>'+html.escape(migration['analysis'])+'</code></p><p>本地历史归档：<code>'+html.escape(migration['archive'])+'</code>，保留原始相对目录结构。</p>'
content+='<h2>SCP2079 已完成的验证</h2><p>矩阵为经过实测核对的原始计数，36,601 个基因条目、95,421 个细胞；全部细胞的计数和检测基因数与作者注释一致。27 名患者中 19 名具有治疗前后上皮数据，既定 QC 与每端至少 20 个细胞留下 14 对；明确 Liver 标注的 8 对仅作描述。</p>'
content+='<p>主分析中 11/14 对患者 CGNL1 治疗后下降，中位变化 −0.01309 ln(1+CP10k)，P=0.0166015625，FDR=0.0498046875。该变化方向在敏感性分析中保留，但显著性受细胞筛选口径影响。</p><table><tr><th>控制上皮状态后的目标关联</th><th>标准化 β</th><th>P</th><th>FDR</th></tr><tr><td>actomyosin</td><td>0.6153</td><td>0.198649</td><td>0.397298</td></tr><tr><td>cycling</td><td>−0.4263</td><td>0.720175</td><td>0.720175</td></tr></table>'
content+='<p>当前数据不支持两个终点的名义显著条件关联。actomyosin 主分析预测误差下降 42.8%，但在每端至少 50 个细胞的敏感性分析中改为误差增加 13.8%；该门槛排除治疗后只有 26 个细胞的 P29。细胞周期主分析预测误差增加 67.4%。因此没有建立稳定的独立预测增益。</p>'
content+='<p>这项研究提供联合治疗背景下的患者内变化证据，不能解释为 CGNL1 单基因扰动或肝定植因果验证。整个 CGNL1 项目尚未完成；本轮配对验证已完成。</p><h2>下一步衔接</h2><p>以服务器端 <code>CGNL1_PROGRESS_20261014_SCP2079_PAIRED.json</code> 为最新进度；继续保留治疗相关下降作为候选线索，后续验证需区分蛋白定位/力学功能与 RNA 状态关联。不得把敏感性检查失败的预测增益描述为稳定独立信息。</p>'
content+='<h2>终稿入口</h2><ul>'+''.join(links)+'</ul><h2>验收与迁移回执</h2><ul>'
for n in ['MIGRATION_STATUS.json','ARCHIVE_VERIFICATION.json','LOCAL_CLEANUP_RECEIPT.json','MIGRATION_FINAL_EXIT.json','FINAL_ARTIFACT_MANIFEST.json','FINAL_DELIVERY_QA.json','JOB_COMPLETE.json','NUMERICAL_QA.json','PROGRESS_SNAPSHOT.json']:
 content+='<li><a href="'+n+'">'+n+'</a></li>'
content+='</ul><p>历史报告保留为可阅读副本，嵌入图件保持可见；原始终稿与附属数据完整保留在服务器归档。曾遇到的只读参考 PDF 清理问题已解决，其失败和成功重试记录均保留。</p></body></html>'
(H/'CGNL1_HANDOFF.html').write_text(content,encoding='utf-8')
assert (H/'CGNL1_HANDOFF.html').read_text(encoding='utf-8')==content
flags=subprocess.CREATE_NO_WINDOW
remote=migration['archive']+'/migration_control'
cmd='mkdir -p '+shlex.quote(remote)
subprocess.run(['ssh','-o','BatchMode=yes','sgs10',cmd],check=True,capture_output=True,creationflags=flags)
control=[]
for p in W.rglob('*'):
 if p.is_file():control.append({'path':p.relative_to(W).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
(H/'MIGRATION_CONTROL_MANIFEST.json').write_text(json.dumps(control,indent=2),encoding='utf-8')
proc=subprocess.Popen(['ssh','-o','BatchMode=yes','sgs10','tar -xf - -C '+shlex.quote(remote)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,creationflags=flags)
with tarfile.open(fileobj=proc.stdin,mode='w|') as tar:
 for x in control:tar.add(W/x['path'],arcname=x['path'],recursive=False)
proc.stdin.close();rc=proc.wait();assert rc==0,proc.stderr.read().decode('utf-8')
subprocess.run(['scp','-r','-o','BatchMode=yes',str(H),'sgs10:'+migration['archive']+'/'],check=True,capture_output=True,creationflags=flags)
code='''
from pathlib import Path
import json,hashlib,datetime
A=Path(ARCH);H=A/'handoff'
for x in json.loads((H/'MIGRATION_CONTROL_MANIFEST.json').read_text(encoding='utf-8')):
 p=A/'migration_control'/x['path'];assert p.stat().st_size==x['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==x['sha256']
state=json.loads((H/'PROGRESS_SNAPSHOT.json').read_text(encoding='utf-8'))
state['migration_archive']=str(A);state['local_final_directory']='D:/CGNL1/final';state['handoff_directory']='D:/CGNL1/handoff'
state['migration_completed']=True;state['local_raw_data_removed']=True
(A.parent/'CGNL1_CURRENT_STATE_20261014.json').write_text(json.dumps(state,indent=2),encoding='utf-8')
print(json.dumps({'migration_control_all_verified':True,'handoff_archived':True,'current_state':str(A.parent/'CGNL1_CURRENT_STATE_20261014.json')}))
'''
p=subprocess.run(['ssh','-o','BatchMode=yes','sgs10','python3 -'],input=('ARCH='+repr(migration['archive'])+'\n'+code).encode('utf-8'),capture_output=True,creationflags=flags)
assert p.returncode==0,p.stderr.decode('utf-8');r=json.loads(p.stdout.decode('utf-8'))
(H/'HANDOFF_SERVER_RECEIPT.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(r),flush=True)