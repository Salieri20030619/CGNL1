from pathlib import Path
import json,hashlib,zipfile,re,html,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'functional_sources'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,z):
 p.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');assert json.loads(p.read_text())==z
expected={'PKAKO_author_table.html':'ab61638dcb8b51c1ab952d15b09887df7c86303b2b902192b1fc95418dda0551','PKAKO_source.xlsx':'fd70d3c5f5e15b9f224cb6461171119920d58b69607458c39f2cf18088bec4f2','PKAKO_article.pdf':'ca5a2584c5eb6c6950e546a5b5460fff856bc468f5b1d3df923c4066ac34ca28'}
for n,h in expected.items():assert sha(D/n)==h
with zipfile.ZipFile(D/'PKAKO_source.xlsx') as z:assert z.testzip() is None
ns={};exec((P/'audit/CGNL1_FOLLOWUP_20261023_v1/followup.py').read_text().split('baseline=read')[0],ns)
rows=[]
for i,r in ns['rows'](D/'PKAKO_source.xlsx','sort by gene symbol'):
 if r.get('B'+str(i))=='Cgnl1':rows.append(dict(row=i,uniprot=r['A'+str(i)],site=r['C'+str(i)],motif=r['E'+str(i)],phospho_log2_source=r['F'+str(i)],source_P=r['G'+str(i)],totalprotein_log2_source=r['H'+str(i)]))
assert len(rows)==10
target=next(x for x in rows if x['site']=='S298');null=next(x for x in rows if x['site']=='S284')
assert target['phospho_log2_source']=='-1.73 ± 0.28' and null['phospho_log2_source']=='0.09 ± 0.31'
alignment=json.loads((D/'HUMAN_MOUSE_ALIGNMENT.json').read_text())
assert all(x['mapping']['298']['human_position']==297 and x['mapping']['284']['human_position']==283 for x in alignment['records'])
s=(D/'PKAKO_author_table.html').read_text(encoding='utf-8');clean=lambda t:re.sub(r'\s+',' ',html.unescape(re.sub('<[^>]+>',' ',t)))
htmlrows=[clean(x) for x in re.findall(r'<tr\b[^>]*>.*?</tr>',s,re.S|re.I) if 'Cgnl1' in x];assert len(htmlrows)==10
assert any('-1.73 ± 0.28' in r and '8.72E-03' in r for r in htmlrows)
checks=[]
for wave in ['wave1','wave2','wave3_source_purity']:
 a=json.loads((O/wave/'ACCEPTANCE.json').read_text())
 for item in a['artifacts']:
  f=Path(item['path']);assert sha(f)==item['sha256'];checks.append(str(f))
utc=datetime.datetime.now(datetime.timezone.utc).isoformat()
z=dict(utc=utc,doi='10.1073/pnas.1709123114',source_url='https://esbl.nhlbi.nih.gov/Databases/PKAKO/',model='mouse kidney mpkCCD cells;three independent PKA-null/control clone pairs',perturbation='CRISPR-Cas9 deletion of both PKA catalytic genes',sources=expected,rows=rows,source_P_not_new_test=True,source_FDR_not_in_table=True,uncertainty_symbol_preserved_without_assuming_SD_or_SEM=True,
 interpretation='PKA deletion is associated with a site-selective Cgnl1 phospho response:mouseS298/humanS297 decreases;mouseS284/humanS283 mean response nearzero and sourceP0.657;totalproteinmean0.08log2.',
 boundaries=['Not a digestive disease cohort','No CGNL1 mutant/rescue or direct kinase-substrate assay established here','Three clone pairs not patient replication','No new P for site difference,protein adjustment or cross-species extrapolation from summary statistics','PKA knockout can act through indirect kinases or clonal adaptation','Do not call all CGNL1 phosphosites decreased;all10rows retained'],new_tests=0,html_workbook_agree=True,baseline_artifacts_hash_verified=len(checks))
save(D/'PKAKO_CGNL1_QUALIFICATION.json',z)
report='''CGNL1持续探究：PKA扰动线索与候选方向更新
本次为公开源数据资格核查，没有新增统计检验或更改既有结果。

新证据：作者公开的PKA双催化亚基敲除源表
研究：Isobe等，PNAS 2017，DOI 10.1073/pnas.1709123114。
设计：小鼠肾集合管上皮mpkCCD细胞，3对独立PKA-null与相应对照克隆。
已取得HTML作者表、XLSX源表和论文PDF，完成文件哈希、XLSX CRC及HTML/XLSX关键行一致性检查。
源表S298：磷酸化log2(KO/对照)为−1.73 ± 0.28；S284为+0.09 ± 0.31；总Cgnl1为+0.08 ± 0.23。保留原表±符号，未擅自指定SD/SEM。源表对应P值分别为0.008721047152391914和0.6572453642563502，属于原作者统计，不是本项目新增检验；表内未提供相应FDR。完整10条Cgnl1记录均保存，未仅保留阳性。
既有全长比对在24条受检最优比对中一致支持：小鼠S298对应人S297，小鼠S284对应人S283。

对研究方向的实际影响
这条证据支持将“不同位点受不同上游调控”作为可检验问题。它为C01的S297端及C03的S283/S297差异增加了扰动先验，说明不能把CGNL1全部磷酸化视为一个统一升降指标。
它没有证明PKA直接磷酸化S297，没有证明该位点改变屏障或肿瘤表型，也不能作为消化道疾病独立重复。来源是敲除克隆，仍有间接激酶调控和克隆适应等解释。未用汇总均值计算新的位点差异P值或蛋白校正P值。
可执行的后续实验候选：先在相关消化系统上皮模型中确认S198、S283、S297、S486的可测性；在总CGNL1和细胞状态可控的条件下，比较PKA急性扰动与现有Mer/Gas6线索；只有确认位点响应后，再开展相应位点突变/救援及屏障或细胞连接功能实验。这是候选实验顺序，不代表已完成实验，也不预设所有疾病共用同一激酶。

候选比较保持开放
C01胰腺癌多位点与C02 S486相对模式仍为条件性共同主候选；C03保留HCC异质性比较；C04保留炎症材料资格。PKA文献不能消除H02阴性、H03插补/重叠限制或C01同疾病独立验证缺口。
本轮复核了前三轮验收清单所列全部结果/代码/输入文件的哈希，均与验收时一致。C04转换工具已完成独立安装及版本检查，但RAW访问仍受限，没有把工具准备视为原始数据试验完成。本轮未重复限流链接。
Goal保持active，接下来继续精确来源核查及可获得材料的单文件试验；最终主线仍需依据独立重复与功能可行性判定。
'''
f=O/'CGNL1_PKA_FUNCTIONAL_UPDATE_v1.txt';f.write_text(report,encoding='utf-8');assert f.read_text()==report
save(D/'PKAKO_ACCEPTANCE.json',dict(utc=utc,source_hashes=expected,html_workbook_agree=True,baseline_hash_count=len(checks),report=str(f),report_sha256=sha(f),new_tests=0,goal_status='active'))
state=json.loads((O/'STATUS.json').read_text());state.update(utc=utc,PKA_functional_lead='mouseS298/humanS297 responds toPKAdKO;S284/humanS283 source-null;notGIor directkinaseproof',latest_functional_report=str(f));save(O/'STATUS.json',state)
print(json.dumps(dict(report=str(f),sha256=sha(f),source_rows=len(rows),baseline_hash_checks=len(checks))))
