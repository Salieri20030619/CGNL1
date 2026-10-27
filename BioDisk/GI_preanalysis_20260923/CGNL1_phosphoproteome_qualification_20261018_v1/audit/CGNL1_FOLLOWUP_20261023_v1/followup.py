from pathlib import Path
import json, hashlib, statistics, datetime, collections, zipfile, xml.etree.ElementTree as ET, math

P = Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
A = P/'audit/CGNL1_FIRST_REVIEW_20261022_v1'
O = P/'audit/CGNL1_FOLLOWUP_20261023_v1'
D = Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576/source_search_593')
O.mkdir(exist_ok=True)
UTC = datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(p): return json.loads(p.read_text(encoding='utf-8'))
def save(name, value):
    f=O/name
    f.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    assert read(f)==value
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

# Read workbook XML without changing or re-saving any source workbook.
NS={'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def sheets(z):
    rel=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    targets={r.attrib['Id']:r.attrib['Target'] for r in rel}
    result={}
    for s in ET.fromstring(z.read('xl/workbook.xml')).findall('s:sheets/s:sheet',NS):
        t=targets[s.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
        result[s.attrib['name']]=t.lstrip('/') if t.startswith('/') else 'xl/'+t
    return result
def rows(path, sheet, wanted=None):
    with zipfile.ZipFile(path) as z:
        ss=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            for _,e in ET.iterparse(z.open('xl/sharedStrings.xml'),events=['end']):
                if e.tag.endswith('}si'):
                    ss.append(''.join(x.text or '' for x in e.iter() if x.tag.endswith('}t')));e.clear()
        loc=sheets(z)[sheet]
        for _,e in ET.iterparse(z.open(loc),events=['end']):
            if not e.tag.endswith('}row'): continue
            rn=int(e.attrib['r'])
            if wanted is None or rn in wanted:
                vals={}
                for c in e:
                    coord=c.attrib['r'];typ=c.attrib.get('t');v=c.find('s:v',NS)
                    if typ=='inlineStr': value=''.join(t.text or '' for t in c.iter() if t.tag.endswith('}t'))
                    elif v is None:value=None
                    elif typ=='s':value=ss[int(v.text)]
                    elif typ in ['str','e']:value=v.text
                    else:
                        value=float(v.text)
                        if value.is_integer():value=int(value)
                    vals[coord]=value
                yield rn,vals
            e.clear()
            if wanted is not None and rn>=max(wanted):break
def colnum(coord):
    n=0
    for c in coord:
        if not c.isalpha():break
        n=n*26+ord(c)-64
    return n

baseline=read(O/'BASELINE_HASHES_BEFORE.json')
source_files=[D/x for x in ['B02_groups_NODE.xlsx','B02_protein_NODE.xlsx','B02_phospho_NODE.xlsx','B02_TableS1.xlsx']]
source_now={str(f):sha(f) for f in source_files}
source_before=read(O/'SOURCE_HASHES_BEFORE.json') if (O/'SOURCE_HASHES_BEFORE.json').exists() else source_now
assert source_now==source_before
if not (O/'SOURCE_HASHES_BEFORE.json').exists():save('SOURCE_HASHES_BEFORE.json',source_before)

# A01: deterministic algebra and an explicit comparison policy; no replacement tests.
components=read(A/'H01_PATIENT_COMPONENTS.json')
values=[x for x in read(P/'results/H01_582/PATIENT_VALUES.json') if x['site']=='S283']
by_case={x['case']:x for x in values}
errors=[]
for x in components:
    assert x['case'] in by_case
    errors.extend([abs(x['source_adjusted']-x['net_subtracted_offset']-x['centered_adjusted']),
                   abs(x['source_adjusted']-by_case[x['case']]['unscaled_adjusted']),
                   abs(x['centered_adjusted']-by_case[x['case']]['adjusted'])])
assert len(components)==len(by_case)==159 and max(errors)<1e-12
summary={k:statistics.mean(x[k] for x in components) for k in ['source_adjusted','centered_adjusted','net_subtracted_offset']}
summary.update(n=159,sign_flips=sum((x['source_adjusted']>0)!=(x['centered_adjusted']>0) for x in components),max_algebra_error=max(errors))
norm=dict(issue='CGNL-A01',utc=UTC,status='partial_resolved_comparison_policy_fixed_source_lineage_open',summary=summary,
 estimands=[{'id':'E_SOURCE','formula':'(site_T-site_N)-(protein_T-protein_N)','meaning':'作者或PDC导出参照下的位点相对总蛋白变化；保留源处理差异。'},
 {'id':'E_CENTERED','formula':'E_SOURCE-[(median_site_T-median_site_N)-(median_protein_T-median_protein_N)]','meaning':'相对于各样本全局位点及蛋白中位数的相对变化；全局磷酸化变化也被移除。'}],
 policy={'baseline_unchanged':True,'both_estimands_reported':True,'select_by_P_or_direction':False,'pool_cross_cohort_effects':False,
 'rule':'跨队列先逐项确认log底、参照样本、源归一化、位点聚合和缺失处理。只有定义相同的估计量才比较效应；当前不合并效应、不以符号反转断言癌种特异性。H01不得作为方向稳定的外部确认。',
 'applies_to':['H01','B03','P01','H02','H03'],
 'known_result_exposure':'本规则在已查看既有结果后制定，是解释约束；不冒称前瞻性注册或新增独立验证。'},
 remaining=['H01历史CDAP导出文件对应的确切流水线版本/commit尚未绑定','与P01及其他源处理的严格等价性未确立'],
 acceptance='已完成算术、估计量区分和不依赖有利结果的报告规则；跨源测量等价性仍未通过。')
save('A01_NORMALIZATION_DECISION.json',norm)

# A04: independently reconstruct channel map from original design XLSX.
with zipfile.ZipFile(D/'B02_groups_NODE.xlsx') as z: design_sheets=list(sheets(z))
design_sheet=next(x for x in design_sheets if 'TMT' in x)
rr=list(rows(D/'B02_groups_NODE.xlsx',design_sheet))
channel_map=[]
for rn,cs in rr:
    batch=next((str(v) for v in cs.values() if str(v).startswith('iCCA_')),None)
    if batch is None:continue
    for coord,v in cs.items():
        if 3<=colnum(coord)<=12:
            assert isinstance(v,(int,float))
            channel_map.append(dict(case=str(int(v)),batch=batch,source_row=rn,source_column=colnum(coord)))
assert channel_map==read(A/'B02_TMT_CHANNEL_MAP.json')
mapping=collections.defaultdict(list)
for x in channel_map:mapping[x['case']].append(x)
cached={}
source_check=[]
for kind,name,sheet,target in [('protein','B02_protein_NODE.xlsx','proteomics_NA','Q0VF96'),('site','B02_phospho_NODE.xlsx','phospho_NA','CGNL1:S283')]:
    ds=next(x for x in read(P/('results/source_search_593/B02/'+name+'.extraction.json')) if x['sheet']==sheet)
    hit=next(x for x in ds['hits'] if x['values'][0]==target)
    observed=dict(rows(D/name,sheet,{1,hit['row']}))
    head=observed[1];vals=observed[hit['row']]
    cached[kind]={str(v):vals.get(''.join(c for c in coord if c.isalpha())+str(hit['row'])) for coord,v in head.items() if isinstance(v,int)}
    for coord,v in vals.items():
        cv=hit['values'][colnum(coord)-1]
        assert v==cv or (isinstance(v,(int,float)) and isinstance(cv,(int,float)) and abs(v-cv)<1e-12)
    with zipfile.ZipFile(D/name) as z: workbook_sheets=list(sheets(z))
    source_check.append(dict(file=name,sheets=workbook_sheets,target_row=hit['row'],source_target_matches_previous_extraction=True))
assert set(cached['site'])==set(cached['protein']) and len(cached['site'])==214
assert all(isinstance(v,(float,int)) and math.isfinite(v) for d in cached.values() for v in d.values())
patients=[]
for case in cached['site']:
    hist=mapping[case]
    patients.append(dict(case=case,channels=hist,unique_batch=hist[0]['batch'] if len(hist)==1 else None,
       batch_status='unambiguous' if len(hist)==1 else 'repeated_patient_aggregation_unresolved',
       S283=cached['site'][case],CGNL1_total_protein=cached['protein'][case],independent_patient_weight=1))
save('B02_PATIENT_BATCH_PROVENANCE.json',patients)
repeat={k:[x['batch'] for x in v] for k,v in mapping.items() if len(v)>1}
assert len(channel_map)==220 and len(mapping)==215 and len(repeat)==5
assert set(mapping)-set(cached['site'])=={'227'}
b02=dict(issue='CGNL-A04',utc=UTC,status='partial_resolved_repeated_design_verified_final_aggregation_unknown',
 total_channels=220,design_patients=215,final_patients=214,unambiguous_patients=209,repeated_patients=repeat,
 design_only_patient=['227'],source_checks=source_check,
 original_methods={'source':'sources/B02_Dong2022.pdf','doi':'10.1016/j.ccell.2021.12.006','first_author':'Dong',
 'pages':['PDF24/e5 clinical specimens','PDF25/e6 TMT design','PDF28/e9 normalization and batch QC'],
 'findings':['10 tumor channels and 131C pooled internal reference per plex','22 TMT plexes','both omics median normalized then log2 transformed','within-patient samples cluster adjacently in global-proteome QC'],
 'not_established':['precise mean/median/selection rule used to collapse repeated patients','omics-specific weights or retained replicate','reason for patient227 absence from final matrices']},
 rules=['保留214个独立患者，各计一次；220通道不能变为220例','5例重复患者不随意指定单一批次，不删除来制造合格主分析','209例是可定义的批次敏感性子集，不能默认为替换主队列','作者全局蛋白批次QC不能证明S283局部信号无批次影响'],
 missing_artifact='作者重复通道合并脚本/说明，或带通道标签的未合并定量矩阵及确定性的处理规则',
 model_freeze_allowed=False,clinical_inferential_tests_run=0)
save('A04_BATCH_DECISION.json',b02)

# A06: source-based evidence levels, explicit unresolved pairwise status.
h03=P/'metadata/systematic_scout_554/fulltexts/PMC10783603.xml'
root=ET.parse(h03).getroot();paras=[' '.join(''.join(x.itertext()).split()) for x in root.iter('p')]
h03_evidence=[x for x in paras if ('160 primary HCC patients' in x or 'To verify the robustness and universality' in x)]
assert len(h03_evidence)==2
save('H03_SOURCE_PARAGRAPHS.json',{'file':str(h03),'sha256':sha(h03),'paragraphs':h03_evidence})
overlap=dict(issue='CGNL-A06',utc=UTC,status='partial_resolved_source_distinction_supported_donor_linkage_unresolved',
 records=[{'cohort':'H01','publication':'Gao et al. Cell 2019','doi':'10.1016/j.cell.2019.08.052','source':'PDC000198/PDC000199','clinical_pairs_registered':160,'analysis_pairs':159,'recruitment_dates':None,'donor_crosswalk_to_other_HCC':None,'source_institution_evidence':'复旦/中山相关研究；本轮未取得可核对完整招募期的Methods，不从作者单位推断全部招募地点'},
 {'cohort':'H02','publication':'Jiang et al. Nature 2019','doi':'10.1038/s41586-019-0987-8','reporting_summary':'sources/H02_MOESM2.pdf','original_early_stage_cases':110,'analysis_eligible_pairs':53,'S283_complete_pairs':38,'recruitment_dates':None,'donor_crosswalk_to_other_HCC':None,'inclusion':'BCLC0/A，无既往化疗或放疗；报告表列108例HBV背景，不能将全部110例默认为HBV阳性'},
 {'cohort':'H03','publication':'Integrated omics landscape of hepatocellular carcinoma suggests proteomic subtypes for precision therapy','pmcid':'PMC10783603','recruitment_institution':'Mengchao Hepatobiliary Hospital of Fujian Medical University','original_patients':160,'S283_complete_pairs':118,'recruitment_dates':None,'donor_crosswalk_to_other_HCC':None,'independence_evidence':'本研究队列与Gao/Jiang外部数据在原文中分列；支持独立招募来源，但不等于患者级排重证明','measurement_limitation':'GRR插补及观测值掩码仍未解决'}],
 comparisons=[{'cohorts':['H01','H02'],'status':'unresolved','allowed_claim':'分别列示外部队列证据，不能认证无重复独立验证'},
 {'cohorts':['H01','H03'],'status':'distinct_recruitment_source_supported_patient_overlap_not_excluded','allowed_claim':'来源区分有支持，患者级独立性仍待验证'},
 {'cohorts':['H02','H03'],'status':'distinct_recruitment_source_supported_patient_overlap_not_excluded','allowed_claim':'来源区分有支持，患者级独立性仍待验证'},
 {'cohorts':['PDAC_phosphosite_357','PDC_processing_sensitivity_399'],'status':'same_44_discovery_patients','independent_cohort_count':1}],
 rule='不相加HCC样本量作为已去重总数；原ID不同不证明无重复，原ID相同也不证明重复。无日期/稳定ID/作者声明时保留unresolved。',
 next_evidence=['原始招募医院及时间段','可追溯去标识化跨研究患者映射或作者明确重叠声明','如有合法现成患者基因型，使用预定身份一致性规则而非表达相似度猜测'])
save('A06_COHORT_OVERLAP_LEDGER.json',overlap)

candidate=dict(utc=UTC,status='draft_not_frozen_not_fitted',whole_project_complete=False,
 population='B02 source-observed S283 and total-CGNL1 tumor samples,214 unique patients',
 candidate_question='S283信号在控制总CGNL1蛋白和预设临床因素后是否与血管侵犯相关？',
 endpoint='Vascular invasion,Yes=1,No=0',candidate_covariates=['S283 continuous','total CGNL1 protein continuous','Age continuous','Sex categorical','tumor size continuous with prespecified scale'],
 no_cutpoint_search=True,primary_target='S283 coefficient',no_outcome_P_values_examined=True,
 statistical_policy='新有效检验按目标系数原始P<0.05识别名义信号；多重家族预先定义并报告FDR。单个预设主要终点不强制FDR。',
 pending=['CGNL-A04重复定量合并来源','T04明确的问题与疾病模型依据','肿瘤大小函数形式、单位及编码正式固定','有依据的批次敏感性方案及模型容量验收'],
 freeze_allowed=False)
save('B02_MODEL_READINESS.json',candidate)
source_status=[]
for f in (O/'sources').glob('*.pdf'):
    valid=f.read_bytes()[:5]==b'%PDF-'
    source_status.append(dict(file=f.name,sha256=sha(f),pdf_magic=valid,usable_as_PDF=valid,
       note='已提取文本' if valid and f.with_suffix('.txt').exists() else 'HTTP200不等于PDF有效；未作为全文证据'))
save('SOURCE_CONTENT_VALIDATION.json',source_status)
save('REQUIRED_SOURCE_EVIDENCE.json',{'B02':['647/321/133/587/823在原批次与iCCA_22之间的最终选值或聚合规则，蛋白及磷酸化分别记录','227被排除的原始理由','对应代码版本或通道级定量文件'],
 'HCC':['H01/H02/H03实际招募医院和起止时间','跨研究患者映射或明确不重叠声明','H03观测值掩码与插补位置记录'],
 'H01':['对应已下载CDAP导出矩阵的历史版本/参数文件'],
 'contacts_sent':False,'scope':'所需证据清单，不代表已向作者发信或获得材料'})
issues=[{'id':'CGNL-A01','status':'部分解决','completed':'159例代数复核与双估计量比较规则落地','remaining':'历史CDAP版本和跨源测量等价性'},
 {'id':'CGNL-A04','status':'部分解决','completed':'原始工作簿复核220通道/215设计ID/214最终患者；5例重复映射','remaining':'作者最终合并/选择规则及227缺失原因'},
 {'id':'CGNL-A06','status':'部分解决','completed':'H03来源证据、H02报告表入选信息及逐对重叠台账','remaining':'患者级跨研究排重或明确无重叠证据'},
 {'id':'CGNL-A08','status':'已解决于新台账','severity':'low','evidence':'B02 DOI与原文第一作者为Dong，旧清单记为Zou；本轮新台账使用Dong，旧清单保留','impact':'书目定位准确性，统计结果不受影响'}]
save('ISSUE_FOLLOWUP_STATUS.json',dict(utc=UTC,issues=issues,untouched_prior_issues=['CGNL-A02','CGNL-A03','CGNL-A05','CGNL-A07']))

report='''CGNL项目续推进核查记录（2026-10-23）

本轮落实CGNL-A01、A04、A06的来源与方法核查，保留此前分析、审核和总计划。未新增临床推断模型，未把数值QA当作机制验证。

CGNL-A01：部分解决。
159例H01的源口径均值为+0.3498475432，额外中心化净偏移为0.4194519818，得到-0.0696044386；62例符号改变。逐患者代数及原分析结果对应通过。已固定双口径并列规则，不按P值或有利方向挑选归一化；当前不合并跨源效应。额外中心化会移除全局位点水平变化，两个口径不能解释为同一生物量。历史CDAP版本和严格跨源可比性仍未确认。

CGNL-A04：部分解决，B02模型尚不能冻结。
已从原始XLSX独立复核220个肿瘤通道、215个设计ID及214个最终患者。647、321、133、587、823各跨两个批次；227存在于设计但不在最终矩阵。S283与总CGNL1源表目标行复核一致、214例皆有数值。原文确认22个TMT批次、131C混合内参、两组学中位数归一化及重复患者QC，但未找到最终重复值的合并或选择规则。保留214例一人一份；不随意赋单批次或删除5例。209例只能作为明确界定的敏感性子集。需要作者处理脚本/说明，或带通道标签的定量矩阵与确定性重建规则。

CGNL-A06：部分解决。
H03原文明确来自福建孟超医院，并将Gao/Jiang数据作为外部队列处理。来源区分得到支持，但患者级无重叠尚未认证。H02公开报告表确认110例早期HCC，其中108例有HBV背景；本项目53对合格、38对S283完整值的分析范围保持原样。H01/H02/H03的招募时间和可跨研究链接的患者标识仍缺失，逐对保留未解决状态，不相加为已去重独立样本总量。PDAC的替代处理仍只计同一44例发现队列一次。

有依据的新增问题CGNL-A08：旧清单B02作者简写为Zou，原文第一作者为Dong；已在新台账补正，原清单未改，数值结果不变。

B02建模：已保存未冻结的准备记录，主要问题限定为控制总CGNL1蛋白与预设临床因素后的S283-血管侵犯关联。尚缺重复通道处理来源、T04疾病问题依据、正式变量尺度及批次敏感性规格。本轮不拟合模型、不查看该关联P值。

来源：B02原文DOI 10.1016/j.ccell.2021.12.006（UEA公开全文，PDF第24、25、28页）；H02 DOI 10.1038/s41586-019-0987-8的公开Reporting Summary；H03 PMC10783603原文；项目原始工作簿、H01患者级结果及首轮审核。

研究状态：仍为候选线索，whole_project_complete=false。上述未决项是具体证据缺口，不以继续等待下载或监控进程存活替代解决。
'''
(O/'CGNL1_FOLLOWUP_REPORT_20261023.txt').write_text(report,encoding='utf-8')
after={f:sha(Path(f)) for f in baseline}
source_after={str(f):sha(f) for f in source_files}
assert baseline==after and source_before==source_after
save('READ_ONLY_ACCEPTANCE.json',dict(utc=UTC,baseline_files_unchanged=len(baseline),source_workbooks_unchanged=len(source_before),
 H01_patient_algebra_pass=True,B02_original_design_matches=True,B02_target_rows_reverified=True,independent_patients=214,
 new_clinical_models=0,new_inferential_tests=0,scientific_mainline_approved=False))
artifacts=[dict(path=str(f),bytes=f.stat().st_size,sha256=sha(f)) for f in O.rglob('*') if f.is_file() and f.name!='FOLLOWUP_ACCEPTANCE.json']
save('FOLLOWUP_ACCEPTANCE.json',dict(utc=UTC,authorized_followup_executed=True,all_issues_resolved=False,
 model_freeze_allowed=False,whole_project_complete=False,artifacts=artifacts))
print(json.dumps(dict(output=str(O),normalization=summary,repeat_patients=repeat,baseline_files_unchanged=len(baseline),artifact_count=len(artifacts)),ensure_ascii=False))
