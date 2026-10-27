from pathlib import Path
import zipfile,hashlib,json,io,re,shutil,datetime
import pandas as pd
ROOT=Path(r'C:\Users\31937\CGNL1_identity_review_20261002')
OUT=ROOT/'external_review_context_exposed_20261002_v1';RAW=OUT/'frozen_original';RAW.mkdir(exist_ok=True)
source=Path(r'C:\Users\31937\Downloads\匿名身份审阅_冻结结果_v1.zip');original_packet=Path(r'C:\Users\31937\Downloads\单细胞身份_匿名盲审_20261002_v1.zip')
with zipfile.ZipFile(source) as z:
 assert z.testzip() is None
 freeze=json.loads(z.read('SHA256_FREEZE.json'))
 assert hashlib.sha256(original_packet.read_bytes()).hexdigest()==freeze['input_sha256']
 allowed={'SHA256_FREEZE.json'}|{f['name'] for f in freeze['files']};assert set(z.namelist())==allowed
 for name in z.namelist():
  assert Path(name).name==name and '/' not in name and '\\' not in name
  blob=z.read(name)
  if name!='SHA256_FREEZE.json':
   rec=next(x for x in freeze['files'] if x['name']==name);assert len(blob)==rec['bytes'] and hashlib.sha256(blob).hexdigest()==rec['sha256']
  (RAW/name).write_bytes(blob)
 tname=next(n for n in z.namelist() if n.endswith('.tsv'));t=pd.read_csv(io.BytesIO(z.read(tname)),sep='\t',dtype=str,keep_default_na=False)
shutil.copy2(source,RAW/source.name)
public=Path(r'C:\Users\31937\Downloads\单细胞身份_匿名盲审_20261002_v1');template=pd.read_csv(public/'05_审阅结果模板.tsv',sep='\t',dtype=str,keep_default_na=False)
assert len(t)==44 and t[['dataset','cluster']].equals(template[['dataset','cluster']]);assert t[['dataset','cluster']].drop_duplicates().shape[0]==44;assert (t!='').all().all()
key=pd.read_csv(ROOT/'blind_handoff_private_20261002/解盲对照_禁止上传.tsv',sep='\t',dtype=str)
m=t.merge(key,left_on=['dataset','cluster'],right_on=['anonymous_dataset','anonymous_cluster'],validate='one_to_one')
immune={'myeloid','T_NK','B','plasma','erythroid'}
def refclass(s):
 if s in immune:return 'immune'
 if s in ['fibroblast','mural']:return 'stroma'
 return s

def reviewer_class(s):
 if s.startswith('谱系未定'):return 'unresolved'
 if '内皮' in s:return 'endothelial'
 if '成纤维' in s:return 'stroma'
 if '肝细胞' in s:return 'hepatocyte'
 if '上皮' in s:return 'epithelial'
 if any(x in s for x in ['T细胞','B细胞','淋巴','髓','单核','巨噬','树突','肥大','中性粒','浆细胞']):return 'immune'
 raise ValueError(s)
m['review_coarse_class']=m.broad_identity.map(reviewer_class);m['disease_reference_coarse']=m.current_disease_reference.map(refclass);m['normal_reference_coarse']=m.current_normal_reference.map(refclass)
def status(row):
 if row.review_coarse_class=='unresolved':return '审阅者未定谱系，不能视为确认自动标签'
 a=row.review_coarse_class==row.disease_reference_coarse;b=row.review_coarse_class==row.normal_reference_coarse
 return '粗谱系与双参考一致' if a and b else ('仅与肿瘤参考一致' if a else ('仅与正常参考一致' if b else '与双参考均有分歧'))
m['coarse_comparison']=m.apply(status,axis=1)
notes={'C006':'自动T/NK标签过于确定；增殖状态不能替代谱系，保留未定。','C011':'pDC类别未在现有粗参考标签中单列，不能将B标签当作最终身份。','C012':'肥大细胞类别缺失，自动T/NK标签需复核。','C013':'关键跨谱系分歧：审阅为上皮，双参考均预测T/NK；不能因双参考一致即认定其正确。','C032':'髓系大类一致，单核/巨噬/DC边界尚未解决。','C039':'支持细胞毒淋巴群，NK与T边界未确定。','C026':'髓系大类一致，需单核/巨噬/DC细分证据。','C038':'上皮来源与完整细胞资格分别判断；低RNA条形码需补QC。','C041':'肥大细胞类别缺失，原参考预测不应当作最终细分身份。','C043':'基质谱系基本相符；成纤维/壁细胞边界需补标记及共表达。'}
m['comparison_note']=m.cluster.map(notes).fillna('保留原冻结判断；粗谱系一致不代表细亚型已确认。')
m['review_blindness_status']='匿名审阅；系统历史上下文暴露已披露；非严格独立盲审'
m.to_csv(OUT/'44簇解盲比较.tsv',sep='\t',index=False)
priority=m[(m.confidence!='高')|m.cluster.isin(notes)|m.coarse_comparison.isin(['与双参考均有分歧','仅与正常参考一致'])].copy();priority.to_csv(OUT/'主要分歧与待补证据.tsv',sep='\t',index=False)
# Verify explicitly quoted support genes and their supplied rounded values.
markers=pd.read_csv(public/'02_top30_markers.tsv',sep='\t');lookup=markers.set_index(['dataset','cluster','gene']);gene_checks=[];value_checks=[]
for row in t.itertuples():
 for gene in re.findall(r'([A-Z][A-Z0-9-]*)[（(]',row.supporting_genes):gene_checks.append((row.dataset,row.cluster,gene,(row.dataset,row.cluster,gene) in lookup.index))
 for gene,pct,delta in re.findall(r'([A-Z][A-Z0-9-]*)[（(]检出\s*([0-9.]+)%[，,]Δ\s*([-0-9.]+)[）)]',row.supporting_genes):
  k=(row.dataset,row.cluster,gene)
  if k in lookup.index:
   v=lookup.loc[k];value_checks.append(dict(dataset=row.dataset,cluster=row.cluster,gene=gene,pct_delta=abs(float(pct)-float(v.fraction_detected)*100),mean_delta=abs(float(delta)-float(v.mean_difference_vs_other_clusters))))
assert all(x[3] for x in gene_checks)
errors=[x for x in value_checks if x['pct_delta']>.051 or x['mean_delta']>.006]
pd.DataFrame(value_checks).to_csv(OUT/'标记数值抽取核对.tsv',sep='\t',index=False)
mito={'dataset_anonymous':'D02','dataset_original':'GSE245552','checked_feature_files':17,'features_per_file':29569,'MT_dash_matching_genes_per_file':0,'compact_canonical_MT_names_per_file':0,'interpretation':'Existing mitochondrial QC values are not evaluable under the current feature-name matching rule. The reported zero is not evidence of biologically low mitochondrial proportion. Author preprocessing or missing feature mapping requires source review.','action':'Record QC caveat; preserve frozen reviewer files and all current analysis outputs unchanged.'}
(OUT/'MITO_QC_AUDIT.json').write_text(json.dumps(mito,ensure_ascii=False,indent=2),encoding='utf-8')
receipt={'imported_at':datetime.datetime.now().astimezone().isoformat(),'input_zip':str(source),'input_zip_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'original_evidence_zip_hash_matches':True,'frozen_files_hashes_match':True,'rows':len(t),'columns':len(t.columns),'template_order_matches':True,'confidence_counts':t.confidence.value_counts().to_dict(),'blindness':'context_exposure_disclosed_not_strict_independent_blind','reviewer_statement_preserved':freeze['blindness'],'coarse_comparison_counts':m.coarse_comparison.value_counts().to_dict(),'support_gene_checks':len(gene_checks),'numeric_support_checks':len(value_checks),'numeric_check_errors':errors,'endothelial_clusters':m[m.review_coarse_class=='endothelial'][['dataset','cluster','original_dataset','original_cluster','confidence']].to_dict('records'),'labels_applied_to_analysis':False}
(OUT/'IMPORT_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
report='''网页匿名审阅解盲比较（冻结原件保留）

1. 输入证据包SHA256与网页冻结清单一致；两份审阅文件哈希和长度全部匹配。44行、10列、匿名编号及模板顺序核对通过。
2. 本版本的正确定位：匿名证据审阅，历史上下文暴露已披露。审阅者明确说明系统附带过历史项目摘要，因此不是严格独立盲审。原声明、报告和TSV保持原样。
3. 两个内皮簇的主要身份得到本次审阅支持：D01/C018→GSE285991原簇14；D02/C028→GSE245552原簇21；均为高置信度。此处是簇级身份支持，不证明簇内每个条形码都纯净。
4. 不能将双参考一致视为最终答案。D01/C013被本次审阅判为上皮，而原双参考均预测T/NK，是需重点解决的跨谱系分歧。D01/C011的pDC、D01/C012与D02/C041的肥大细胞，提示原参考标签范围不足。D01/C006应保留增殖、谱系未定。NK/T、单核/巨噬/DC、成纤维/壁细胞及低RNA上皮条形码仍需补证据。
5. 已只读核查线粒体QC疑点：D02对应队列17份特征表均无当前规则匹配的线粒体基因名，原零值不能解释为线粒体比例低。应按此项不可评估记录，并核查作者预处理和特征映射。没有据此改动原矩阵、当前簇或冻结审阅答案。
6. 建议由未接触项目资料的另一环境，对同一匿名证据独立审阅；不要提供本文件、原参考标签或本次审阅结果。可附单独生成的匿名QC说明，防止把D02零值误读为实测低比例。等第二份答案冻结后再三方比较。
7. 当前仅完成接收、校验与比较；没有把网页标签直接回填为最终分析注释，也没有重跑目标关联。
'''
(OUT/'解盲比较说明.txt').write_text(report,encoding='utf-8')
qc_note='''匿名证据补充说明：线粒体QC字段
日期：2026-10-02

本说明不提供原始数据集、项目身份、细胞类型标签、自动预测或其他审阅结论。

对匿名队列D02的输入特征表核查发现：目前计算规则无法匹配线粒体基因特征。因此，原概览表和逐簇证据中D02的“线粒体比例中位数=0”不能解释为实测比例为零，也不能据此认定其质量更好。

本轮身份审阅请将D02这一字段一律按“不可评估/NA”处理；它不构成低线粒体比例或该QC步骤有效的证据。其他标记、表达差异、检出比例、供者计数、匿名编号及UMAP保持原样。不据这一发现自动删除任何簇，也不把它当作支持某种细胞身份的证据。

请在审阅结果中注明：收到这份匿名QC附注。完成并冻结全部簇的独立判断后，再接收任何解盲信息。
'''
qpath=Path(r'C:\Users\31937\Downloads\匿名QC补充说明_D02_20261002.txt');qpath.write_text(qc_note,encoding='utf-8')
assert not any(x in qc_note for x in ['CGNL1','GSE245552','GSE285991','内皮','免疫','C018','C028'])
print(json.dumps(receipt,ensure_ascii=False,indent=2));print('QC_NOTE',qpath)