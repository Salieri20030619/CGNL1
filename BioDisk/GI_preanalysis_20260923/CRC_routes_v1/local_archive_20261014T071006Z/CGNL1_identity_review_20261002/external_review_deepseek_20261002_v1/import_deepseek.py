from pathlib import Path
import re,json,hashlib,shutil,datetime
import pandas as pd
ROOT=Path(r'C:\Users\31937\CGNL1_identity_review_20261002');OUT=ROOT/'external_review_deepseek_20261002_v1';RAW=OUT/'frozen_original';RAW.mkdir(exist_ok=True)
source=Path(r'C:\Users\31937\Downloads\deepseek_markdown_20261002_56c7cf.md');blob=source.read_bytes();text=blob.decode('utf-8-sig');(RAW/source.name).write_bytes(blob)
columns=['dataset','cluster','broad_identity','optional_subtype','supporting_genes','conflicting_evidence','confidence','mixed_or_doublet_concern','additional_evidence_needed','reasoning']
rows=[]
for line in text.splitlines():
 if re.match(r'^\|\s*D0[12]\s*\|',line):
  row=[x.strip() for x in line.strip().strip('|').split('|')];assert len(row)==10,(len(row),row[:2]);rows.append(row)
t=pd.DataFrame(rows,columns=columns);assert len(t)==44 and len(t[['dataset','cluster']].drop_duplicates())==44 and (t!='').all().all()
public=Path(r'C:\Users\31937\Downloads\单细胞身份_匿名盲审_20261002_v1');template=pd.read_csv(public/'05_审阅结果模板.tsv',sep='\t',dtype=str,keep_default_na=False);assert t[['dataset','cluster']].equals(template[['dataset','cluster']]);t.to_csv(OUT/'DeepSeek_44簇结构化原判.tsv',sep='\t',index=False)
web=pd.read_csv(ROOT/'external_review_context_exposed_20261002_v1/frozen_original/匿名身份审阅_44簇_v1.tsv',sep='\t',dtype=str,keep_default_na=False)
key=pd.read_csv(ROOT/'blind_handoff_private_20261002/解盲对照_禁止上传.tsv',sep='\t',dtype=str)
m=t.merge(web,on=['dataset','cluster'],suffixes=('_deepseek','_web'),validate='one_to_one').merge(key,left_on=['dataset','cluster'],right_on=['anonymous_dataset','anonymous_cluster'],validate='one_to_one')
# Explicitly expand compact gene notation, not assume it was an observed single symbol.
markers=pd.read_csv(public/'02_top30_markers.tsv',sep='\t');gene_lookup={(ds,cl):set(g.gene) for (ds,cl),g in markers.groupby(['dataset','cluster'])}
def expand(token):
 parts=token.split('/');base=parts[0];out=[base]
 for s in parts[1:]:
  if s.isdigit():out.append(re.sub(r'\d+$','',base)+s)
  elif re.fullmatch(r'[A-Z]+\d+',s):out.append(re.sub(r'[A-Z]+\d+$','',base)+s)
  elif re.fullmatch(r'[A-Z]+',s):out.append(re.sub(r'[A-Z]+$','',base)+s)
  else:out.append(s)
 return out
checks=[]
for row in t.itertuples():
 supplied=gene_lookup[(row.dataset,row.cluster)]
 for token in re.split(r'[,，、]\s*',row.supporting_genes):
  token=token.strip();uncertain='?' in token or '？' in token;token=token.replace('?','').replace('？','').strip()
  if not re.fullmatch(r'[A-Z][A-Z0-9/-]*',token):
   checks.append(dict(dataset=row.dataset,cluster=row.cluster,original_token=token,gene='',status='requires_manual_token_review'));continue
  for gene in expand(token):checks.append(dict(dataset=row.dataset,cluster=row.cluster,original_token=token,gene=gene,status='explicitly_uncertain_not_observed_evidence' if uncertain else ('in_supplied_top30' if gene in supplied else 'not_in_supplied_top30')))
gc=pd.DataFrame(checks);gc.to_csv(OUT/'支持基因来源核对.tsv',sep='\t',index=False);unprovided=gc[gc.status=='not_in_supplied_top30'];unprovided.to_csv(OUT/'支持列表中未由材料提供的基因.tsv',sep='\t',index=False)
def family(s,sub):
 s=s.replace(' ','')
 if s.startswith('谱系未定') or s=='增殖细胞':return 'lineage_unresolved_cycling'
 if '内皮' in s:return 'endothelial'
 if '成纤维' in s or '基质细胞' in s:return 'fibroblast_stromal'
 if '肝细胞' in s:return 'hepatocyte'
 if '上皮' in s:return 'epithelial'
 if '浆细胞样树突' in s or ('树突' in s and 'pDC' in sub):return 'pDC'
 if '肥大' in s:return 'mast'
 if '中性粒' in s:return 'neutrophil'
 if '单核' in s or '巨噬' in s or '髓系' in s:return 'mononuclear_phagocyte'
 if '浆细胞' in s:return 'plasma'
 if 'NK' in s or '细胞毒淋巴' in s:return 'cytotoxic_lymphocyte_NK_tendency'
 if 'T细胞' in s:return 'T_cell'
 if 'B细胞' in s:return 'B_cell'
 raise ValueError((s,sub))
m['harmonized_deepseek']=[family(a,b) for a,b in zip(m.broad_identity_deepseek,m.optional_subtype_deepseek)];m['harmonized_web']=[family(a,b) for a,b in zip(m.broad_identity_web,m.optional_subtype_web)];m['broad_agreement']=m.harmonized_deepseek==m.harmonized_web
m['deepseek_unsupported_support_genes']=m.apply(lambda r:';'.join(unprovided[(unprovided.dataset==r.dataset)&(unprovided.cluster==r.cluster)].gene),axis=1)
notes={'C006':'两份审阅均缺乏确定谱系证据；增殖是状态，保留谱系未定。','C003':'主要身份标签不一致：DeepSeek写增殖细胞，网页写T细胞；所给CD3D/E检出支持T谱系，增殖应另列为状态。','C039':'可保留NK倾向；未进前30的CD3不能当阴性，不能据此确定纯NK。','C035':'T身份一致；Treg细分仍为倾向，FOXP3?不能记作实测支持。','C026':'髓系大类一致；单核/巨噬/DC边界仍需补证据，不能因模型高置信度升级。','C038':'上皮来源得到支持；完整细胞、低RNA及背景问题另行核查。','C042':'粒细胞可天然低RNA，低UMI不等于质量差或双细胞；保留中性粒身份与QC问题分开。','C013':'两份审阅均判上皮，原双参考均预测T/NK；需修正对自动参考一致性的信赖，保留分歧证据。','C011':'两份审阅支持pDC，原双参考B标签没有覆盖该类别，不能原样继承。','C012':'两份审阅支持肥大细胞，原双参考T/NK不应作为最终身份。','C041':'两份审阅支持肥大细胞，原参考标签粒度不足。','C043':'成纤维/基质一致；壁细胞/肌成纤维细分仍需完整标记。'}
m['adjudication_note']=m.cluster.map(notes).fillna('主要谱系一致时可形成簇级候选共识；亚型及单细胞纯度不自动升级。')
m.to_csv(OUT/'44簇三方解盲比较.tsv',sep='\t',index=False)
# Provisional consensus is a new table, not an overwrite of either reviewer or analysis annotations.
consensus=[]
name_map={'endothelial':'血管内皮细胞','fibroblast_stromal':'成纤维/基质细胞','hepatocyte':'肝细胞','epithelial':'上皮细胞','pDC':'浆细胞样树突细胞','mast':'肥大细胞','neutrophil':'中性粒细胞','mononuclear_phagocyte':'单核吞噬细胞（单核/巨噬边界待定）','plasma':'浆细胞谱系','cytotoxic_lymphocyte_NK_tendency':'细胞毒淋巴细胞（NK倾向）','T_cell':'T细胞','B_cell':'B细胞','lineage_unresolved_cycling':'谱系未定；增殖状态'}
for row in m.itertuples():
 proposed=name_map[row.harmonized_web] if row.broad_agreement else 'T细胞；增殖状态（基于已提供CD3D/E证据的待核定建议）'
 if row.cluster=='C038':proposed='上皮来源细胞/条形码；完整细胞资格待核查'
 consensus.append(dict(dataset=row.dataset,cluster=row.cluster,original_dataset=row.original_dataset,original_cluster=row.original_cluster,proposed_broad_identity=proposed,reviewer_broad_agreement=bool(row.broad_agreement),decision_basis=row.adjudication_note,status='簇级建议；未自动覆盖细胞级注释；目标关联未重测'))
pd.DataFrame(consensus).to_csv(OUT/'簇级共识建议_未回填.tsv',sep='\t',index=False)
receipt={'imported_at':datetime.datetime.now().astimezone().isoformat(),'source_file':str(source),'source_sha256':hashlib.sha256(blob).hexdigest(),'rows':len(t),'columns':10,'all_cluster_ids_and_order_match':True,'reviewer_blindness_statement':'declares no prior labels, automatic predictions or prior conclusions','independence_verification':'reviewer declaration only; source is Markdown without original input hash or external freeze manifest','QC_addendum_acknowledgement':'not explicit in the supplied Markdown','confidence_counts_original':t.confidence.value_counts().to_dict(),'harmonized_major_identity_agreement':int(m.broad_agreement.sum()),'harmonized_comparison_n':44,'note_on_agreement':'Allows synonyms and coarse resolution; does not imply matching confidence or subtype certainty. C006 remains unresolved; C003 differs in identity/state labeling.','support_token_checks':len(gc),'support_genes_not_in_supplied_top30':len(unprovided),'clusters_with_unprovided_support_genes':len(unprovided[['dataset','cluster']].drop_duplicates()),'uncertain_support_tokens':int((gc.status=='explicitly_uncertain_not_observed_evidence').sum()),'labels_applied_to_analysis':False,'target_association_tests_run':False}
(OUT/'IMPORT_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(receipt,ensure_ascii=False,indent=2));print('UNPROVIDED_SUPPORT',unprovided.to_string(index=False));print('CORE_ENDO',m[m.harmonized_web=='endothelial'][['dataset','cluster','confidence_deepseek','confidence_web','deepseek_unsupported_support_genes']].to_string(index=False));print('DISAGREEMENTS',m[~m.broad_agreement][['dataset','cluster','broad_identity_deepseek','broad_identity_web']].to_string(index=False))