from pathlib import Path
import json,gzip,re,datetime
import pandas as pd,numpy as np
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/reassessment_20261008_v1');D=R/'perturbation_inputs';O=R/'results';spec=json.loads((R/'REASSESSMENT_SPEC.json').read_text());genes=['CGNL1']+spec['R4_real_perturbation']['positive_controls'];rows=[];samples=[];qa=[]
with gzip.open(D/'GSE163568_family.soft.gz','rt') as f:soft=f.read()
records=[];cur=None
for line in soft.splitlines():
 if line.startswith('^SAMPLE = '):
  cur={'accession':line.split(' = ',1)[1]};records.append(cur)
 elif cur is not None and line.startswith('!Sample_') and ' = ' in line:
  k,v=line.split(' = ',1);cur.setdefault(k,[]).append(v)
for cell,rawfile in [('HUVEC','GSE163568_TLA_HUVEC_SE_counts.txt.gz'),('HREC','GSE163568_TLA-HRECcounts.txt.gz')]:
 raw=pd.read_csv(D/rawfile,sep=r'\s+',index_col=0);assert raw.index.is_unique
 cols=['Fc1','Fc2','Fc3','DLL4_1','DLL4_2','DLL4_3'];assert raw.columns.tolist()==cols;assert np.isfinite(raw.to_numpy()).all() and (raw.to_numpy()>=0).all() and np.equal(raw.to_numpy(),np.round(raw.to_numpy())).all()
 de=pd.read_excel(D/f'GSE163568_DEG_{cell}_TLA.xlsx',engine='openpyxl')
 # Author spreadsheets converted some MARCH/SEPT-like symbols to dates. Quarantine,
 # never guess replacements. CGNL1 and the fixed controls must retain unique mapping.
 bad=~de.Gene.map(lambda x:isinstance(x,str))|de.Gene.duplicated(keep=False)
 issues=de[bad].copy();issues.to_csv(O/f'R4_{cell}_author_gene_annotation_issues.tsv',sep='\t',index=False)
 assert not de.loc[bad,'Gene'].isin(genes).any()
 de=de[~bad].set_index('Gene');assert de.index.is_unique
 mapping={}
 for rec in records:
  title=rec.get('!Sample_title',[''])[0]
  if title.startswith(cell+'_TLA_'):
   replicate=int(re.search(r'rep(\d+)',title).group(1));column=('Fc' if 'Control' in title else 'DLL4_')+str(replicate);mapping[column]=rec
 assert set(mapping)==set(cols)
 common=raw.index.intersection(de.index);factor={};maxcv=0
 for col in cols:
  rr=raw.loc[common,col].to_numpy(float);nn=de.loc[common,col].to_numpy(float);good=(rr>0)&(nn>0)&np.isfinite(nn);rat=rr[good]/nn[good];med=float(np.median(rat));factor[col]=med;maxcv=max(maxcv,float(np.quantile(abs(rat/med-1),.999)))
 assert maxcv<1e-6,('raw_author_normalized_alignment_failed',cell,maxcv)
 for gene in genes:
  if gene not in raw.index or gene not in de.index:
   rows.append(dict(dataset='GSE163568',cell_type=cell,gene=gene,status='not_measured'));continue
  r=de.loc[gene];v=r[cols].to_numpy(float);fc=float(np.log2((v[3:].mean()+1)/(v[:3].mean()+1)))
  rows.append(dict(dataset='GSE163568',cell_type=cell,gene=gene,status='qualified_exploratory_assay_response',treatment='tethered_DLL4_6h',n_control_assay_replicates=3,n_treated_assay_replicates=3,donor_material='two_donor_pools' if cell=='HUVEC' else 'single_donor',control_mean_normalized=float(v[:3].mean()),treated_mean_normalized=float(v[3:].mean()),log2_ratio_means_plus1=fc,author_log2FoldChange=float(r.log2FoldChange),author_P=float(r.pvalue),author_FDR=float(r.padj),treated_above_control_mean=int((v[3:]>v[:3].mean()).sum()),all_treated_above_all_controls=bool(v[3:].min()>v[:3].max()),inference_scope='author_DESeq2_assay_level_not_independent_patient_or_donor_inference'))
  for col,val in zip(cols,v):samples.append(dict(dataset='GSE163568',cell_type=cell,gene=gene,column=col,GSM=mapping[col]['accession'],title=mapping[col]['!Sample_title'][0],condition='control' if col.startswith('Fc') else 'DLL4',raw_count=float(raw.loc[gene,col]),author_normalized_count=float(val)))
 qa.append(dict(cell_type=cell,features=len(raw),sample_mapping_verified=True,raw_integer_counts=True,quarantined_author_gene_annotation_rows=len(issues),target_controls_mapping_unaffected=True,author_normalization_constant_factors=factor,normalization_ratio_relative_error_p999=maxcv,reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC8804000/',donor_independence_limit='triplicate assays from two HUVEC pools or one HREC donor; no new donor-level inference'))
pd.DataFrame(rows).to_csv(O/'R4_DLL4_target_and_controls.tsv',sep='\t',index=False);pd.DataFrame(samples).to_csv(O/'R4_DLL4_individual_values.tsv',sep='\t',index=False);(O/'R4_DLL4_QUALIFICATION.json').write_text(json.dumps(qa,indent=2));(O/'DLL4_COMPLETE.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'cell_systems':2,'raw_counts_and_author_normalization_verified':True,'new_inferential_tests':False,'author_statistics_preserved':True},indent=2));print(pd.DataFrame(rows).query("gene=='CGNL1'").to_json(orient='records'))
