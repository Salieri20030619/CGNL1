from common import *
import gzip,re
p=R/'raw_small/GSE228955_Normalized_counts_cpm.txt.gz'
t=pd.read_csv(p,sep='\t',index_col=0);t.index=t.index.astype(str)
mapping={}
for col in t.columns:
 if col.startswith('CGNL1_KO_22'):mapping[col]='CGNL1_KO_22'
 elif col.startswith('CGN_CGNL1_dKO_21D3'):mapping[col]='double_KO_21D3'
 elif col.startswith('CGN_KO_1G8'):mapping[col]='CGN_KO_1G8'
 elif col.startswith('WT_'):mapping[col]='WT'
assert len(mapping)==12 and all(list(mapping.values()).count(x)==3 for x in set(mapping.values()))
metadata=[];current={}
with gzip.open(R/'raw_small/GSE228955.soft.gz','rt') as f:
 for line in f:
  if line.startswith('^SAMPLE'):
   if current:metadata.append(current)
   current={'accession':line.split(' = ',1)[1].strip(),'characteristics':[]}
  elif line.startswith('!Sample_title'):current['title']=line.split(' = ',1)[1].strip()
  elif line.startswith('!Sample_characteristics_ch1'):current['characteristics'].append(line.split(' = ',1)[1].strip())
 if current:metadata.append(current)
(R/'metadata/GSE228955_samples.json').write_text(json.dumps(metadata,indent=2))
wanted=S['genes']+['CGNL1','CGN','ARHGEF2','TIAM1','RHOA','RAC1','TJP1','PLEKHA7']
valid=[];issues=[]
for g in wanted:
 match=np.flatnonzero(t.index.str.upper()==g)
 if len(match)==1:valid.append((g,int(match[0])))
 else:issues.append({'gene':g,'matching_rows':len(match)})
rows=[]
for g,i in valid:
 for sample,group in mapping.items():
  value=float(t.iloc[i][sample]);assert np.isfinite(value) and value>=0
  rows.append(dict(gene=g,sample=sample,group=group,CPM=value,ln1pCP10k=float(np.log1p(value/100))))
v=pd.DataFrame(rows);v.to_csv(O/'D5_KO_individual_gene_values.tsv',sep='\t',index=False)
summ=[]
for g,q in v.groupby('gene'):
 wt=q[q.group.eq('WT')]
 for group in ['CGNL1_KO_22','CGN_KO_1G8','double_KO_21D3']:
  ko=q[q.group.eq(group)]
  summ.append(dict(gene=g,group=group,WT_mean_CPM=wt.CPM.mean(),KO_mean_CPM=ko.CPM.mean(),log2FC_CPMplus1=float(np.log2((ko.CPM.mean()+1)/(wt.CPM.mean()+1))),delta_ln1pCP10k=float(ko.ln1pCP10k.mean()-wt.ln1pCP10k.mean()),all_KO_above_WT=bool(ko.CPM.min()>wt.CPM.max()),all_KO_below_WT=bool(ko.CPM.max()<wt.CPM.min()),independent_KO_clones=1,new_inferential_P='not_performed_clone_genotype_confounded'))
pd.DataFrame(summ).to_csv(O/'D5_KO_gene_effects.tsv',sep='\t',index=False)
core=[x for x in S['genes'] if x in set(v.gene)]
modules=[]
for sample,group in mapping.items():
 values=v[v['sample'].eq(sample)&v.gene.isin(core)].ln1pCP10k
 modules.append(dict(sample=sample,group=group,score=float(values.mean()),covered_genes=len(core),full_seven_gene_score=len(core)==7))
mo=pd.DataFrame(modules);mo.to_csv(O/'D5_KO_module_individual.tsv',sep='\t',index=False)
effect=[]
for group in ['CGNL1_KO_22','CGN_KO_1G8','double_KO_21D3']:
 wt=mo[mo.group.eq('WT')].score;ko=mo[mo.group.eq(group)].score
 effect.append(dict(group=group,WT_mean=float(wt.mean()),KO_mean=float(ko.mean()),delta=float(ko.mean()-wt.mean()),covered_genes=len(core),full_seven_gene_score=len(core)==7,independent_KO_clones=1))
pd.DataFrame(effect).to_csv(O/'D5_KO_module_effects.tsv',sep='\t',index=False)
qual={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'dataset':'GSE228955','primary_file':p.name,'full_12_columns_verified':True,'source_species':'Canis lupus familiaris','tissue':'kidney MDCKII','counts_type':'author_normalized_CPM','KO_clones_per_genotype':1,'culture_samples_per_group':3,'title_is_stale_WT_only_but_samples_include_KO':True,'xlsx_WT_only_not_combined_with_full_matrix':True,'symbol_issues':issues,'core_symbols_used':core,'inferential_P_performed':False,'interpretation':'descriptive real genetic perturbation with clone/genotype confounding; not CRC causality'}
(O/'D5_KO_QUALIFICATION.json').write_text(json.dumps(qual,indent=2))
print('D5_COMPLETE',len(core),len(issues),flush=True)
