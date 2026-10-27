from pathlib import Path
import json,csv,itertools,datetime,collections
import numpy as np,pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
p=Path.cwd();out=p/'results/paired_pdx_310';out.mkdir(exist_ok=False);spec=json.loads((p/'metadata/STATISTICAL_SPECIFICATION_309.json').read_text());pair=json.loads((p/'metadata/SAMPLE_PAIRING_296b.json').read_text());a=np.load(p/'processed/matrix_input_306.npz',allow_pickle=False);ids=list(a['probe_ids']);accessions=list(a['accessions']);values=a['values'];lookup={v:i for i,v in enumerate(ids)}
platform={};inside=False;h=None
with (p/'raw/GPL10332_self_297.soft').open(encoding='utf-8') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row:continue
  if row[0]=='!platform_table_begin':inside=True;continue
  if row[0]=='!platform_table_end':break
  if not inside:continue
  if h is None:h=row;continue
  d=dict(zip(h,row));platform[d['ID']]=d
experimental=[i for i,k in enumerate(ids) if platform[k]['CONTROL_TYPE']=='FALSE'];negative=[i for i,k in enumerate(ids) if platform[k]['CONTROL_TYPE']=='neg'];assert len(experimental)==43118 and len(negative)==153;assert np.isfinite(values).all() and (values>0).all()
log=np.log2(values[experimental]);ordered=np.sort(log,axis=0);reference=ordered.mean(axis=1);normal=np.empty_like(log)
for j in range(log.shape[1]):
 order=np.argsort(log[:,j],kind='mergesort');x=log[order,j];_,starts,counts=np.unique(x,return_index=True,return_counts=True);means=np.add.reduceat(reference,starts)/counts;normal[order,j]=np.repeat(means,counts)
assert np.max(np.abs(normal.mean(axis=0)-reference.mean()))<1e-10
exp_index={ids[v]:i for i,v in enumerate(experimental)};manifest=spec['probe_manifest'];selected={r['ID']:normal[exp_index[r['ID']]] for r in manifest};pd.DataFrame(selected,index=accessions).rename_axis('accession').to_csv(out/'normalized_selected_features.tsv',sep='\t')
genes=sorted({r['GENE_SYMBOL'] for r in manifest});gene_values={}
for gene in genes:
 groups=collections.defaultdict(list)
 for r in manifest:
  if r['GENE_SYMBOL']==gene:groups[r['NAME']].append(selected[r['ID']])
 gene_values[gene]=np.mean([np.mean(v,axis=0) for v in groups.values()],axis=0)
df=pd.DataFrame(gene_values,index=accessions).rename_axis('accession').reset_index();meta=pd.DataFrame(pair['samples']);df=df.merge(meta[['accession','patient_id','condition','tumour_id','original_title']],on='accession',validate='one_to_one');df.to_csv(out/'array_gene_values.tsv',sep='\t',index=False)
tumours=df.groupby(['patient_id','condition','tumour_id'],sort=True)[genes].mean().reset_index();counts=df.groupby(['patient_id','condition','tumour_id']).size().rename('n_arrays').reset_index();tumours=tumours.merge(counts,validate='one_to_one');tumours.to_csv(out/'tumour_gene_values.tsv',sep='\t',index=False)
patient=tumours.groupby(['patient_id','condition'],sort=True)[genes].mean().reset_index();assert len(patient)==30 and patient.patient_id.nunique()==10;patient.to_csv(out/'patient_condition_values.tsv',sep='\t',index=False)
contrasts={'resistant_minus_acute':('resistant','acute'),'acute_minus_untreated':('acute','untreated'),'resistant_minus_untreated':('resistant','untreated')};deltas=[]
for pid,g in patient.groupby('patient_id'):
 indexed=g.set_index('condition')
 for name,(end,start) in contrasts.items():
  for gene in genes:deltas.append({'patient_id':pid,'contrast':name,'gene':gene,'difference':float(indexed.loc[end,gene]-indexed.loc[start,gene])})
delta=pd.DataFrame(deltas);delta.to_csv(out/'patient_deltas.tsv',sep='\t',index=False)
def exact(values):
 v=np.round(np.asarray(values),12);v=v[v!=0]
 if not len(v):return 1.
 ranks=np.rint(2*rankdata(abs(v),method='average')).astype(int);observed=abs(np.dot(np.sign(v).astype(int),ranks));return sum(abs(np.dot(signs,ranks))>=observed for signs in itertools.product([-1,1],repeat=len(v)))/(2**len(v))
def result(name,v,role):
 v=np.asarray(v);return dict(contrast=name,role=role,n_patients=len(v),mean_difference=float(v.mean()),median_difference=float(np.median(v)),positive=int((v>0).sum()),negative=int((v<0).sum()),P=exact(v))
results=[]
for name in contrasts:
 vals=delta[(delta.gene=='CGNL1')&(delta.contrast==name)].difference.to_numpy();results.append(result(name,vals,'primary' if name=='resistant_minus_acute' else 'secondary'))
secondary=results[1:];q=multipletests([r['P'] for r in secondary],method='fdr_bh')[1]
for r,fdr in zip(secondary,q):r['FDR']=float(fdr);r['family']='two_secondary_CGNL1_contrasts'
pivot=delta[delta.contrast=='resistant_minus_acute'].pivot(index='patient_id',columns='gene',values='difference');relative=pivot.CGNL1-pivot[['EPCAM','CDH1','KRT8']].mean(axis=1);results.append(result('resistant_minus_acute_target_minus_epithelial_mean',relative.to_numpy(),'lineage_sensitivity'));relative.rename('difference').to_csv(out/'lineage_relative_patient_deltas.tsv',sep='\t')
pd.DataFrame(results).to_csv(out/'target_contrasts.tsv',sep='\t',index=False)
context=[]
for (gene,name),g in delta.groupby(['gene','contrast']):
 v=g.difference.to_numpy();context.append({'gene':gene,'contrast':name,'mean_difference':float(v.mean()),'positive':int((v>0).sum()),'negative':int((v<0).sum())})
pd.DataFrame(context).to_csv(out/'context_description.tsv',sep='\t',index=False)
target_id=next(r['ID'] for r in manifest if r['GENE_SYMBOL']=='CGNL1');bg=np.quantile(values[negative],.99,axis=0);raw=values[lookup[target_id]];qc=pd.DataFrame({'accession':accessions,'CGNL1_raw_signal':raw,'negative_control_99pct':bg,'signal_to_background99':raw/bg,'above_negative99':raw>bg}).merge(meta[['accession','patient_id','condition']],on='accession',validate='one_to_one');qc.to_csv(out/'target_background_QC.tsv',sep='\t',index=False)
pos=delta[(delta.gene=='DUSP6')&(delta.contrast=='acute_minus_untreated')];positive_control={'patients':len(pos),'decreased':int((pos.difference<0).sum()),'patient_deltas':dict(zip(pos.patient_id,pos.difference)),'published_expectation':'decrease in all10models'}
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'arrays':len(df),'patients':10,'resistant_lesions':int((tumours.condition=='resistant').sum()),'normalization_features':len(experimental),'primary_and_secondary':results,'target_above_negative99_by_condition':qc.groupby('condition').above_negative99.agg(['sum','count']).to_dict('index'),'positive_control':positive_control,'test_rounding_decimals':12,'independent_QA_pending':True,'causal_validation':False,'whole_project_complete':False};(out/'ANALYSIS_COMPLETE_310.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x),flush=True)
