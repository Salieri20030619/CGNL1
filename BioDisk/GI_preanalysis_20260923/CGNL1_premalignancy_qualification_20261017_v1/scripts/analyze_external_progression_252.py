from pathlib import Path
import json,csv,gzip,re,datetime
import numpy as np,pandas as pd,statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
p=Path.cwd();q=p/'processed/external_251';out=p/'results/external_252';out.mkdir(exist_ok=False)
assert json.loads((q/'QUALIFICATION_COMPLETE_251.json').read_text())['samples']==45
ann=json.loads((q/'CGNL1_annotation_candidates.json').read_text());eligible=[];mapping=[];diagnostic=[]
for row in ann['rows']:
 assert len(row)==len(ann['header'])
 d=dict(zip(ann['header'],row));syms={v.strip() for v in re.split(r'///|;',d['Gene Symbol']) if v.strip()};ids=set(re.findall(r'ENSG[0-9]+',d['Ensembl']))
 annotated=syms=={'CGNL1'} and ids=={'ENSG00000128849'}
 if annotated:diagnostic.append(d['ID'])
 ok=annotated and d['ID']=='ADXEC.12191.C1_s_at' and 'transcript_biotype:protein_coding' in d['Transcript Assignments'] and 'NM_032866' in d['Transcript Assignments'] and d['Alignments'].startswith('chr15:')
 mapping.append(dict(probe=d['ID'],symbols=d['Gene Symbol'],ensembl=d['Ensembl'],eligible=ok,reason='coding CGNL1 probe qualified in amendment255' if ok else 'not primary; see specificity amendment255'))
 if ok:eligible.append(d['ID'])
assert eligible and len(eligible)==len(set(eligible)),mapping
(out/'probe_mapping.json').write_text(json.dumps(mapping,indent=2),encoding='utf-8')
meta=json.loads((q/'sample_metadata_rows.json').read_text());single={};chars=[]
for row in meta:
 if row[0]=='!Sample_characteristics_ch1':chars.append(row[1:])
 else:single[row[0]]=row[1:]
records=[]
for i,gsm in enumerate(single['!Sample_geo_accession']):
 d=dict(accession=gsm,title=single['!Sample_title'][i])
 for row in chars:
  k,v=row[i].split(': ',1);d[k.lower()]=v
 records.append(d)
df=pd.DataFrame(records).set_index('accession');assert df.title.nunique()==45
assert df.status.value_counts().to_dict()=={'Healthy Control':15,'CIM Control':9,'CIM-GC':8,'IIM Control':7,'IIM-GC':6}
values={};header=None
with gzip.open(p/'raw/external_250/GSE78523_series_matrix.txt.gz','rt',encoding='utf-8',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row or row[0].startswith('!'):continue
  if header is None:header=row;continue
  if row[0] in diagnostic:
   assert row[0] not in values;values[row[0]]=np.asarray(row[1:],dtype=float);assert np.isfinite(values[row[0]]).all()
assert set(values)==set(diagnostic)
for probe in diagnostic:df[probe]=pd.Series(values[probe],index=header[1:])
assert df[eligible].notna().all().all();df['CGNL1_mean_log2']=df[eligible].mean(axis=1);df['age']=df['age'].astype(float);df['man']=df['sex'].map({'Man':1.,'Woman':0.});assert df.man.notna().all()
df['progression']=df.status.isin(['IIM-GC','CIM-GC']).astype(float);df['IIM']=df.status.str.startswith('IIM').astype(float);df.to_csv(out/'sample_source_values.tsv',sep='\t')
results=[]
def fit(label,d,cols,target,y='CGNL1_mean_log2',role='primary'):
 X=pd.DataFrame({'intercept':np.ones(len(d))},index=d.index)
 for col in cols:X[col]=d[col]-d[col].mean() if col=='age' else d[col]
 assert np.linalg.matrix_rank(X.to_numpy())==X.shape[1]
 m=sm.OLS(d[y],X).fit(cov_type='HC3',use_t=True);ci=m.conf_int().loc[target];hat=np.sum(X.to_numpy()*(X.to_numpy()@np.linalg.inv(X.T@X)),axis=1)
 r=dict(contrast=label,role=role,outcome=y,n=len(d),model_columns=','.join(X.columns),effect=float(m.params[target]),CI_low=float(ci.iloc[0]),CI_high=float(ci.iloc[1]),P=float(m.pvalues[target]),df_resid=float(m.df_resid),max_leverage=float(max(hat)),target=target,n_group1=int(d[target].sum()),n_group0=int(len(d)-d[target].sum()))
 results.append(r);pd.concat([d[[y]],X],axis=1).to_csv(out/(label+'_model_values.tsv'),sep='\t');return r
im=df[df.status!='Healthy Control'].copy();fit('primary_progression_subtype',im,['progression','IIM'],'progression');fit('sensitivity_progression_demographic',im,['progression','IIM','age','man'],'progression',role='sensitivity')
for subtype in ['IIM','CIM']:
 d=im[im.status.str.startswith(subtype)].copy();fit('secondary_'+subtype+'_progression',d,['progression'],'progression',role='secondary')
d=df[~df.status.isin(['IIM-GC','CIM-GC'])].copy();d['IM']= (d.status!='Healthy Control').astype(float);fit('secondary_IM_vs_healthy',d,['IM','age','man'],'IM',role='secondary')
for i,probe in enumerate(diagnostic):fit('probe_'+str(i+1),im,['progression','IIM'],'progression',y=probe,role='probe_sensitivity')
for role in ['secondary','probe_sensitivity']:
 selected=[r for r in results if r['role']==role];adj=multipletests([r['P'] for r in selected],method='fdr_bh')[1]
 for r,v in zip(selected,adj):r['FDR']=float(v);r['family']=role
pd.DataFrame(results).to_csv(out/'target_contrasts.tsv',sep='\t',index=False)
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary':results[0],'all_results':results,'eligible_probes':eligible,'diagnostic_probes':diagnostic,'QA_pending':True,'batch_provenance':'Source processing mentions ComBat; no batch IDs or unequivocal stage of deposited matrix, so this cannot establish batch-independent validation.','causal_validation':False,'whole_project_complete':False}
(out/'ANALYSIS_COMPLETE_252.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x),flush=True)
