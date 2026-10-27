
from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,gzip,csv,itertools,datetime,subprocess
import numpy as np,pandas as pd
from scipy import stats
p=Path.cwd();old=p.parent/'CGNL1_gastric_peritoneal_validation_20261016_v1'
launch=json.loads((p/'metadata/SPHEROID_DOWNLOAD_LAUNCH_231.json').read_text())
if not Path(launch['receipt']).exists():subprocess.run(['/usr/bin/python3',str(old/'scripts/cgnl_wait_pidfd_system.py'),str(launch['wrapper_pid'])],check=True)
assert Path(launch['receipt']).exists() and json.loads(Path(launch['receipt']).read_text())['exit_code']==0
assert json.loads((p/'metadata/DOWNLOAD_COMPLETE_231.json').read_text())['all_pass']
spec=json.loads((p/'SPHEROID_ANALYSIS_SPECIFICATION_233.json').read_text());tasks=json.loads((p/'metadata/DOWNLOAD_PLAN_231.json').read_text())['tasks'];needed=set(spec['context_genes']+['CGNL1'])
processed=p/'processed/spheroid_233';processed.mkdir(parents=True,exist_ok=True);out=p/'results/spheroid_233';out.mkdir(parents=True,exist_ok=True)
samples=[];totals_by_sample={};vectors={};cycle_rows=[]
for task in tasks:
 f=p/'raw'/task['name']
 with gzip.open(f,'rt',encoding='utf-8',newline='') as h:
  rd=csv.reader(h,delimiter=chr(9));header=next(rd);first=next(rd)
  cells=header[1:] if len(first)==len(header) else header
  assert len(first)==len(cells)+1 and len(set(cells))==len(cells) and all(cells)
  n=len(cells);lib=np.zeros(n,dtype=np.int64);ng=np.zeros(n,dtype=np.int64);gene_totals={};selected={}
  for row in itertools.chain([first],rd):
   assert len(row)==n+1
   gene=row[0];assert gene not in gene_totals
   v=np.asarray(row[1:],dtype=float);assert np.isfinite(v).all() and (v>=0).all() and np.equal(v,np.floor(v)).all()
   v=v.astype(np.int64);lib+=v;ng+=(v>0);gene_totals[gene]=int(v.sum())
   if gene in needed:selected[gene]=v
 assert (lib>0).all()
 acc=task['accession'];title=task['title'];totals_by_sample[acc]=gene_totals;vectors[acc]=selected
 np.savez_compressed(processed/(acc+'_cell_counts.npz'),cell_id=np.asarray(cells,dtype=str),provided_gene_UMI=lib,n_genes=ng,**selected)
 with gzip.open(processed/(acc+'_gene_totals.json.gz'),'wt',encoding='utf-8') as h:json.dump(gene_totals,h)
 target=selected.get('CGNL1')
 rec={'accession':acc,'title':title,'patient':title.split('-')[0],'condition':'spheroid' if title.endswith('-spheroid') else 'parent','cells':n,'provided_genes':len(gene_totals),'provided_UMI':int(lib.sum()),'target_measured':target is not None,'target_UMI':int(target.sum()) if target is not None else None,'target_detected_cells':int((target>0).sum()) if target is not None else None}
 samples.append(rec);print('PARSED',title,n,len(gene_totals),'target_measured',target is not None,flush=True)
 if target is not None and all(g in selected for g in ['MKI67','TOP2A']):
  positive=(selected['MKI67']>0)|(selected['TOP2A']>0)
  for label,mask in [('MKI67_or_TOP2A_detected',positive),('neither_detected',~positive)]:
   total=int(lib[mask].sum());targetsum=int(target[mask].sum())
   cycle_rows.append({'accession':acc,'patient':rec['patient'],'condition':rec['condition'],'group':label,'cells':int(mask.sum()),'CGNL1_UMI':targetsum,'provided_UMI':total,'CGNL1_log2CPM1':float(np.log2(1+targetsum/total*1e6)) if total else None})
common=set.intersection(*(set(d) for d in totals_by_sample.values()));assert len(common)>1000
common_sorted=sorted(common);accessions=[r['accession'] for r in samples]
gmat=np.array([[totals_by_sample[a][g] for a in accessions] for g in common_sorted],dtype=float)
sf_genes=[i for i,g in enumerate(common_sorted) if g!='CGNL1' and np.all(gmat[i]>0)]
gm=np.exp(np.log(gmat[sf_genes]).mean(axis=1));sf=np.median(gmat[sf_genes]/gm[:,None],axis=0);sf/=np.exp(np.log(sf).mean())
contexts=[]
for i,r in enumerate(samples):
 acc=r['accession'];denom=sum(totals_by_sample[acc][g] for g in common);r['common_genes']=len(common);r['common_UMI']=denom;r['median_ratio_size_factor']=float(sf[i])
 if r['target_measured']:
  r['CGNL1_common_log2CPM1']=float(np.log2(1+r['target_UMI']/denom*1e6));r['CGNL1_provided_log2CPM1']=float(np.log2(1+r['target_UMI']/r['provided_UMI']*1e6));r['CGNL1_medianratio_log2count1']=float(np.log2(1+r['target_UMI']/sf[i]))
 for g in spec['context_genes']:
  value=totals_by_sample[acc].get(g)
  contexts.append({'accession':acc,'patient':r['patient'],'condition':r['condition'],'gene':g,'measured':value is not None,'UMI':value,'common_log2CPM1':float(np.log2(1+value/denom*1e6)) if value is not None else None})
df=pd.DataFrame(samples);assert len(df)==6 and df.patient.nunique()==3 and not df[['patient','condition']].duplicated().any()
summaries=[];deltas=[]
if df.target_measured.all() and 'CGNL1' in common:
 for measure in ['CGNL1_common_log2CPM1','CGNL1_provided_log2CPM1','CGNL1_medianratio_log2count1']:
  wide=df.pivot(index='patient',columns='condition',values=measure);d=wide.spheroid-wide.parent
  row={'measure':measure,'n_independent_patients':3,'mean_difference':float(d.mean()),'median_difference':float(d.median()),'positive':int((d>0).sum()),'negative':int((d<0).sum()),'P':None,'role':'primary' if measure=='CGNL1_common_log2CPM1' else 'normalization_sensitivity'}
  if measure=='CGNL1_common_log2CPM1' and (d!=0).any():
   nz=d[d!=0].to_numpy();ranks=stats.rankdata(abs(nz));observed=abs(np.dot(np.sign(nz),ranks));row['P']=float(np.mean([abs(np.dot(sign,ranks))>=observed-1e-12 for sign in itertools.product([-1,1],repeat=len(nz))]))
  summaries.append(row)
  deltas.extend({'patient':patient,'measure':measure,'difference':float(value)} for patient,value in d.items())
else:
 summaries.append({'role':'primary','status':'not_evaluable_due_to_source_filtered_target_absence','missing_accessions':df.loc[~df.target_measured,'accession'].tolist(),'no_zero_imputation':True})
df.to_csv(out/'sample_target_summary.tsv',sep=chr(9),index=False,encoding='utf-8')
pd.DataFrame(summaries).to_csv(out/'paired_contrasts.tsv',sep=chr(9),index=False,encoding='utf-8')
pd.DataFrame(deltas).to_csv(out/'patient_deltas.tsv',sep=chr(9),index=False,encoding='utf-8')
pd.DataFrame(contexts).to_csv(out/'context_gene_summary.tsv',sep=chr(9),index=False,encoding='utf-8')
pd.DataFrame(cycle_rows).to_csv(out/'marker_detection_strata.tsv',sep=chr(9),index=False,encoding='utf-8')
(out/'COMMON_GENE_UNIVERSE.json').write_text(json.dumps(common_sorted),encoding='utf-8')
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_files':6,'independent_patients':3,'common_measured_genes':len(common),'all_target_measured':bool(df.target_measured.all()),'contrasts':summaries,'independent_numeric_QA_pending':True,'causal_validation':False,'whole_project_complete':False}
(out/'ANALYSIS_COMPLETE_233.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2),flush=True)
