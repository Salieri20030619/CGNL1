from pathlib import Path
import json,gzip,datetime,sys
import numpy as np,pandas as pd
from scipy import sparse,io,stats
from threadpoolctl import threadpool_limits
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[1];B=R.parent;D=R/'05_detection_controls';C=json.loads((R/'frozen/EXECUTION_CONFIG.json').read_text());old=B/'cohort_verification_20260924';reg=pd.read_csv(old/'qualified_endothelial_patients.tsv',sep='\t');qualified=set(reg.donor)
rows=[];matched=[];inputs=[]
def frozen_endo(X,sy):
 lib=np.asarray(X.sum(axis=0)).ravel();ng=X.getnnz(axis=0);mt=np.asarray(X[np.char.startswith(sy.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)*100;qc=(lib>=500)&(ng>=200)&(mt<=25);scores=[];hits=[];labels=list(C['lineage_markers'])
 for gs in C['lineage_markers'].values():
  z=X[np.isin(sy,gs)].toarray();scores.append(np.log1p(z/np.maximum(lib,1)*1e4).mean(axis=0));hits.append((z>0).sum(axis=0))
 scores=np.array(scores).T;hits=np.array(hits).T;best=scores.argmax(axis=1);margin=np.sort(scores,axis=1)[:,-1]-np.sort(scores,axis=1)[:,-2];j=labels.index('endothelial');return qc&(best==j)&(hits[:,j]>=2)&(margin>=.4)
cache=[]
def profile(dataset,donor,sample,X,sy,endo):
 if donor not in qualified:return
 k=np.flatnonzero(endo)
 if len(k)<50:return
 Y=X[:,k].tocsr();Y.eliminate_zeros();assert np.allclose(Y.data,np.round(Y.data))
 lib=np.asarray(Y.sum(axis=0)).ravel();ng=Y.getnnz(axis=0);name=dataset+'_'+sample;f=D/(name+'.npz');sparse.save_npz(f,Y);pd.DataFrame({'gene':sy}).to_csv(D/(name+'_genes.tsv'),sep='\t',index=False)
 cache.append(dict(dataset=dataset,donor=donor,sample=sample,matrix=str(f),genes=str(D/(name+'_genes.tsv')),cells=len(k)));print('CACHED',donor,flush=True)
def load_matrix(paths):
 f={k:next(v for v in paths if k in v.name) for k in ['features','barcodes','matrix']};sy=pd.read_csv(f['features'],sep='\t',header=None)[1].to_numpy(dtype=str);bc=pd.read_csv(f['barcodes'],sep='\t',header=None)[0].to_numpy(dtype=str)
 with gzip.open(f['matrix'],'rb') as h,threadpool_limits(limits=2):X=io.mmread(h).tocsr()
 X.eliminate_zeros();assert X.shape==(len(sy),len(bc));assert len(set(bc))==len(bc)
 inputs.extend(map(str,paths));return X,sy,bc
# discovery, original count matrix and author metadata order from earlier validated conversion
f=B/'deep_v2/results/cnv/counts_gene_by_cell.npz';X=sparse.load_npz(f).tocsr();sy=pd.read_csv(B/'deep_v2/results/cnv/gene_names.tsv',sep='\t').gene.to_numpy(dtype=str);meta=pd.read_csv(B/'metadata/nonimmune_meta.txt.gz',sep='\t',index_col=0);assert X.shape==(len(sy),len(meta));endo=frozen_endo(X,sy);inputs.append(str(f))
for donor,g in meta[meta.organs=='LCT'].groupby('patients'):
 ix=meta.index.get_indexer(g.index);profile('GSE225857','GSE225857:'+donor,donor,X[:,ix],sy,endo[ix])
del X
# Previously frozen external cell identities; no reannotation on target expression.
cfg=json.loads((B/'spatial_v3/background_validation/final_exploration_20260924_v2/CONFIG.json').read_text());I=Path(cfg['input_run']);S=Path(cfg['source_run'])
for s in json.loads((I/'SAMPLES.json').read_text()):
 if s['patient_label_provisional'] not in qualified:continue
 if s['dataset']=='GSE231559' and not (s['title'].startswith('L') and not s['title'].endswith('N')):continue
 X,sy,bc=load_matrix([I/f for f in s['files']]);lib=np.asarray(X.sum(axis=0)).ravel();ng=X.getnnz(axis=0);mt=np.asarray(X[np.char.startswith(sy,'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)*100;qc=(lib>=500)&(ng>=200)&(mt<=25);X=X[:,qc];bc=bc[qc]
 ann=pd.read_csv(S/'cells'/(s['sample']+'.tsv.gz'),sep='\t').set_index('barcode').reindex(bc);assert ann.label.notna().all()
 profile(s['dataset'],s['patient_label_provisional'],s['sample'],X,sy,ann.label.to_numpy()=='endothelial')
# Additional prior-qualified three-patient cohort, using the same frozen markers.
for s in json.loads((old/'GSE245552_LM_SAMPLES.json').read_text()):
 donor='GSE245552:'+s['donor']
 if donor not in qualified:continue
 X,sy,bc=load_matrix([old/'GSE245552'/u.rsplit('/',1)[1] for u in s['urls']]);profile('GSE245552',donor,s['sample'],X,sy,frozen_endo(X,sy))

(D/'CACHE_INDEX.json').write_text(json.dumps(cache,indent=2))
excluded=set(['CGNL1']+C['junction_panel']+sum(list(C['lineage_markers'].values()),[]));discover=[]
for v in cache:
 if v['dataset']!='GSE225857':continue
 Y=sparse.load_npz(v['matrix']);sy=pd.read_csv(v['genes'],sep='\t').gene.to_numpy(dtype=str);assert len(set(sy))==len(sy)
 discover.append(pd.DataFrame({'gene':sy,'mean_count':np.asarray(Y.mean(axis=1)).ravel(),'detection':Y.getnnz(axis=1)/Y.shape[1]}).set_index('gene'))
avg=sum(discover)/len(discover);cg=avg.loc['CGNL1'];assert (cg>0).all();ratios=np.log2(avg[['mean_count','detection']].replace(0,np.nan)/cg);ok=(ratios.abs()<=1).all(axis=1);ok&=~avg.index.isin(excluded);ok&=~avg.index.str.match(r'^(MT-|RPL|RPS)')
selected=avg.loc[ok].copy();selected['distance']=np.square(ratios.loc[ok]).sum(axis=1);selected=selected.sort_values(['distance']).head(20);selected.to_csv(D/'discovery_selected_controls.tsv',sep='\t');controls=selected.index.tolist();assert len(controls)>0
rows=[]
for v in cache:
 Y=sparse.load_npz(v['matrix']);sy=pd.read_csv(v['genes'],sep='\t').gene.to_numpy(dtype=str);idx={g:i for i,g in enumerate(sy)};lib=np.asarray(Y.sum(axis=0)).ravel();ng=Y.getnnz(axis=0);feat=np.c_[np.log1p(lib),np.log1p(ng)];panel=np.log1p(Y[[idx[g] for g in C['junction_panel']]].toarray()/lib*1e4).mean(axis=0)
 for gene in ['CGNL1']+controls:
  if gene not in idx:continue
  z=Y[idx[gene]].toarray().ravel();pos=np.flatnonzero(z>0);neg=np.flatnonzero(z==0);used=set();pairs=[]
  for i in sorted(pos,key=lambda k:(z[k],k)):
   can=np.array([k for k in neg if k not in used],int)
   if not len(can):break
   delta=np.abs(feat[can]-feat[i]);valid=(delta<=.25).all(axis=1)
   if valid.any():
    can=can[valid];delta=delta[valid];j=can[np.argmin((delta**2).sum(axis=1))];used.add(j);pairs.append((i,j))
  row=dict(dataset=v['dataset'],donor=v['donor'],sample=v['sample'],gene=gene,pairs=len(pairs),eligible=len(pairs)>=10,effect=np.nan)
  if len(pairs)>=10:
   a,b=np.array(pairs).T;row['effect']=float(np.mean(panel[a]-panel[b]));row['mean_abs_logUMI_difference']=float(np.abs(feat[a,0]-feat[b,0]).mean())
  rows.append(row)
t=pd.DataFrame(rows);t.to_csv(D/'all_control_lesion_results.tsv',sep='\t',index=False);eligible=t[t.eligible];pat=eligible.groupby(['dataset','donor','gene'],as_index=False).effect.mean();pat.to_csv(D/'patient_control_effects.tsv',sep='\t',index=False)
old=pd.read_csv(R/'04_patient_state/matched_patient_effects_descriptive.tsv',sep='\t').set_index('donor').effect;new=pat[pat.gene=='CGNL1'].set_index('donor').effect;assert set(old.index)==set(new.index);assert np.allclose(old.sort_index(),new.sort_index())
comp=[]
for g in controls:
 v=pat[pat.gene==g].set_index('donor').effect;common=new.index.intersection(v.index)
 if len(common):comp.append(dict(control_gene=g,common_patients=len(common),CGNL1_mean=float(new[common].mean()),control_mean=float(v[common].mean()),target_minus_control=float((new[common]-v[common]).mean()),target_greater_patients=int((new[common]>v[common]).sum()),all_target_patients=len(common)==len(new)))
comp=pd.DataFrame(comp);comp.to_csv(D/'common_patient_control_comparisons.tsv',sep='\t',index=False)
summary={'exit_code':0,'selected_control_genes':len(controls),'target_eligible_patients':len(new),'CGNL1_mean_matched_difference':float(new.mean()),'CGNL1_positive_patients':int((new>0).sum()),'comparators_with_all_target_patients':int(comp.all_target_patients.sum()),'comparators_with_at_least3_common_patients':int((comp.common_patients>=3).sum()),'target_above_control_among_ge3':int(((comp.common_patients>=3)&(comp.target_minus_control>0)).sum()),'scope':'descriptive technical specificity comparison, no functional causality; same-subtype sensitivity pending','finished':datetime.datetime.now().astimezone().isoformat()};(D/'COMPLETE.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)
