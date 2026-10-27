from pathlib import Path
import json,datetime,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
R=Path(__file__).resolve().parents[1];B=R.parent;D=R/'02_multiome';C=json.loads((R/'frozen/EXECUTION_CONFIG.json').read_text())
a=ad.read_h5ad(B/'expanded_screen_20260925/multiome/CRCLM_multiome_GEX_decontaminated.h5ad');X=a.X.tocsr().astype(np.float64);X.eliminate_zeros();m=a.obs.copy();m['cell_id']=m.index;genes=a.var_names.to_numpy(dtype=str);lib=np.asarray(X.sum(axis=1)).ravel();assert (lib>0).all()
qc=(m.nCount_RNA.to_numpy()>=500)&(m.nFeature_RNA.to_numpy()>=200)&(m['percent.mt'].to_numpy()<=25)
labels=list(C['lineage_markers']);scores=[];hits=[]
for label in labels:
 ix=np.flatnonzero(np.isin(genes,C['lineage_markers'][label]));assert len(ix)>0
 z=X[:,ix].toarray();scores.append(np.log1p(z/lib[:,None]*1e4).mean(axis=1));hits.append((z>=1).sum(axis=1))
scores=np.array(scores).T;hits=np.array(hits).T;best=scores.argmax(axis=1);margin=np.sort(scores,axis=1)[:,-1]-np.sort(scores,axis=1)[:,-2];supported=(margin>=.4)&(hits[np.arange(len(m)),best]>=2)
m['frozen_identity']=np.where(supported,np.array(labels)[best],'unresolved');m['QC']=qc;m['identity_margin']=margin;m['endo_markers_ge1']=hits[:,labels.index('endothelial')];m['adjusted_total']=lib
m['author_endo']=m.Cell_type.astype(str).eq('Endothelial');m['high_confidence_endo']=m.author_endo&qc&(m.frozen_identity=='endothelial')
m.to_csv(D/'cell_identity_QC.tsv.gz',sep='\t',index=False)
gs=['CGNL1']+C['junction_panel'];ix=[]
for g in gs:
 k=np.flatnonzero(genes==g);assert len(k)==1;ix.append(k[0])
vals=X[:,ix].toarray();norm=np.log1p(vals/lib[:,None]*1e4)
rows=[];source=[];subtypes=[]
for pt in sorted(m.Patient.astype(str).unique()):
 pm=m.Patient.astype(str).eq(pt).to_numpy();base={'patient':'CRCLM_multiome:'+pt,'sample':';'.join(m.loc[pm,'Sample'].astype(str).unique()),'therapy':';'.join(m.loc[pm,'Therapy'].astype(str).unique()),'tissue':'LM','independence':'cross_study_pending'}
 for definition,mask in [('author_endo_QC',pm&qc&m.author_endo.to_numpy()),('high_confidence_endo',pm&m.high_confidence_endo.to_numpy())]:
  n=int(mask.sum());row=dict(base,definition=definition,cells=n,coverage50=n>=50,unit='log1p_adjusted_CP10K',A06_raw_detection='unavailable')
  if n:
   pb=np.log1p(vals[mask].sum(axis=0)/lib[mask].sum()*1e4)
   for j,g in enumerate(gs):row[g+'_pseudobulk']=pb[j];row[g+'_mean_cell']=float(norm[mask,j].mean())
   row['junction_pseudobulk']=float(pb[1:].mean());row['junction_mean_cell']=float(norm[mask,1:].mean());row['median_original_UMI']=float(m.loc[mask,'nCount_RNA'].median())
  rows.append(row)
 for ct in sorted(m.Cell_type.astype(str).unique()):
  mask=pm&qc&m.Cell_type.astype(str).eq(ct).to_numpy();n=int(mask.sum())
  if n:source.append(dict(base,lineage=ct,cells=n,CGNL1_mean_cell=float(norm[mask,0].mean()),junction_mean_cell=float(norm[mask,1:].mean()),annotation='author_not_malignancy_validated'))
 for st,g in m.loc[pm&m.author_endo.to_numpy()&qc].groupby('Cell_subtype',observed=True):
  subtypes.append(dict(base,author_subtype=str(st),cells=len(g),author_endo_fraction=len(g)/int((pm&m.author_endo.to_numpy()&qc).sum())))
pd.DataFrame(rows).to_csv(D/'patient_continuous_profiles.tsv',sep='\t',index=False);pd.DataFrame(source).to_csv(D/'patient_lineage_expression.tsv',sep='\t',index=False);pd.DataFrame(subtypes).to_csv(D/'patient_author_endo_subtypes.tsv',sep='\t',index=False)
pd.crosstab(m.loc[qc,'Cell_type'],m.loc[qc,'frozen_identity']).to_csv(D/'author_vs_frozen_identity.tsv',sep='\t')
s=m.loc[m.high_confidence_endo].copy();out=a[s.index].copy();out.obs=s;out.write_h5ad(D/'high_confidence_endo_adjusted_expression.h5ad')
t=pd.DataFrame(rows);summary={'status':'completed_continuous_expression_qualification','matrix_scale':'decontX_adjusted_counts','original_raw_UMI_detection_A06':'not_available','patients_total':15,'author_endo_cells':int(m.author_endo.sum()),'author_endo_QC_cells':int((m.author_endo&qc).sum()),'high_confidence_endo_cells':int(m.high_confidence_endo.sum()),'author_endo_QC_patients_ge50':int(((t.definition=='author_endo_QC')&t.coverage50).sum()),'high_confidence_patients_ge50':int(((t.definition=='high_confidence_endo')&t.coverage50).sum()),'target_P_values_run':False,'pending':'patient independence, final export-chain/source QC audit; no claim of independent validation','finished':datetime.datetime.now().astimezone().isoformat()}
(D/'COMPLETE.json').write_text(json.dumps(dict(summary,exit_code=0),indent=2));print(json.dumps(summary),flush=True)
