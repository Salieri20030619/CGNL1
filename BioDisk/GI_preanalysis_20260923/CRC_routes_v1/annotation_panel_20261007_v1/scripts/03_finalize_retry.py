from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8',NUMBA_NUM_THREADS='8')
import json,hashlib,importlib.util,datetime,time
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'annotation_panel_20261007_v1';D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_annotation_20261007_v1');O=B/'v3_pilot_20261006';OD=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
q=ad.read_h5ad(D/'cache/refined_log_normalized.h5ad');a=ad.read_h5ad(D/'inputs/full_counts.h5ad');spec=json.loads((R/'standards/ANNOTATION_SPECIFICATION.json').read_text());review=json.loads((R/'review/cluster_annotation_review.json').read_text());mapping={str(k):v['label'] for k,v in review['clusters'].items()}
assert set(mapping)==set(q.obs.review_cluster.astype(str))
o=q.obs.copy();o['cluster_identity']=o.review_cluster.astype(str).map(mapping);o['annotation']=o.cluster_identity.copy();o['annotation_confidence']='marker_cluster_supported';o['annotation_reason']='reviewed expression cluster and at least two measured core markers; confidence is qualitative'
panels={k:v[0].split() for k,v in spec['panels'].items()};panels['T_NK']=panels['T']+panels['NK'];panels['myeloid']=panels['myeloid']+panels['dendritic'];panels['epithelial_cholangiocyte']=panels['epithelial']
panels['endothelial']+=['PLVAP','CD34'];panels['hepatocyte_lineage']=panels['hepatocyte'];excluded=set(spec['excluded_targets']);H={}
for label,genes in panels.items():
 gs=[g for g in genes if g in a.var_names and g not in excluded];H[label]=np.asarray((q[:,gs].X>0).sum(1)).ravel()
def broad(x):return {'T':'T_NK','NK':'T_NK','hepatocyte':'hepatocyte_lineage','dendritic':'myeloid','cholangiocyte':'epithelial_cholangiocyte','epithelial':'epithelial_cholangiocyte'}.get(x,x)
refs=[o[k+'_prediction'].map(broad).to_numpy() for k in ['normal_liver','multiome']]
corehits=np.zeros(len(o),int)
for label in set(o.cluster_identity):
 mask=o.cluster_identity.eq(label).to_numpy()
 if label in H:corehits[mask]=H[label][mask]
o['core_markers_detected']=corehits
identity=o.cluster_identity.map(broad).to_numpy();agree1=refs[0]==identity;agree2=refs[1]==identity
# Transfer failed domain validation; retain conflicting calls only as audit flags, never as a veto.
conflict=(refs[0]==refs[1])&(refs[0]!=identity)&(o.normal_liver_margin.to_numpy()>=.1)&(o.multiome_margin.to_numpy()>=.1)
low=(corehits<2)|o.cluster_identity.eq('uncertain').to_numpy()
o.loc[low,'annotation']='uncertain';o.loc[low,'annotation_confidence']='unresolved';o.loc[low,'annotation_reason']='insufficient measured core-marker support or unresolved cluster'
o['reference_conflict_flag']=conflict
witnesses={'T_NK':('CD3D CD3E',2),'B':('MS4A1 CD19 CD79A',2),'plasma':('MZB1 TNFRSF17',2),'myeloid':('CD68 CD163 C1QA FCGR1A',3),'endothelial':('FLT1 CLEC14A PTPRB ROBO4 PLVAP',3),'mural':('RGS5 MYH11 ACTA2 TAGLN',3),'fibroblast':('COL1A1 COL1A2 DCN PDGFRA',3)}
mixed=np.zeros(len(o),bool)
for alternative,(gs,nmin) in witnesses.items():
 hit=np.asarray((q[:,gs.split()].X>0).sum(1)).ravel()>=nmin
 compatible=(identity==alternative)|np.isin(identity,['uncertain'])
 if alternative in ['mural','fibroblast']:compatible|=np.isin(identity,['mural','fibroblast'])
 mixed|=hit&~compatible
o['mixed_marker_flag']=mixed;o.loc[mixed,'annotation']='uncertain';o.loc[mixed,'annotation_confidence']='unresolved';o.loc[mixed,'annotation_reason']='multiple-lineage marker overlap; segmentation/ambient/true state unresolved'
for c in ['cluster03','cluster06','cluster10','review_cluster','cluster_identity','annotation','annotation_confidence','annotation_reason','core_markers_detected','reference_conflict_flag','mixed_marker_flag','normal_liver_prediction','normal_liver_score','normal_liver_margin','multiome_prediction','multiome_score','multiome_margin']:
 if c=='core_markers_detected':a.obs[c]=-1
 elif c in ['normal_liver_score','normal_liver_margin','multiome_score','multiome_margin']:a.obs[c]=np.nan
 elif c in ['reference_conflict_flag','mixed_marker_flag']:a.obs[c]=False
 else:a.obs[c]='QC_excluded'
 a.obs.loc[o.index,c]=o[c].to_numpy()
a.obs['sample']='GSM9815454';a.obs['patient']='M242';a.uns['annotation_specification_json']=json.dumps(spec);a.uns['cluster_review_json']=json.dumps(review)
a.obs.to_csv(R/'results/cell_annotations.tsv.gz',sep='\t',index_label='cell_id');a.write_h5ad(D/'GSM9815454_panel_annotated_v1.h5ad',compression='gzip')
counts=a.obs.annotation.value_counts().rename_axis('annotation').reset_index(name='cells');counts.to_csv(R/'results/annotation_counts.tsv',sep='\t',index=False)
pd.crosstab(o.cluster06,o.annotation).to_csv(R/'results/cluster_final_labels.tsv',sep='\t')
for key in ['cluster03','cluster10']:pd.crosstab(o[key],o.annotation).to_csv(R/f'results/{key}_annotation_composition.tsv',sep='\t')
old=ad.read_h5ad(OD/'inputs/GSM9815454_quiche_technical.h5ad');p=a[old.obs_names].copy();p.obs['candidate_identity']=pd.Categorical(p.obs.annotation);assert np.all(p.obs.qc_pass) and p.n_obs==6000;assert np.array_equal(p.obsm['spatial'],old.obsm['spatial'])
labels=p.obs.annotation.astype(str).to_numpy();pd.crosstab(old.obs.candidate_identity,p.obs.annotation).to_csv(R/'results/pilot_old_new_labels.tsv',sep='\t');p.obs.to_csv(R/'results/pilot6000_annotations.tsv.gz',sep='\t',index_label='cell_id')
labelhash=hashlib.sha256('\n'.join(i+'\t'+v for i,v in zip(p.obs_names,labels)).encode()).hexdigest();source=O/'vendor/quiche_graph.py';cfg=dict(method='QUICHE official unmodified graph.py spatial_niches_khop',sample='GSM9815454',biological_samples=1,cells=6000,selection='exact original pilot cell IDs and order',radii_um=[100,200,300],k=10,khop=3,min_cell_threshold=3,max_neighbors=30,labels_sha256=labelhash,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),annotation='panel-adapted cluster review plus per-cell evidence/conflict flags',differential_abundance_run=False,CGNL1='not measured in this panel',frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(R/'standards/QUICHE_RETRY_PARAMETERS.json').write_text(json.dumps(cfg,indent=2))
module=importlib.util.spec_from_file_location('quiche_native_graph',source);qu=importlib.util.module_from_spec(module);module.loader.exec_module(qu)
xy=p.obsm['spatial'];dist,ix=NearestNeighbors(n_neighbors=30,algorithm='kd_tree').fit(xy).kneighbors();runs=[];oldmap={'endothelial_identity':'endothelial','tcell':'T_NK','epithelial':'epithelial_cholangiocyte','hepatocyte':'hepatocyte_lineage','unclassified':'uncertain'}
for rad in [100,200,300]:
 t=time.time();niches,nn=qu.spatial_niches_khop(p,radius=rad,k=10,khop=3,labels_key='candidate_identity',fov_key='sample',min_cell_threshold=3,n_jobs=8)
 assert np.isfinite(niches.to_numpy()).all() and np.allclose(niches.sum(axis=1),1)
 for i in range(20):
  keep=ix[i][dist[i]<rad]
  if len(keep)>=3:
   expected=pd.Series(labels[keep]).value_counts(normalize=True).reindex(niches.columns,fill_value=0);assert np.allclose(niches.loc[p.obs_names[i]].to_numpy(),expected.to_numpy())
 niches.to_csv(R/f'results/QUICHE_niches_radius{rad}.tsv.gz',sep='\t');previous=pd.read_csv(O/f'results/QUICHE_niches_radius{rad}.tsv.gz',sep='\t',index_col=0).rename(columns=oldmap);cols=sorted(set(previous.columns)|set(niches.columns));before=previous.reindex(index=niches.index,columns=cols,fill_value=0);after=niches.reindex(columns=cols,fill_value=0);tv=.5*np.abs(before-after).sum(axis=1)
 runs.append(dict(radius_um=rad,cells_retained=len(niches),classes=len(niches.columns),row_sum_valid=True,independent_20_rows_passed=True,mean_neighborhood_total_variation=float(tv.mean()),median_neighborhood_total_variation=float(tv.median()),cells_with_30_neighbors_within_radius=int((dist[:,-1]<rad).sum()),elapsed_s=time.time()-t));pd.DataFrame({'old_mean_fraction':before.mean(),'new_mean_fraction':after.mean()}).to_csv(R/f'results/QUICHE_comparison_radius{rad}.tsv',sep='\t');print('RETRY_SCALE_COMPLETE',runs[-1],flush=True)
pd.DataFrame(runs).to_csv(R/'results/QUICHE_retry_summary.tsv',sep='\t',index=False)
# Diagnostic figures, separate descriptive cell counts from biological replication.
palette={'hepatocyte_lineage':'#c49a42','endothelial':'#247ba0','mural':'#6b4c9a','fibroblast':'#6e9c42','epithelial_cholangiocyte':'#e07b39','T_NK':'#d14b78','T':'#d14b78','NK':'#882255','B':'#4477aa','plasma':'#66ccee','myeloid':'#aa4499','mast':'#44aa99','neutrophil':'#999933','uncertain':'#aaaaaa','QC_excluded':'#dddddd'}
fig,axes=plt.subplots(1,2,figsize=(14,6))
for ax,obj,title in [(axes[0],a,'M242T: all 80,036 cells'),(axes[1],p,'Original 6,000-cell QUICHE window')]:
 for label in sorted(obj.obs.annotation.unique()):
  m=obj.obs.annotation.eq(label).to_numpy();ax.scatter(*obj.obsm['spatial'][m].T,s=.4 if obj.n_obs>6000 else 2,c=palette.get(label,'#111111'),label=f'{label} ({m.sum():,})',rasterized=True)
 ax.set_title(title);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xlabel('x (um)');ax.set_ylabel('y (um)');ax.legend(fontsize=6,markerscale=5,loc='best')
fig.tight_layout();fig.savefig(R/'figures/annotation_spatial.png',dpi=170);plt.close(fig)
means=pd.read_csv(R/'results/refined_cluster_mean_log.tsv',sep='\t',dtype=str);means=means.set_index(means.columns[0]).astype(float);fracs=pd.read_csv(R/'results/refined_cluster_detection_fraction.tsv',sep='\t',dtype=str);fracs=fracs.set_index(fracs.columns[0]).astype(float)
genes='HNF1A G6PC CYP2E1 APOB FGB PECAM1 VWF KDR CLEC4G RGS5 MYH11 ACTA2 COL1A1 DCN PDGFRA EPCAM KRT19 KRT7 CD3D CD3E GNLY NKG7 MS4A1 CD79A MZB1 CD68 CD163 C1QA LYZ CD1C TPSB2 FCGR3B'.split();genes=[g for g in genes if g in means and g not in excluded];standard=(means[genes]-means[genes].mean())/means[genes].std().replace(0,1)
fig,ax=plt.subplots(figsize=(14,max(9,len(means)*.3)))
for j,g in enumerate(genes):ax.scatter(np.repeat(j,len(means)),np.arange(len(means)),s=fracs[g].to_numpy()*100,c=standard[g].clip(-2,2),cmap='RdBu_r',vmin=-2,vmax=2)
ax.set_xticks(range(len(genes)),genes,rotation=90);ax.set_yticks(range(len(means)),[f'{cl}: {mapping[str(cl)]}' for cl in means.index]);ax.invert_yaxis();ax.set_title('Measured markers: size = detection fraction; color = gene-wise standardized mean');fig.tight_layout();fig.savefig(R/'figures/marker_dotplot.png',dpi=180);plt.close(fig)
summary=dict(exit_code=0,annotation_completed=True,quiche_retry_completed=True,all_cells=a.n_obs,qc_pass_cells=int(a.obs.qc_pass.sum()),annotation_counts=counts.set_index('annotation').cells.to_dict(),pilot_counts=p.obs.annotation.value_counts().to_dict(),old_pilot_unclassified=int(old.obs.candidate_identity.eq('unclassified').sum()),new_pilot_uncertain=int(p.obs.annotation.eq('uncertain').sum()),per_cell_reference_conflicts=int(conflict.sum()),runs=runs,CGNL1_status='not measured',biological_samples=1,no_inferential_tests=True,labels_sha256=labelhash,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(R/'results/ANNOTATION_AND_RETRY_COMPLETE.json').write_text(json.dumps(summary,indent=2));print('ALL_COMPLETE',json.dumps(summary),flush=True)
