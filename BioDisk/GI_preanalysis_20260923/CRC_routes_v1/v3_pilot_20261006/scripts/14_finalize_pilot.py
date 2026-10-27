from pathlib import Path
import os,json,hashlib,datetime,sys
os.environ.update(OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8')
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/v3_pilot_20261006');D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
refs=json.loads((R/'results/TECHNICAL_INPUTS.json').read_text());arr=[ad.read_h5ad(x['path']) for x in refs]
obs=pd.concat([a.obs for a in arr]);xy=np.vstack([a.obsm['spatial'] for a in arr]);A=sparse.block_diag([a.obsp['spatial_connectivities'] for a in arr],format='csr')
deg=np.asarray(A.sum(axis=1)).ravel();W=sparse.diags(1/np.maximum(deg,1))@A
keys=['endothelial_identity','mural','fibroblast','epithelial','hepatocyte','tcell','myeloid']
base=np.hstack([obs[keys].to_numpy(),np.asarray(W@obs[keys].to_numpy())]);base=StandardScaler().fit_transform(base)
baseline=[KMeans(n_clusters=8,n_init=10,random_state=x).fit_predict(base) for x in [20261006,20261007]]
nc=[pd.read_csv(R/'results'/f'NicheCompass_clusters_{x}.tsv',sep='\t').set_index('spot').loc[obs.index,'technical_cluster'].to_numpy() for x in [20261006,20261007]]
rr,cc=A.nonzero();rng=np.random.default_rng(20261006)
def metrics(labels):
 same=float(np.mean(labels[rr]==labels[cc]));chance=float(sum((np.bincount(labels)/len(labels))**2));return same,chance
m=[]
for name,labs in [('marker_neighborhood_baseline',baseline),('NicheCompass_reduced',nc)]:
 same,chance=metrics(labs[0]);m.append({'method':name,'cross_seed_ARI':float(adjusted_rand_score(*labs)),'neighbor_same_cluster_fraction':same,'cluster_size_expected_same_fraction':chance,'CGNL1_used_in_clustering':False,'ground_truth_accuracy':'not_measured','inference_P':'not_run'})
pd.DataFrame(m).to_csv(R/'results/TECHNICAL_BASELINE_COMPARISON.tsv',sep='\t',index=False)
fig,axes=plt.subplots(2,3,figsize=(11,7.5))
for i,a in enumerate(arr):
 mask=obs['sample'].values==a.obs['sample'].iloc[0];p=xy[mask]
 for j,(name,lab) in enumerate([('Marker-neighborhood baseline',baseline[0]),('NicheCompass seed 1',nc[0]),('NicheCompass seed 2',nc[1])]):
  ax=axes[i,j];ax.scatter(p[:,0],-p[:,1],c=lab[mask],cmap='tab10',vmin=0,vmax=9,s=7);ax.set_aspect('equal');ax.axis('off');ax.set_title(a.obs['sample'].iloc[0]+'\n'+name,fontsize=10)
fig.suptitle('Unsupervised technical clusters; colors are arbitrary across runs',fontsize=12);fig.tight_layout();fig.savefig(R/'results/nichecompass_baseline_maps.png',dpi=160);plt.close(fig)
Q=ad.read_h5ad(D/'inputs/GSM9815454_quiche_technical.h5ad');q=Q.obs.candidate_identity.astype(str)
fig,axes=plt.subplots(1,2,figsize=(10,4.6));cats=sorted(set(q))
for i,c in enumerate(cats):
 mask=q==c;axes[0].scatter(Q.obsm['spatial'][mask,0],-Q.obsm['spatial'][mask,1],s=3,label=c)
axes[0].legend(fontsize=6,loc='upper left',bbox_to_anchor=(1,1));axes[0].set_aspect('equal');axes[0].set_title('Provisional marker labels')
n=pd.read_csv(R/'results/QUICHE_niches_radius200.tsv.gz',sep='\t',index_col=0)
im=axes[1].scatter(Q.obsm['spatial'][:,0],-Q.obsm['spatial'][:,1],c=n.loc[Q.obs_names,'endothelial_identity'],s=4,cmap='viridis');axes[1].set_aspect('equal');axes[1].set_title('QUICHE: endothelial-neighbor fraction');fig.colorbar(im,ax=axes[1],shrink=.7)
for ax in axes:ax.axis('off')
fig.suptitle('One-patient technical pilot; CGNL1 not provided by the feature panel',fontsize=11);fig.tight_layout();fig.savefig(R/'results/quiche_technical_maps.png',dpi=160);plt.close(fig)
# Transport success does not imply a valid h5ad file.
p=D/'downloads/iss_cells.h5ad';content={'file':str(p),'transport_success':True,'hdf5_magic_valid':p.read_bytes()[:8]==b'\x89HDF\r\n\x1a\n','status':'invalid_input_html_viewer_response','used_in_analysis':False}
(R/'results/ISS_CONTENT_VALIDATION.json').write_text(json.dumps(content,indent=2))
status=[
{'method':'Spatial boundary gradient','status':'measurement_complete_region_model_not_run','evidence':'four patient-labelled sections;9353QC spots;independent pathology masks absent'},
{'method':'LIANA','status':'native_run_complete','evidence':'two sections;8 and3 pairs passed expression coverage;spatial shuffle changes scores'},
{'method':'NicheCompass','status':'native_run_complete_technical_only','evidence':'3000spots;812genes;8prior+8addonGP;two30epoch seeds;ARI0.7428;no claim of superiority'},
{'method':'QUICHE','status':'native_partial_run','evidence':'official unmodified graph module;6000cells;100/200/300um;81.8%unclassified;noDA;onepatient;CGNL1absent'},
{'method':'SMURF','status':'input_not_qualified_not_run','evidence':'oldZenodo19043057 and latest21431802 restricted;no publicHDfile list'},
{'method':'FlowSig','status':'conditional_qualification_only_not_run','evidence':'cell-specific inflow-GEM-target and suitable conditions not established;LIANA completion does not establish these'},
{'method':'OSDR','status':'conditional_qualification_only_not_run','evidence':'division and migration/death assumptions not qualified'},
{'method':'SpaCEy','status':'input_not_qualified_not_run','evidence':'no adequate independent spatial-protein clinical cohort'}
]
pd.DataFrame(status).to_csv(R/'results/METHOD_STATUS.tsv',sep='\t',index=False)
import importlib.metadata
versions={}
for pkg in ['python','numpy','pandas','scipy','scanpy','anndata','liana','nichecompass','torch','torch-geometric','torch-sparse','mlflow-skinny']:
 try:versions[pkg]=importlib.metadata.version(pkg) if pkg!='python' else sys.version
 except:versions[pkg]='unknown'
(R/'results/FINAL_ENVIRONMENT.json').write_text(json.dumps(versions,indent=2))
freeze=json.loads((R/'PLAN_FREEZE.json').read_text());checks={f:hashlib.sha256((R/f).read_bytes()).hexdigest()==h for f,h in freeze['files'].items()}
assert all(checks.values())
manifest={'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'v3_plan_and_initial_method_data_validation','full_research_complete':False,'pilot_complete':True,'network_mode':'direct_only','plan_freeze_verified':checks,'native_finished':['LIANA','NicheCompass'],'native_partial':['QUICHE official spatial module'],'boundary_validation_not_complete':True,'files':{}}
for folder in ['results','scripts','metadata']:
 for p in sorted((R/folder).rglob('*')):
  if p.is_file() and not p.name.endswith('.part'):
   manifest['files'][str(p.relative_to(R))]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
(R/'PILOT_COMPLETION.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({'baseline_comparison':m,'method_status':status,'freeze_verified':True}),flush=True)
