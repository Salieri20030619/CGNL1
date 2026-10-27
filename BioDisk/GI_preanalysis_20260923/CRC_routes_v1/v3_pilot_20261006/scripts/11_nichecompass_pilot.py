from pathlib import Path
import os,json,datetime,time,traceback
for k in list(os.environ):
 if 'proxy' in k.lower():os.environ.pop(k,None)
os.environ.update(NO_PROXY='*',no_proxy='*',OMP_NUM_THREADS='8',OPENBLAS_NUM_THREADS='8',MKL_NUM_THREADS='8',DO_NOT_TRACK='1',MLFLOW_ENABLE_SYSTEM_METRICS_LOGGING='false')
import numpy as np,pandas as pd,anndata as ad,scanpy as sc,torch
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from nichecompass.models import NicheCompass
from nichecompass.utils import add_gps_from_gp_dict_to_adata
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/v3_pilot_20261006');D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_v3_pilot_20261006')
os.environ['MLFLOW_TRACKING_URI']='file://'+str(D/'cache/mlruns')
torch.set_num_threads(8)
plan=json.loads((R/'plan_v3.json').read_text())
excl=set(['CGNL1']+sum([plan['panels'][x] for x in ['notch_response','adhesion','historical_junction']],[]))
pairs=[('DLL4','NOTCH1'),('DLL4','NOTCH4'),('JAG1','NOTCH1'),('JAG1','NOTCH4'),('VEGFA','KDR'),('PDGFB','PDGFRB'),('CXCL12','CXCR4'),('TGFB1','TGFBR2')]
arr=[ad.read_h5ad(x['path']) for x in json.loads((R/'results/TECHNICAL_INPUTS.json').read_text())]
a=ad.concat(arr,join='inner',pairwise=True)
# Explicit block diagonal graph; no edges connect sections.
a.obsp['spatial_connectivities']=sparse.block_diag([x.obsp['spatial_connectivities'] for x in arr],format='csr')
a.obs_names_make_unique()
genes=set(a.var_names);gps={}
for ligand,receptor in pairs:
 if ligand in genes and receptor in genes and ligand not in excl and receptor not in excl:
  gps[ligand+'_'+receptor]={'sources':[ligand],'targets':[receptor],'sources_categories':['ligand'],'targets_categories':['receptor']}
norm=a.X.copy().tocsr();tot=np.asarray(norm.sum(axis=1)).ravel();norm=sparse.diags(1e4/np.maximum(tot,1))@norm;norm.data=np.log1p(norm.data)
var=np.asarray(norm.power(2).mean(axis=0)).ravel()-np.asarray(norm.mean(axis=0)).ravel()**2
ranking=np.argsort(-var,kind='stable')
chosen=[str(a.var_names[i]) for i in ranking if str(a.var_names[i]) not in excl][:800]
chosen=sorted(set(chosen)|{g for v in gps.values() for g in v['sources']+v['targets']})
a=a[:,chosen].copy();a.layers['counts']=a.X.copy()
assert not excl.intersection(set(a.var_names))
add_gps_from_gp_dict_to_adata(gps,a,min_source_genes_per_gp=1,min_target_genes_per_gp=1,filter_genes_not_in_masks=False)
params={'method':'NicheCompass','source_version':'0.3.3','freeze_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'genes':a.n_vars,'spots':a.n_obs,'prior_programs':list(gps),'target_exclusions':sorted(excl),'epochs':30,'seeds':[20261006,20261007],'scope':'reduced_real_data_technical_fit_not_biological_validation','model_overrides':{'n_addon_gp':8,'n_hidden_encoder':32,'conv_layer_encoder':'gcnconv','lambda_edge_recon':100,'lambda_gene_expr_recon':1,'lambda_l1_addon':0,'n_sampled_neighbors':10},'new_patient_n':0}
(R/'results/NICHECOMPASS_PARAMETERS.json').write_text(json.dumps(params,indent=2))
runs=[];clusters=[]
for seed in params['seeds']:
 aa=a.copy();np.random.seed(seed);torch.manual_seed(seed);t=time.time()
 model=NicheCompass(adata=aa,counts_key='counts',cat_covariates_keys=['sample'],cat_covariates_no_edges=[True],n_addon_gp=8,n_hidden_encoder=32,conv_layer_encoder='gcnconv',seed=seed)
 model.train(n_epochs=30,n_epochs_all_gps=30,lambda_edge_recon=100,lambda_gene_expr_recon=1,lambda_l1_addon=0,edge_batch_size=2048,node_batch_size=512,use_cuda_if_available=False,n_sampled_neighbors=10)
 emb=np.asarray(model.adata.obsm['nichecompass_latent'])
 assert np.isfinite(emb).all() and emb.shape[0]==a.n_obs
 pd.DataFrame(emb,index=a.obs_names).to_csv(R/'results'/f'NicheCompass_latent_{seed}.tsv.gz',sep='\t')
 labels=KMeans(n_clusters=8,n_init=10,random_state=seed).fit_predict(emb);clusters.append(labels)
 pd.DataFrame({'spot':a.obs_names,'sample':a.obs['sample'].values,'technical_cluster':labels}).to_csv(R/'results'/f'NicheCompass_clusters_{seed}.tsv',sep='\t',index=False)
 torch.save(model.model.state_dict(),D/'inputs'/f'NicheCompass_weights_{seed}.pt')
 runs.append({'seed':seed,'elapsed_s':time.time()-t,'latent_shape':list(emb.shape),'finite':True})
 print('SEED_COMPLETE',runs[-1],flush=True)
ari=float(adjusted_rand_score(*clusters))
out={'exit_code':0,'status':'native_run_complete','runs':runs,'cross_seed_cluster_ARI':ari,'same_data_seed_stability_not_external_replication':True,'requires_further_tuning_before_biology':ari<.7,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(R/'results/NICHECOMPASS_COMPLETE.json').write_text(json.dumps(out,indent=2));print(json.dumps(out),flush=True)
