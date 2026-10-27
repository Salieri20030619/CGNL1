import os
os.environ['OMP_NUM_THREADS']='2';os.environ['OPENBLAS_NUM_THREADS']='2';os.environ['NUMBA_NUM_THREADS']='2'
import pathlib,sys,json,gzip,re,numpy as np,pandas as pd
from scipy import sparse
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2');B=O.parent;R=O/'results/cnv';R.mkdir(exist_ok=True);sys.path.insert(0,str(O/'vendor'));import anndata as ad,scanpy as sc,infercnvpy as cnv
plan=json.loads((O/'PLAN.json').read_text());held=set(g.replace('.','-') for g in plan['frozen_genes']);held.update(['ALB','FGA','FGB','FGG','APOA1','APOA2','TTR','HP','AGT','GC','CRP']);(R/'CNV_PLAN.json').write_text(json.dumps({'method':'infercnvpy0.6.1; window100,step10,lfc_clip3,dynamic_threshold1.5; per-patient inference','QC':'nFeature>=200; mitochondrial<=25%; author metadata singlets','reference':'Fibroblast and endothelial author labels; each category>=40 cells, seeded random50% train50% held-out calibration; if unavailable category not used','cutoff':'Mean absolute RNA-inferred CNV deviation; tumor candidates above maximum within-reference-type95th percentile in held-out calibration;99th sensitivity. Not DNA-confirmed malignancy.','heldout_genes':sorted(held),'gene_filter':'autosomes; mapped NCBI GRCh38.p14; original gene mean counts>=0.1; candidate and hepatic-proxy genes excluded','min_subset_cells':20,'seed':24},indent=2))
m=pd.read_csv(B/'metadata/nonimmune_meta.txt.gz',sep='\t',index_col=0);m.index=m.index.str.replace('-','.',regex=False);assert m.index.is_unique
cache=R/'counts_gene_by_cell.npz'
if cache.exists():X=sparse.load_npz(cache);genes=pd.read_csv(R/'gene_names.tsv',sep='\t').gene.values
else:
 parts=[];names=[]
 for k,c in enumerate(pd.read_csv(B/'raw/nonimmune_counts.txt.gz',sep='\t',index_col=0,chunksize=128)):
  x=c.loc[:,m.index].values;assert (x>=0).all() and np.all(x==np.floor(x));parts.append(sparse.csr_matrix(x,dtype=np.float32));names.extend(c.index.astype(str))
  if k%50==0:print('SPARSE_READ',len(names),flush=True)
 X=sparse.vstack(parts,format='csr');genes=np.array(names);sparse.save_npz(cache,X);pd.DataFrame({'gene':genes}).to_csv(R/'gene_names.tsv',sep='\t',index=False)
assert np.array_equal(np.asarray(X.sum(axis=0)).ravel(),m.nCount_RNA.values)
coords=[]
with gzip.open(O/'metadata/GRCh38.p14.ncbi.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  c=line.rstrip().split('\t')
  if c[2]!='gene':continue
  mat=re.fullmatch(r'NC_0000(\d{2})\.\d+',c[0])
  if not mat or not 1<=int(mat.group(1))<=22:continue
  sym=re.search(r'(?:^|; )gene "([^"]+)"',c[8]);sym=sym.group(1) if sym else re.search(r'gene_id "([^"]+)"',c[8]).group(1);coords.append({'symbol':sym,'chromosome':'chr'+str(int(mat.group(1))),'start':int(c[3]),'end':int(c[4])})
co=pd.DataFrame(coords);co=co[~co.symbol.duplicated(keep=False)].set_index('symbol');canon=np.array([g.replace('.','-') for g in genes]);mean=np.asarray(X.mean(axis=1)).ravel();keep=np.isin(canon,co.index)&~np.isin(canon,list(held))&(mean>=.1);vi=np.flatnonzero(keep);var=co.reindex(canon[keep]).copy();var.index=genes[keep];var.to_csv(R/'genomic_gene_order.tsv',sep='\t');assert not set(canon[keep])&held;print('CNV_GENES',len(vi),flush=True)
rng=np.random.default_rng(24);outputs=[];profiles=[];coverage=[]
for pt in sorted(m.patients.unique()):
 ix=np.flatnonzero((m.patients==pt)&(m['percent.mt']<=25)&(m.nFeature_RNA>=200));mm=m.iloc[ix].copy();mm['broad']=np.where(mm.cluster.str.startswith('Tu'),'Tumor',np.where(mm.cluster.str.startswith('F'),'Fibroblast','Endothelial'));mm['reference_group']='nonreference';mm['reference_calibration']=False;refs=[]
 for ct in ['Fibroblast','Endothelial']:
  idx=np.flatnonzero(mm.broad==ct)
  if len(idx)>=40:
   idx=rng.permutation(idx);train=idx[:len(idx)//2];cal=idx[len(idx)//2:];name=ct+'_train';mm.iloc[train,mm.columns.get_loc('reference_group')]=name;mm.iloc[cal,mm.columns.get_loc('reference_calibration')]=True;refs.append(name)
 assert refs,pt
 XX=X[vi][:,ix].T.tocsr().multiply((1e4/m.iloc[ix].nCount_RNA.values)[:,None]).tocsr();XX.data=np.log1p(XX.data);aa=ad.AnnData(XX,obs=mm,var=var);print('CNV_PATIENT',pt,aa.shape,refs,flush=True);cnv.tl.infercnv(aa,reference_key='reference_group',reference_cat=refs,window_size=100,step=10,dynamic_threshold=1.5,chunksize=500,n_jobs=2);C=aa.obsm['X_cnv'];score=np.asarray(abs(C).mean(axis=1)).ravel();mm['cnv_score']=score;cuts={}
 for quant in [.95,.99]:
  cuts[quant]=max(np.quantile(score[(mm.broad==ct)&mm.reference_calibration],quant) for ct in ['Fibroblast','Endothelial'] if np.sum((mm.broad==ct)&mm.reference_calibration)>0);mm['cnv_supported_'+str(int(quant*100))]=(mm.broad=='Tumor')&(score>cuts[quant]);mm['cnv_cutoff_'+str(int(quant*100))]=cuts[quant]
 for (site,ct),ii in mm.groupby(['organs','broad']).indices.items():
  ii=np.asarray(ii);profiles.append({'patient':pt,'site':site,'broad':ct,'n_cells':len(ii),'profile':np.asarray(C[ii].mean(axis=0)).ravel().tolist()})
 mm.index.name='cell_id';outputs.append(mm);sparse.save_npz(R/(pt+'_cnv_profiles.npz'),C);(R/(pt+'_cnv_chromosomes.json')).write_text(json.dumps(aa.uns['cnv'],default=lambda x:int(x) if isinstance(x,np.integer) else str(x)));coverage.append({'patient':pt,'cells':len(mm),'ref_train':int(mm.reference_group.ne('nonreference').sum()),'ref_calibration':int(mm.reference_calibration.sum()),'cutoff95':float(cuts[.95]),'cutoff99':float(cuts[.99]),'tumor_author_QC':int((mm.broad=='Tumor').sum()),'supported95':int(mm.cnv_supported_95.sum()),'supported99':int(mm.cnv_supported_99.sum())});print(coverage[-1],flush=True)
cm=pd.concat(outputs);cm.to_csv(R/'cell_CNV_classification.tsv.gz',sep='\t');pd.DataFrame(coverage).to_csv(R/'patient_CNV_QC.tsv',sep='\t',index=False);(R/'group_profiles.json').write_text(json.dumps(profiles));counts=cm.groupby(['patients','organs']).agg(author_tumor=('broad',lambda x:int((x=='Tumor').sum())),CNV95=('cnv_supported_95','sum'),CNV99=('cnv_supported_99','sum')).reset_index();counts.to_csv(R/'paired_CNV_cell_coverage.tsv',sep='\t',index=False);(R/'COMPLETE.json').write_text(json.dumps({'exit_code':0,'patients':len(coverage),'CNV_genes':len(vi),'QC_cells':len(cm),'method_version':cnv.__version__,'scope':'RNA-inferred CNV support only; thresholds calibrated using held-out nonmalignant references; all126 candidates held out from inference'},indent=2));print(counts.to_string(index=False))
