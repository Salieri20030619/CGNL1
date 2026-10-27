from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import numpy as np,pandas as pd,json,gzip,hashlib,datetime
from scipy import io,sparse
from threadpoolctl import threadpool_limits
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'new_cohorts_20261013_v1';O=R/'results'
S=json.loads((R/'SPECIFICATION.json').read_text());assert hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text())['sha256']
needed=sorted(set(['CGNL1']+S['ordinary_markers']+sum(S['panels'].values(),[])+sum(S['outcomes'].values(),[])))
rows=[];cover=[];source=[]
def extract(X,genes,cells,meta):
 X=X.tocsr();X.eliminate_zeros()
 assert (X.data>=0).all() and np.allclose(X.data,np.round(X.data))
 umi=np.asarray(X.sum(0)).ravel();ng=X.getnnz(0);mt=np.asarray(X[np.char.startswith(genes.astype(str),'MT-')].sum(0)).ravel()/np.maximum(umi,1)*100
 good=(umi>=500)&(ng>=200)&(mt<=25)
 scores=[];hits=[];labels=list(S['lineage_markers'])
 for label,gs in S['lineage_markers'].items():
  a=X[np.isin(genes,gs)].toarray();assert len(a)
  scores.append(np.log1p(a/np.maximum(umi,1)*1e4).mean(0));hits.append((a>0).sum(0))
 scores=np.array(scores).T;hits=np.array(hits).T;j=labels.index('epithelial')
 margin=np.sort(scores,axis=1)[:,-1]-np.sort(scores,axis=1)[:,-2]
 epi=good&(scores.argmax(1)==j)&(hits[:,j]>=2)&(margin>=.4)
 cover.append(dict(meta,raw_cells=len(cells),qc_cells=int(good.sum()),epithelial_cells=int(epi.sum()),eligible=int(epi.sum())>=20))
 pd.DataFrame(dict(cell=cells,UMI=umi,genes=ng,mito_pct=mt,QC=good,frozen_epithelial=epi)).to_csv(O/(meta['dataset']+'_'+meta['sample']+'_cell_QC.tsv.gz'),sep='\t',index=False)
 if epi.sum()<20:return
 tot=np.asarray(X[:,epi].sum(1)).ravel();den=tot.sum()
 vals={g:float(np.log1p(tot[genes==g].sum()/den*1e4)) for g in needed if np.any(genes==g)}
 for name,gs in {**S['panels'],**S['outcomes']}.items():
  vv=[vals[g] for g in gs if g in vals];vals[name]=float(np.mean(vv)) if len(vv)>=.8*len(gs) else np.nan;vals[name+'_coverage']=len(vv)
 rows.append(dict(meta,definition='frozen_marker_QC',cells=int(epi.sum()),total_UMI=float(den),depth=float(np.log1p(umi[epi]).mean()),**vals))
D=B/'cohort_verification_20260924/GSE178318'
genes=pd.read_csv(D/'GSE178318_genes.tsv.gz',sep='\t',header=None)[1].astype(str).to_numpy(dtype=str)
bc=pd.read_csv(D/'GSE178318_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str).to_numpy()
tags=np.array([c.split('_',1)[1] for c in bc]);print('READ_GSE178318',flush=True)
with gzip.open(D/'GSE178318_matrix.mtx.gz','rb') as f,threadpool_limits(limits=4):X=io.mmread(f).tocsr()
assert X.shape==(len(genes),len(bc))
for tag in sorted(set(tags)):
 if tag.endswith('_PBMC'):continue
 ix=np.flatnonzero(tags==tag)
 extract(X[:,ix],genes,bc[ix],dict(dataset='GSE178318',sample=tag,patient='GSE178318:'+tag.split('_')[0],site='liver_metastasis' if tag.endswith('_LM') else 'primary_CRC',patient_mapping='explicit_COL_patient_barcode'))
 print('EXTRACT',tag,cover[-1]['epithelial_cells'],flush=True)
del X
D=B/'spatial_v3/background_validation/run_20260924T052010Z'
samples=json.loads((D/'SAMPLES.json').read_text())
for rec in samples:
 if rec['dataset']!='GSE231559' or not any('tumor' in s.lower() for s in rec['characteristics']):continue
 paths={('features' if 'features' in p else 'barcodes' if 'barcodes' in p else 'matrix'):D/p for p in rec['files']}
 genes=pd.read_csv(paths['features'],sep='\t',header=None)[1].astype(str).to_numpy(dtype=str)
 bc=pd.read_csv(paths['barcodes'],sep='\t',header=None)[0].astype(str).to_numpy()
 with gzip.open(paths['matrix'],'rb') as f,threadpool_limits(limits=4):X=io.mmread(f).tocsr()
 site='liver_metastasis' if any('liver tumor' in s.lower() for s in rec['characteristics']) else 'primary_CRC'
 extract(X,genes,bc,dict(dataset='GSE231559',sample=rec['sample'],patient=rec['patient_label_provisional'],site=site,patient_mapping='author_L_or_C_label; no C-L pairing assumed'))
 print('EXTRACT',rec['title'],cover[-1]['epithelial_cells'],flush=True)
 del X
pd.DataFrame(rows).to_csv(O/'liver_library_profiles.tsv',sep='\t',index=False)
pd.DataFrame(cover).to_csv(O/'liver_cell_coverage.tsv',sep='\t',index=False)
# Collapse repeated L8 tumor samples using count-scale gene amounts and UMI denominators.
pt=[]
for key,g in pd.DataFrame(rows).groupby(['dataset','patient','site']):
 den=g.total_UMI.sum();d=dict(zip(['dataset','patient','site'],key));d.update(cells=int(g.cells.sum()),libraries=len(g),definition='frozen_marker_QC',depth=float(np.average(g.depth,weights=g.cells)))
 for gene in needed:
  if gene in g:d[gene]=float(np.log1p(np.sum(np.expm1(g[gene])*g.total_UMI)/den))
 for name,gs in {**S['panels'],**S['outcomes']}.items():
  values=[d[x] for x in gs if x in d];d[name]=float(np.mean(values)) if len(values)>=.8*len(gs) else np.nan
 pt.append(d)
pd.DataFrame(pt).to_csv(O/'liver_patient_profiles.tsv',sep='\t',index=False)
(R/'LIVER_EXTRACTION_COMPLETE.json').write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),libraries=len(rows),patient_site_profiles=len(pt),all_pass=True),indent=2))
print('LIVER_EXTRACTION_COMPLETE',flush=True)
