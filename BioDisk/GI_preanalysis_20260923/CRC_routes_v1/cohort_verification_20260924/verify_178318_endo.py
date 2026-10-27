import pathlib,json,gzip,numpy as np,pandas as pd
from scipy import io,sparse
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');D=B/'cohort_verification_20260924';R=D/'GSE178318'
cfg=json.loads((B/'spatial_v3/background_validation/queue_20260924_reference_v3/CONFIG.json').read_text());rows=[]
import hashlib
manifest=json.loads((R/'MANIFEST.json').read_text())
for rec in manifest:
 h=hashlib.sha256()
 with (R/rec['file']).open('rb') as fh:
  for block in iter(lambda:fh.read(1048576),b''):h.update(block)
 assert h.hexdigest()==rec['sha256'],rec['file']
print('transfer_sha256_verified',flush=True)
ft=pd.read_csv(R/'GSE178318_genes.tsv.gz',sep='\t',header=None);sy=ft[1].astype(str).to_numpy(dtype=str);allbc=pd.read_csv(R/'GSE178318_barcodes.tsv.gz',sep='\t',header=None)[0].astype(str).to_numpy()
tags=np.array([x.split('_',1)[1] for x in allbc]);print(pd.Series(tags).value_counts().to_dict(),flush=True)
assert set(tags)=={p+'_'+t for p in ['COL07','COL12','COL15','COL16','COL17','COL18'] for t in ['CRC','LM']}|{'COL12_PBMC','COL17_PBMC','COL18_PBMC'}
with gzip.open(R/'GSE178318_matrix.mtx.gz','rb') as f,threadpool_limits(limits=2):ALL=io.mmread(f).tocsr()
assert ALL.shape==(len(sy),len(allbc)) and pd.Index(allbc).is_unique
for tag in sorted(set(tags)):
 if not tag.endswith('_LM'):continue
 s={'sample':tag,'title':tag,'donor':tag.split('_')[0]};ind=np.flatnonzero(tags==tag);X=ALL[:,ind].tocsr();bc=allbc[ind]
 X.eliminate_zeros();assert X.shape==(len(sy),len(bc)) and pd.Index(bc).is_unique;assert (X.data>=0).all() and np.allclose(X.data,np.round(X.data))
 lib=np.asarray(X.sum(axis=0)).ravel();ng=X.getnnz(axis=0);mt=np.asarray(X[np.char.startswith(sy,'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)*100;good=(lib>=500)&(ng>=200)&(mt<=25);scores=[];hits=[];labels=list(cfg['lineage_markers'])
 for gs in cfg['lineage_markers'].values():
  z=X[np.isin(sy,gs)].toarray();scores.append(np.log1p(z/np.maximum(lib,1)*1e4).mean(axis=0));hits.append((z>0).sum(axis=0))
 scores=np.array(scores).T;hits=np.array(hits).T;best=scores.argmax(axis=1);margin=np.sort(scores,axis=1)[:,-1]-np.sort(scores,axis=1)[:,-2];j=labels.index('endothelial');endo=good&(best==j)&(hits[:,j]>=2)&(margin>=.4);a=np.flatnonzero(endo);g=np.flatnonzero(sy=='CGNL1');assert len(g)==1;cg=X[g[0]].toarray().ravel();pos=a[cg[a]>0];neg=a[cg[a]==0];feat=np.c_[np.log1p(lib),np.log1p(ng)];used=set();pairs=[]
 for ip in sorted(pos,key=lambda i:(cg[i],i)):
  can=np.array([k for k in neg if k not in used],int)
  if len(can)==0:break
  dis=np.abs(feat[can]-feat[ip]);ok=(dis<=.25).all(axis=1)
  if not ok.any():continue
  can=can[ok];dis=dis[ok];k=can[np.argmin((dis**2).sum(axis=1))];used.add(k);pairs.append((ip,k))
 rows.append({'dataset':'GSE178318','sample':s['sample'],'title':s['title'],'donor':'GSE178318:'+s['donor'],'raw_cells':len(bc),'QC_cells':int(good.sum()),'endo_cells':len(a),'CGNL1_positive':len(pos),'CGNL1_total_UMI_endo':int(cg[a].sum()),'matched_pairs':len(pairs),'eligible_pseudobulk50':len(a)>=50,'eligible_matched_state':len(a)>=50 and len(pairs)>=10})
 pd.DataFrame(rows).to_csv(D/'GSE178318_library_coverage.tsv',sep='\t',index=False);print(s['sample'],len(a),len(pos),len(pairs),flush=True)
r=pd.DataFrame(rows);d=r.groupby('donor',as_index=False).agg(libraries=('sample','nunique'),endo_total=('endo_cells','sum'),CGNL1_positive_total=('CGNL1_positive','sum'),eligible_pseudobulk50=('eligible_pseudobulk50','max'),eligible_matched_state=('eligible_matched_state','max'),best_lesion_endo=('endo_cells','max'),best_lesion_pairs=('matched_pairs','max'));d.to_csv(D/'GSE178318_donor_coverage.tsv',sep='\t',index=False)
res={'status':'verified_counts_and_coverage','libraries':len(r),'donors':len(d),'donors_with_any_lesion_50_endo':int(d.eligible_pseudobulk50.sum()),'donors_with_any_lesion_10_matched_pairs':int(d.eligible_matched_state.sum()),'donors_with_any_CGNL1_positive_endo':int((d.CGNL1_positive_total>0).sum()),'scope':'Frozen marker/QC rules, no association tests; multiple lesions collapsed by donor; primary qualification requires one lesion >=50 endothelial cells, not pooling low-cell lesions to pass.'}
(D/'GSE178318_COVERAGE.json').write_text(json.dumps(res,indent=2));(D/'GSE178318_ENDO_EXIT.json').write_text(json.dumps({'exit_code':0}));print(json.dumps(res))
