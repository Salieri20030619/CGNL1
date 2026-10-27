import pathlib,json,gzip,numpy as np,pandas as pd
from scipy import io,sparse
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');D=B/'cohort_verification_20260924';R=D/'GSE245552'
cfg=json.loads((B/'spatial_v3/background_validation/queue_20260924_reference_v3/CONFIG.json').read_text());samples=json.loads((D/'GSE245552_LM_SAMPLES.json').read_text());rows=[]
for s in samples:
 p={('features' if 'features' in u else 'barcodes' if 'barcodes' in u else 'matrix'):R/u.rsplit('/',1)[1] for u in s['urls']}
 ft=pd.read_csv(p['features'],sep='\t',header=None);sy=ft[1].astype(str).to_numpy(dtype=str);bc=pd.read_csv(p['barcodes'],sep='\t',header=None)[0].to_numpy()
 with gzip.open(p['matrix'],'rb') as f,threadpool_limits(limits=2):X=io.mmread(f).tocsr()
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
 rows.append({'dataset':'GSE245552','sample':s['sample'],'title':s['title'],'donor':'GSE245552:'+s['donor'],'raw_cells':len(bc),'QC_cells':int(good.sum()),'endo_cells':len(a),'CGNL1_positive':len(pos),'CGNL1_total_UMI_endo':int(cg[a].sum()),'matched_pairs':len(pairs),'eligible_pseudobulk50':len(a)>=50,'eligible_matched_state':len(a)>=50 and len(pairs)>=10})
 pd.DataFrame(rows).to_csv(D/'GSE245552_library_coverage.tsv',sep='\t',index=False);print(s['sample'],len(a),len(pos),len(pairs),flush=True)
r=pd.DataFrame(rows);d=r.groupby('donor',as_index=False).agg(libraries=('sample','nunique'),endo_total=('endo_cells','sum'),CGNL1_positive_total=('CGNL1_positive','sum'),eligible_pseudobulk50=('eligible_pseudobulk50','max'),eligible_matched_state=('eligible_matched_state','max'),best_lesion_endo=('endo_cells','max'),best_lesion_pairs=('matched_pairs','max'));d.to_csv(D/'GSE245552_donor_coverage.tsv',sep='\t',index=False)
res={'status':'verified_counts_and_coverage','libraries':len(r),'donors':len(d),'donors_with_any_lesion_50_endo':int(d.eligible_pseudobulk50.sum()),'donors_with_any_lesion_10_matched_pairs':int(d.eligible_matched_state.sum()),'donors_with_any_CGNL1_positive_endo':int((d.CGNL1_positive_total>0).sum()),'scope':'Frozen marker/QC rules, no association tests; multiple lesions collapsed by donor; primary qualification requires one lesion >=50 endothelial cells, not pooling low-cell lesions to pass.'}
(D/'GSE245552_COVERAGE.json').write_text(json.dumps(res,indent=2));(D/'NEW_ENDO_EXIT.json').write_text(json.dumps({'exit_code':0}));print(json.dumps(res))
