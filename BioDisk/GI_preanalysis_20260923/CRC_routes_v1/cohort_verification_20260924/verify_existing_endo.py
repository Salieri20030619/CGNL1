import pathlib,json,numpy as np,pandas as pd
from scipy import sparse
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');D=B/'cohort_verification_20260924';D.mkdir(exist_ok=True);R=B/'deep_v2/results/cnv'
X=sparse.load_npz(R/'counts_gene_by_cell.npz');sy=pd.read_csv(R/'gene_names.tsv',sep='\t').gene.astype(str).to_numpy();m=pd.read_csv(B/'metadata/nonimmune_meta.txt.gz',sep='\t',index_col=0);assert X.shape[1]==len(m) and X.shape[0]==len(sy)
cfg=json.loads((B/'spatial_v3/background_validation/queue_20260924_reference_v3/CONFIG.json').read_text());labels=list(cfg['lineage_markers']);lib=np.asarray(X.sum(axis=0)).ravel();ng=X.getnnz(axis=0);mt=np.asarray(X[np.char.startswith(sy.astype(str),'MT-')].sum(axis=0)).ravel()/np.maximum(lib,1)*100;QC=(lib>=500)&(ng>=200)&(mt<=25);scores=[];hits=[]
for gs in cfg['lineage_markers'].values():
 z=X[np.isin(sy,gs)].toarray();scores.append(np.log1p(z/np.maximum(lib,1)*1e4).mean(axis=0));hits.append((z>0).sum(axis=0))
scores=np.array(scores).T;hits=np.array(hits).T;best=scores.argmax(axis=1);margin=np.sort(scores,axis=1)[:,-1]-np.sort(scores,axis=1)[:,-2];j=labels.index('endothelial');endo=QC&(best==j)&(hits[:,j]>=2)&(margin>=.4);author=QC&m.cluster.str.startswith('E').to_numpy();cg=X[np.flatnonzero(sy=='CGNL1')[0]].toarray().ravel();feat=np.c_[np.log1p(lib),np.log1p(ng)]
rows=[]
for donor,z in m[m.organs=='LCT'].groupby('patients'):
 ix=m.index.get_indexer(z.index);a=ix[endo[ix]];auth=ix[author[ix]];pos=a[cg[a]>0];neg=a[cg[a]==0];used=set();pairs=[]
 for ip in sorted(pos,key=lambda i:(cg[i],i)):
  cand=np.array([k for k in neg if k not in used],int)
  if len(cand)==0:break
  dis=np.abs(feat[cand]-feat[ip]);ok=(dis<=.25).all(axis=1)
  if not ok.any():continue
  cand=cand[ok];dis=dis[ok];k=cand[np.argmin((dis**2).sum(axis=1))];used.add(k);pairs.append((ip,k))
 rows.append({'dataset':'GSE225857','donor':'GSE225857:'+donor,'liver_metastasis_cells':len(ix),'author_endo_before_QC':int(m.iloc[ix].cluster.str.startswith('E').sum()),'author_endo_after_QC':len(auth),'frozen_marker_endo':len(a),'marker_author_overlap':len(set(a)&set(auth)),'CGNL1_positive_frozen':len(pos),'CGNL1_detected_pct':len(pos)/len(a)*100 if len(a) else np.nan,'matched_pairs':len(pairs),'eligible_pseudobulk50':len(a)>=50,'eligible_matched_state':len(a)>=50 and len(pairs)>=10})
pd.DataFrame(rows).to_csv(D/'GSE225857_endothelial_coverage.tsv',sep='\t',index=False)
old=pd.read_csv(B/'spatial_v3/background_validation/final_exploration_20260924_v2/F01_SOURCE_STATE/eligibility.tsv',sep='\t');old=old[(old.lineage=='endothelial')&(old.target=='CGNL1')].copy();old=old.rename(columns={'cells':'frozen_marker_endo','positive_cells':'CGNL1_positive_frozen','positive_pct':'CGNL1_detected_pct','eligible':'eligible_matched_state'});old['eligible_pseudobulk50']=old.frozen_marker_endo>=50;allrows=pd.concat([old, pd.DataFrame(rows)],ignore_index=True);allrows.to_csv(D/'existing_endothelial_coverage.tsv',sep='\t',index=False)
result={'status':'coverage_verified','GSE225857':rows,'combined_qualified_pseudobulk_donor_labels':int(allrows.eligible_pseudobulk50.sum()),'combined_qualified_positive_negative_donor_labels':int(allrows.eligible_matched_state.sum()),'scope':'Raw matrices and existing fixed QC/markers; no target-association tests or threshold changes. GSE225857 already discovery, not independent external validation. Individual labels require study/sample provenance review.'}
(D/'EXISTING_COVERAGE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
