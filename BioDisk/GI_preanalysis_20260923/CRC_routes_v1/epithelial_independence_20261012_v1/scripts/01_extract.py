from pathlib import Path
import os
os.environ.update(OMP_NUM_THREADS='4',OPENBLAS_NUM_THREADS='4',MKL_NUM_THREADS='4')
import numpy as np,pandas as pd,json,gzip,tarfile,hashlib,datetime
from scipy import sparse,io
from threadpoolctl import threadpool_limits
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1')
R=B/'epithelial_independence_20261012_v1';O=R/'results';C=R/'cache'
spec=json.loads((R/'SPECIFICATION.json').read_text())
assert hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()==json.loads((R/'SPEC_FREEZE.json').read_text())['sha256']
panels={**spec['panels'],**spec['outcomes']}
needed=sorted({'CGNL1'}|set(spec['ordinary_marker_comparators'])|set(sum(panels.values(),[])))
def collapse(X,genes):
 u,inv=np.unique(genes,return_inverse=True)
 if len(u)!=len(genes):
  X=sparse.csr_matrix((np.ones(len(genes)),(inv,np.arange(len(genes)))),shape=(len(u),len(genes)))@X
  return X.tocsr(),u
 return X,genes
def scores(vals):
 out={}
 for name,gs in panels.items():
  ok=[vals[g] for g in gs if g in vals]
  out[name]=float(np.mean(ok)) if len(ok)>=.8*len(gs) else np.nan
  out[name+'_coverage']=len(ok)
 return out
reg=json.loads((B/'review_consensus_retest_20261002_v1/CACHE_INDEX.json').read_text())
anns={d:pd.read_csv(B/'review_consensus_retest_20261002_v1/results'/f'{d}_accepted_annotation.tsv.gz',sep='\t') for d in sorted({x['dataset'] for x in reg})}
rows=[];arrays=[];detects=[];gnames=[];fingerprints=[]
for rec in reg:
 print('SC_START',rec['sample'],flush=True)
 source=Path(rec['source_matrix']);stem=str(source).replace('matrix.mtx.gz','')
 genes=pd.read_csv(stem+'features.tsv.gz',sep='\t',header=None)[1].astype(str).to_numpy()
 bc=pd.read_csv(stem+'barcodes.tsv.gz',sep='\t',header=None)[0].astype(str)
 a=anns[rec['dataset']];a=a[a['sample'].eq(rec['sample'])].copy()
 ix=pd.Index(bc).get_indexer(a.barcode);assert (ix>=0).all() and a.barcode.is_unique
 with gzip.open(source,'rb') as f,threadpool_limits(limits=4):X=io.mmread(f).tocsr()
 X=X[:,ix];umi=np.asarray(X.sum(0)).ravel();assert np.array_equal(umi,a.UMI_QC)
 X,genes=collapse(X,genes)
 mask=a.accepted_major_lineage.eq('epithelial').to_numpy();n=int(mask.sum())
 tot=np.asarray(X[:,mask].sum(1)).ravel();detection=np.asarray((X[:,mask]>0).sum(1)).ravel()/max(n,1)
 den=float(tot.sum());cp=tot/max(den,1)*1e4
 vals={g:float(np.log1p(cp[genes==g].sum())) for g in needed if np.any(genes==g)}
 row=dict(dataset=rec['dataset'],sample=rec['sample'],patient=rec['donor'],cells=n,eligible=n>=20,total_UMI=den,depth=float(np.log1p(umi[mask]).mean()) if n else np.nan,**vals,**scores(vals))
 rows.append(row);arrays.append(tot);detects.append(detection);gnames.append(genes)
 fingerprints.append(dict(path=str(source),size=source.stat().st_size,mtime=source.stat().st_mtime,matched_annotation_cells=len(a),epi_cells=n))
 del X
 print('SC_DONE',rec['sample'],'epi_cells',n,flush=True)
common=sorted(set.intersection(*[set(x) for x in gnames]));cts=[];det=[]
for t,d,g in zip(arrays,detects,gnames):
 ix=pd.Index(g).get_indexer(common);cts.append(t[ix]);det.append(d[ix])
lib=pd.DataFrame(rows);lib.to_csv(O/'A_library_profiles.tsv',sep='\t',index=False)
np.savez_compressed(C/'A_all_gene_pseudobulk.npz',genes=np.array(common),counts=np.array(cts),detection=np.array(det),samples=lib['sample'].to_numpy(),denominators=lib.total_UMI.to_numpy())
pt=lib[lib.eligible].groupby(['dataset','patient'],as_index=False).mean(numeric_only=True)
pt.to_csv(O/'A_patient_profiles.tsv',sep='\t',index=False)
# Abundance and detection matching is frozen before functional association models.
cp=np.array(cts)/np.maximum(lib.total_UMI.to_numpy()[:,None],1)*1e4
temp=[];dtemp=[]
for (dataset,patient),ix in lib[lib.eligible].groupby(['dataset','patient']).groups.items():
 temp.append(cp[list(ix)].mean(0));dtemp.append(np.array(det)[list(ix)].mean(0))
means=np.mean(temp,0);dete=np.mean(dtemp,0);j=common.index('CGNL1')
excluded=set(needed);dist=np.c_[np.log(means+1e-6),np.log(np.clip(dete,1e-5,1-1e-5)/(1-np.clip(dete,1e-5,1-1e-5)))]
scale=dist.std(0);distance=np.linalg.norm((dist-dist[j])/np.maximum(scale,1e-12),axis=1)
tbl=pd.DataFrame({'gene':common,'mean_CP10k':means,'detection':dete,'distance':distance})
pool=tbl[(tbl.mean_CP10k>=means[j]/2)&(tbl.mean_CP10k<=means[j]*2)&(abs(tbl.detection-dete[j])<=.15)&~tbl.gene.isin(excluded)].sort_values(['distance','gene'])
chosen=pool.head(20) if len(pool)>=20 else pool.iloc[:0]
chosen.to_csv(O/'A_abundance_matched_controls.tsv',sep='\t',index=False)
tbl.to_csv(O/'A_matching_inputs.tsv.gz',sep='\t',index=False)
(O/'A_CONTROL_FREEZE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'eligible_controls':len(pool),'controls':chosen.gene.tolist(),'target_mean_CP10k':float(means[j]),'target_detection':float(dete[j]),'association_results_examined':False},indent=2))
# Supplementary all-gene patient profiles support frozen comparator fits.
logcp=np.log1p(cp);allpt=[]
for (dataset,patient),ix in lib[lib.eligible].groupby(['dataset','patient']).groups.items():
 allpt.append(pd.Series(logcp[list(ix)].mean(0),index=common,name=patient))
pd.DataFrame(allpt).to_csv(O/'A_patient_allgenes_logCP10k.tsv.gz',sep='\t',index_label='patient')
(R/'SOURCE_FINGERPRINTS_A.json').write_text(json.dumps(fingerprints,indent=2))
(O/'A_EXTRACTION_COMPLETE.json').write_text(json.dumps({'libraries':len(lib),'eligible_libraries':int(lib.eligible.sum()),'patients':len(pt),'controls':len(chosen)},indent=2))
print('A_EXTRACTION_COMPLETE',len(pt),flush=True)
# Existing paired primary/liver gene-level aggregates, without reannotating cells.
agg=np.load(B/'results/route1/aggregates.npz',allow_pickle=True)
meta=pd.read_csv(B/'results/route1/group_metadata.tsv',sep='\t');assert np.array_equal(meta.group_id,agg['group_ids'])
genes=agg['genes'].astype(str);cts=agg['counts'];den=cts.sum(0);L=np.log1p(cts/np.maximum(den,1)[None,:]*1e4);br=[]
for i,a in meta.iterrows():
 vals={g:float(L[genes==g,i].sum()) for g in needed if np.any(genes==g)}
 br.append(dict(a,**vals,**scores(vals)))
pd.DataFrame(br).to_csv(O/'B_author_group_profiles.tsv',sep='\t',index=False)
# Frozen target-excluded RNA-CNV sensitivity matrices.
Xd=B/'deep_v2/results/cnv';X=sparse.load_npz(Xd/'counts_gene_by_cell.npz')
gn=pd.read_csv(Xd/'gene_names.tsv',sep='\t').gene.to_numpy()
m=pd.read_csv(B/'metadata/nonimmune_meta.txt.gz',sep='\t',index_col=0);m.index=m.index.str.replace('-','.',regex=False)
cnv=pd.read_csv(Xd/'cell_CNV_classification.tsv.gz',sep='\t',index_col=0);assert X.shape[1]==len(m)
for q in [95,99]:
 subset=cnv[cnv['cnv_supported_'+str(q)]&~cnv.cluster.str.startswith('Tu10')]
 for (patient,site),a in subset.groupby(['patients','organs']):
  ix=m.index.get_indexer(a.index);assert (ix>=0).all()
  tot=np.asarray(X[:,ix].sum(1)).ravel();den=tot.sum()
  vals={g:float(np.log1p(tot[gn==g].sum()/max(den,1)*1e4)) for g in needed if np.any(gn==g)}
  br.append(dict(patient=patient,site=site,celltype='CNV'+str(q)+'_noTu10',n_cells=len(ix),**vals,**scores(vals)))
pd.DataFrame(br).to_csv(O/'B_all_group_profiles.tsv',sep='\t',index=False)
del X
(O/'B_EXTRACTION_COMPLETE.json').write_text(json.dumps({'group_rows':len(br)},indent=2))
print('B_EXTRACTION_COMPLETE',flush=True)
# Spatial measurement cache only. Pathology eligibility never uses CGNL1.
ann=pd.read_csv(B/'results/spatial_author_metadata.tsv.gz',sep='\t');ann['sample']=ann.Barcode.str.split('~').str[0];ann['barcode']=ann.Barcode.str.split('~').str[1]
jobs=[x for x in json.loads((B/'spatial_v3/metadata/LIVER_DOWNLOAD_PLAN.json').read_text()) if 'Recurrence' not in x['sample_type']]
qa=[]
for job in jobs:
 sample=job['sample'];print('SPATIAL_START',sample,flush=True)
 source=B/'spatial_v3/raw/liver'/f'{sample}.tar.gz'
 with tarfile.open(source,'r|gz') as tar:
  for member in tar:
   if not member.isfile():continue
   if member.name.endswith('filtered_feature_bc_matrix/features.tsv.gz'):
    with gzip.GzipFile(fileobj=tar.extractfile(member)) as f:genes=pd.read_csv(f,sep='\t',header=None)[1].astype(str).to_numpy()
   elif member.name.endswith('filtered_feature_bc_matrix/barcodes.tsv.gz'):
    with gzip.GzipFile(fileobj=tar.extractfile(member)) as f:bc=pd.read_csv(f,sep='\t',header=None)[0].astype(str).to_numpy()
   elif member.name.endswith('filtered_feature_bc_matrix/matrix.mtx.gz'):
    with gzip.GzipFile(fileobj=tar.extractfile(member)) as f,threadpool_limits(limits=4):X=io.mmread(f).tocsr()
 a=ann[ann['sample'].eq(sample)].copy().reset_index(drop=True)
 ix=pd.Index(bc).get_indexer(a.barcode);assert (ix>=0).all()
 X=X[:,ix];X.eliminate_zeros();umi=np.asarray(X.sum(0)).ravel();ng=X.getnnz(0)
 consistent=a.Layer2.eq('Liver_metastatic carcinoma').to_numpy()
 roi=a.Layer3.isin(spec['C']['regions']).to_numpy()
 keep=(umi>=500)&(ng>=200)&consistent&roi
 X,genes=collapse(X[:,keep],genes);umi=umi[keep];a=a[keep].reset_index(drop=True)
 present=[g for g in needed if g in genes];gi=pd.Index(genes).get_indexer(present)
 cnt=X[gi].toarray().T
 np.savez_compressed(C/f'C_{sample}.npz',genes=np.array(present),counts=cnt,umi=umi,region=a.Layer3.to_numpy(),barcode=a.barcode.to_numpy(),patient=np.array([job['patient']]),sample=np.array([sample]))
 row=dict(patient=job['patient'],sample=sample,macro=int((a.Layer3=='Liver macrometastasis tumor').sum()),micro=int((a.Layer3=='Liver micrometastasis tumor').sum()),excluded_conflicting_pathology=int((roi&~consistent).sum()),target_measured='CGNL1' in present,target_total=int(cnt[:,present.index('CGNL1')].sum()) if 'CGNL1' in present else -1)
 qa.append(row);pd.DataFrame(qa).to_csv(O/'C_region_qualification.tsv',sep='\t',index=False)
 del X
 print('SPATIAL_DONE',sample,'macro',row['macro'],'micro',row['micro'],flush=True)
(O/'EXTRACTION_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'A_patients':len(pt),'spatial_samples':len(qa),'source_spec_sha256':hashlib.sha256((R/'SPECIFICATION.json').read_bytes()).hexdigest()},indent=2))
print('EXTRACTION_COMPLETE',flush=True)
