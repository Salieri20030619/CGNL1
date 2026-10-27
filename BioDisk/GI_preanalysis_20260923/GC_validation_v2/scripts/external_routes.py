import pathlib,json,tarfile,gzip,numpy as np,pandas as pd
from scipy import io,stats
from threadpoolctl import threadpool_limits
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');V=B/'GC_validation_v2';O=V/'candidate_routes';plan=json.loads((O/'PLAN.json').read_text());candidates=plan['genes'];programs=plan['programs'];wanted=sorted(set(candidates+sum(programs.values(),[])+['COPZ1']));assert not set(candidates)&set(sum(programs.values(),[]))
rows=[];coverage=[];qc=[]
with tarfile.open(V/'02_raw/GSE228598_RAW.tar') as tar:
 for prefix in sorted(n[:-len('_matrix.mtx.gz')] for n in tar.getnames() if n.endswith('_matrix.mtx.gz')):
  def f(n):return gzip.GzipFile(fileobj=tar.extractfile(n))
  ft=pd.read_csv(f(prefix+'_features.tsv.gz'),sep='\t',header=None);sy=ft.iloc[:,1].astype(str).values;bc=pd.read_csv(f(prefix+'_barcodes.tsv.gz'),sep='\t',header=None)[0]
  with threadpool_limits(limits=2):X=io.mmread(f(prefix+'_matrix.mtx.gz')).tocsr()
  m=pd.read_csv(V/'03_processed'/(prefix+'_gating_refined.tsv.gz'),sep='\t');ix=pd.Index(bc).get_indexer(m.barcode);assert (ix>=0).all();X=X[:,ix];umi=np.asarray(X.sum(axis=0)).ravel();assert np.array_equal(umi,m.umi.values);gi=np.flatnonzero(np.isin(sy,wanted));Y=X[gi];gg=sy[gi];assert len(set(gg))==len(gg)
  for pname,gset in programs.items():coverage.append({'sample':prefix,'program':pname,'planned':len(gset),'measured':sum(g in gg for g in gset),'missing':';'.join(g for g in gset if g not in gg)})
  for lab,ix in m.groupby('label').groups.items():
   ix=np.asarray(list(ix));cnt=np.asarray(Y[:,ix].sum(axis=1)).ravel();det=np.asarray((Y[:,ix]>0).sum(axis=1)).ravel();lib=umi[ix].sum()
   for k,g in enumerate(gg):rows.append({'patient':m.patient.iloc[0],'sample':prefix,'sample_type':m.sample_type.iloc[0],'celltype':lab,'gene':g,'n_cells':len(ix),'sample_cells':len(m),'cell_fraction':len(ix)/len(m),'median_umi':float(np.median(umi[ix])),'library_size':lib,'counts':cnt[k],'n_detected':int(det[k]),'fraction_detected':det[k]/len(ix),'CPM':cnt[k]/max(lib,1)*1e6})
  qc.append({'sample':prefix,'patient':m.patient.iloc[0],'cells':len(m),'genes_requested':len(wanted),'genes_found':len(gg)});print(prefix,len(m),flush=True)
a=pd.DataFrame(rows);a['log2CPM']=np.log2(a.CPM+1);a.to_csv(O/'external_gene_compartment_values.tsv.gz',sep='\t',index=False);pd.DataFrame(coverage).to_csv(O/'program_coverage.tsv',sep='\t',index=False)
# Source distributions refer only to captured cells and UMI counts.
src=a[a.gene.isin(candidates)].copy();src['gene_count_share']=src.counts/src.groupby(['patient','gene']).counts.transform('sum').replace(0,np.nan);src.to_csv(O/'external_candidate_source_shares.tsv',sep='\t',index=False)
sumsrc=src.groupby(['gene','celltype']).agg(n_patients=('patient','nunique'),patients_ge20_cells=('n_cells',lambda x:int((x>=20).sum())),median_cell_fraction=('cell_fraction','median'),median_transcript_share=('gene_count_share','median'),median_detection=('fraction_detected','median')).reset_index();sumsrc.to_csv(O/'external_source_share_summary.tsv',sep='\t',index=False)
# Include zero-absent compartments in abundance summaries, but do not infer gene expression there.
patients=sorted(a.patient.unique());ctypes=sorted(a.celltype.unique());ab=a[['patient','sample_type','celltype','n_cells','sample_cells','cell_fraction']].drop_duplicates();ab.to_csv(O/'external_captured_abundance.tsv',sep='\t',index=False)
wide=a.pivot(index=['patient','sample_type','celltype','n_cells','median_umi'],columns='gene',values='log2CPM').reset_index();ps=[]
for _,r in wide.iterrows():
 for pname,gset in programs.items():
  measured=[g for g in gset if g in wide.columns and pd.notna(r[g])];assert len(measured)>=.8*len(gset);ps.append({'patient':r.patient,'sample_type':r.sample_type,'celltype':r.celltype,'n_cells':r.n_cells,'median_umi':r.median_umi,'program':pname,'score':float(r[measured].mean())})
ps=pd.DataFrame(ps);ps.to_csv(O/'external_program_scores.tsv',sep='\t',index=False)
res=[]
for ct in plan['compartments']:
 for g in candidates:
  for prog in programs:
   x=a[(a.celltype==ct)&(a.gene==g)&(a.n_cells>=20)][['patient','sample_type','n_cells','median_umi','log2CPM']];z=x.merge(ps[(ps.celltype==ct)&(ps.program==prog)][['patient','score']],on='patient',validate='one_to_one');n=len(z);r={'celltype':ct,'gene':g,'program':prog,'n_patients':n,'rho':np.nan,'P':np.nan,'partial_rho':np.nan,'P_partial':np.nan,'partial_df':np.nan}
   if n>=10 and z.log2CPM.nunique()>2 and z.score.nunique()>2:
    co=stats.spearmanr(z.log2CPM,z.score);r.update(rho=float(co.statistic),P=float(co.pvalue));cov=pd.get_dummies(z.sample_type,drop_first=True,dtype=float).values;C=np.column_stack([np.ones(n),cov,stats.rankdata(np.log1p(z.n_cells)),stats.rankdata(np.log1p(z.median_umi))]);rank=np.linalg.matrix_rank(C);rx=stats.rankdata(z.log2CPM);ry=stats.rankdata(z.score);ex=rx-C@np.linalg.lstsq(C,rx,rcond=None)[0];ey=ry-C@np.linalg.lstsq(C,ry,rcond=None)[0];df=n-rank-1
    if min(np.std(ex),np.std(ey))>1e-10 and df>=3:
     rr=float(np.corrcoef(ex,ey)[0,1]);t=rr*np.sqrt(df/max(1-rr*rr,1e-15));r.update(partial_rho=rr,P_partial=float(2*stats.t.sf(abs(t),df)),partial_df=df)
   res.append(r)
r=pd.DataFrame(res);assert len(r)==54
for pcol,qcol in [('P','FDR_54'),('P_partial','FDR_partial_54')]:
 r[qcol]=stats.false_discovery_control(r[pcol].fillna(1).values,method='bh');r.loc[r[pcol].isna(),qcol]=np.nan
r.to_csv(O/'external_state_program_correlations.tsv',sep='\t',index=False)
(O/'EXTERNAL_COMPLETE.json').write_text(json.dumps({'exit_code':0,'patients':len(patients),'samples':len(qc),'planned_tests':54,'valid_correlations':int(r.P.notna().sum()),'valid_partial':int(r.P_partial.notna().sum()),'significant_raw_FDR':int((r.FDR_54<.05).sum()),'significant_partial_FDR':int((r.FDR_partial_54<.05).sum()),'scope':'cross-patient program associations, no metastatic-state control; not replication of paired changes'},indent=2));print(r[(r.FDR_54<.1)|(r.FDR_partial_54<.1)].to_string(index=False));print(sumsrc[sumsrc.celltype.isin(['Macrophage_enriched','Mesothelial_like','Epithelial_like'])].to_string(index=False))
