from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4')
import json,datetime,hashlib
import numpy as np,pandas as pd
from scipy import stats
import statsmodels.api as sm
O=Path.cwd();I=json.loads((O/'IMPLEMENTATION_FREEZE.json').read_text());panels=I['bulk_panels']
plan={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'task':'reproduction of published source, not external validation','target_effects_in_CHP_not_yet_seen':True,'labels':'use supplied filenames12with_PC/19without_PC consistent with paper;folderconflictretained;internalCHPheaderslackclinicalgroup;nooutcome-basedrelabelling','measurement':'CHP declared log2 scaled-RPM; actualrange0-18 consistent; useasprovided,nosecondlog','primary':'CGNL1 with_PC-minus-without_PC,Welch two-sided t test','secondary_models':['M0 PC-only HC3','M1 PC plus processing-version categorical covariates','M2 M1 plus six fixed source proxies from preceding spec'],'FDR':'three HC3 model PC-effect tests as sensitivity family; genome-wide Welch BH provided as contextual evidence, not sole candidate gate','healthy':'23 healthyfiles descriptiveonly;noassumed pairing','gene_duplicates':'retain assayrowidentity;targetrequireexactly1row;panelduplicatesmeanwithinsymbol;retain mappingissues','group_quali_limit':'filename-definedobservationalgroups,individualclinicalcovariatesnotverified'}
(O/'ORIGINAL31_ANALYSIS_FREEZE.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
X=pd.read_csv(O/'processed/original31_vendor_quantification.tsv.gz',sep='\t',index_col=0);meta=pd.read_csv(O/'metadata/original31_sample_manifest.tsv',sep='\t');dups=X.index[X.index.duplicated(False)];pd.Series(dups).value_counts().rename_axis('symbol').to_csv(O/'metadata/CHP_duplicate_symbols.tsv',sep='\t');target=X.loc[X.index=='CGNL1'];assert len(target)==1
meta['CGNL1']=target[meta.filename].to_numpy().ravel();cov=X.groupby(level=0).mean()
for name,genes in panels.items():
 avail=[g for g in genes if g in cov.index];assert len(avail)>=3;meta[name]=cov.loc[avail,meta.filename].mean(0).to_numpy()
meta.to_csv(O/'processed/original31_sample_values.tsv',sep='\t',index=False);tum=meta[meta.group_from_filename.isin(['with_PC','without_PC'])].copy();a=tum[tum.group_from_filename=='with_PC'].CGNL1.to_numpy();b=tum[tum.group_from_filename=='without_PC'].CGNL1.to_numpy();assert len(a)==12 and len(b)==19;tt=stats.ttest_ind(a,b,equal_var=False);se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b));df=(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))**2/((a.var(ddof=1)/len(a))**2/(len(a)-1)+(b.var(ddof=1)/len(b))**2/(len(b)-1));ci=stats.t.interval(.95,df,loc=a.mean()-b.mean(),scale=se)
with_cols=tum[tum.group_from_filename=='with_PC'].filename.tolist();without_cols=tum[tum.group_from_filename=='without_PC'].filename.tolist();tests=stats.ttest_ind(X[with_cols].to_numpy(),X[without_cols].to_numpy(),axis=1,equal_var=False);P=np.nan_to_num(tests.pvalue,nan=1);FDR=stats.false_discovery_control(P);gs=pd.DataFrame({'feature_row':np.arange(len(X)),'symbol':X.index,'effect':X[with_cols].mean(axis=1).to_numpy()-X[without_cols].mean(axis=1).to_numpy(),'P_Welch':P,'FDR_all_features':FDR});gs.to_csv(O/'results/original31_all_feature_context.tsv.gz',sep='\t',index=False)
result={'n_with_PC':12,'n_without_PC':19,'with_mean':float(a.mean()),'without_mean':float(b.mean()),'difference_log2_scaledRPM':float(a.mean()-b.mean()),'geometric_fold':float(2**(a.mean()-b.mean())),'P_Welch':float(tt.pvalue),'P_equal_variance_t':float(stats.ttest_ind(a,b,equal_var=True).pvalue),'P_MannWhitney':float(stats.mannwhitneyu(a,b,alternative='two-sided').pvalue),'CI_low':float(ci[0]),'CI_high':float(ci[1]),'FDR_all_features':float(gs.loc[gs.symbol=='CGNL1','FDR_all_features'].iloc[0]),'family_size':len(X),'published_FC':4.82,'published_P':.0081,'labels':'filename defined, folder mismatch flagged','LOO_effect_min':float(min([np.delete(a,i).mean()-b.mean() for i in range(len(a))]+[a.mean()-np.delete(b,i).mean() for i in range(len(b))])),'LOO_effect_max':float(max([np.delete(a,i).mean()-b.mean() for i in range(len(a))]+[a.mean()-np.delete(b,i).mean() for i in range(len(b))]))}
models=[];ver=pd.get_dummies(tum.program_version,drop_first=True,dtype=float);pc=tum.group_from_filename.eq('with_PC').astype(float).to_numpy()
for name,cc in [('M0',[]),('M1',list(ver.columns)),('M2',list(ver.columns)+list(panels))]:
 columns=[ver[c].to_numpy() if c in ver else tum[c].to_numpy() for c in cc];Z=np.column_stack([np.ones(len(tum)),pc]+columns);rank=np.linalg.matrix_rank(Z);df=len(tum)-rank
 if rank<Z.shape[1] or df<5:models.append({'model':name,'status':'not_estimable','P':np.nan});continue
 fit=sm.OLS(tum.CGNL1.to_numpy(),Z).fit(cov_type='HC3',use_t=True);ci=fit.conf_int()[1];models.append({'model':name,'status':'estimated','n':len(tum),'df':df,'effect':float(fit.params[1]),'P':float(fit.pvalues[1]),'CI_low':float(ci[0]),'CI_high':float(ci[1]),'condition':float(np.linalg.cond(Z))})
mod=pd.DataFrame(models);mod['FDR_three_models']=stats.false_discovery_control(mod.P.fillna(1).to_numpy());mod.to_csv(O/'results/original31_models.tsv',sep='\t',index=False)
(O/'results/ORIGINAL31_RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print('RESULT',json.dumps(result));print('MODELS',mod.to_string(index=False));print('VERSIONS',pd.crosstab(tum.group_from_filename,tum.program_version).to_string())
# Verify newly downloaded author expected-count file against existing analysis input, without altering either.
p=O/'raw/GSE314812_author_counts.download';old=O.parent/'02_raw/GSE314812_counts.tsv.gz';h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();rec={'new_file':str(p),'old_file':str(old),'new_sha256':h(p),'old_sha256':h(old),'byte_identical':h(p)==h(old)}
if not rec['byte_identical']:
 import gzip,io
 def rd(p):
  z=p.read_bytes();return pd.read_csv(io.BytesIO(gzip.decompress(z) if z[:2]==b'\x1f\x8b' else z),sep='\t')
 fresh=rd(p);prev=rd(old);rec.update(shape_new=list(fresh.shape),shape_old=list(prev.shape),table_identical=fresh.equals(prev));assert rec['table_identical']
else:rec['table_identical']=True
(O/'metadata/NEW_BULK_DOWNLOAD_VERIFICATION.json').write_text(json.dumps(rec,indent=2),encoding='utf-8');print('DOWNLOAD_VERIFY',rec)
