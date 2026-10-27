import pathlib,json,numpy as np,pandas as pd
from scipy import stats
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');V=B/'GC_validation_v2';O=V/'candidate_routes';plan=json.loads((O/'PLAN.json').read_text());(O/'BULK_ROUTE_B_PLAN.json').write_text(json.dumps({'scope':'paired gene changes versus four non-overlapping composition proxy changes; 3x4=12 correlations; descriptive source association, not measured proportions','family':'12 bulk paired-change correlations; adjusted-for-preservation-difference sensitivity as separate12 family'},indent=2))
m=pd.read_csv(V/'candidate_discovery/anatomy_proxy_sample_scores.tsv',sep='\t');x=pd.read_csv(B/'03_processed/GC314812_TMM_logCPM.tsv.gz',sep='\t');genes=plan['genes'];proxies=['adipocyte','myeloid','endothelial','mesothelial']
for g in genes:
 z=x[x.symbol==g];assert len(z)==1;m[g]=z[m['sample']].to_numpy().ravel()
rows=[]
for pt,z in m.groupby('patient'):
 p=z[z.site=='GCPT'].iloc[0];a=z[z.site=='GCPM'].iloc[0];r={'patient':pt,'preservation_difference':int(a.preservation=='FFPE')-int(p.preservation=='FFPE')}
 for g in genes+proxies:r[g]=a[g]-p[g]
 rows.append(r)
d=pd.DataFrame(rows);assert len(d)==17;d.to_csv(O/'bulk_gene_proxy_paired_changes.tsv',sep='\t',index=False);out=[]
for g in genes:
 for q in proxies:
  co=stats.spearmanr(d[g],d[q]);C=np.column_stack([np.ones(len(d)),d.preservation_difference]);rank=np.linalg.matrix_rank(C);xx=stats.rankdata(d[g]);yy=stats.rankdata(d[q]);xx-=C@np.linalg.lstsq(C,xx,rcond=None)[0];yy-=C@np.linalg.lstsq(C,yy,rcond=None)[0];pr=np.corrcoef(xx,yy)[0,1];df=len(d)-rank-1;pp=2*stats.t.sf(abs(pr)*np.sqrt(df/max(1-pr*pr,1e-15)),df);out.append({'gene':g,'proxy':q,'n':17,'rho':co.statistic,'P':co.pvalue,'partial_rho_preservation':pr,'P_partial':pp})
r=pd.DataFrame(out);r['FDR_12']=stats.false_discovery_control(r.P.values,method='bh');r['FDR_partial_12']=stats.false_discovery_control(r.P_partial.values,method='bh');r.to_csv(O/'bulk_gene_proxy_correlations.tsv',sep='\t',index=False);print(r.to_string(index=False));(O/'BULK_ROUTE_B_COMPLETE.json').write_text(json.dumps({'exit_code':0,'patients':17,'hypotheses':12},indent=2))
