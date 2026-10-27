from pathlib import Path
import numpy as np,pandas as pd,json,datetime
from scipy import stats
from statsmodels.stats.multitest import multipletests
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/epithelial_independence_20261012_v1');O=R/'results'
def adjust(t,group,pcol,outcol):
 t[outcol]=np.nan
 for _,ix in t.groupby(group,dropna=False).groups.items():
  ok=[i for i in ix if np.isfinite(t.loc[i,pcol])]
  if ok:t.loc[ok,outcol]=multipletests(t.loc[ok,pcol],method='fdr_bh')[1]
 return t
# Paired disease-site profiles and fixed marker-relative contrasts.
bp=pd.read_csv(O/'B_all_group_profiles.tsv',sep='\t')
metrics=['CGNL1','CGNL1_minus_epithelial','CGNL1_minus_mature']
bp[metrics[1]]=bp.CGNL1-bp.epithelial_identity;bp[metrics[2]]=bp.CGNL1-bp.mature_colonocyte
pair=[]
for group in ['Tumor_no_Tu10','CNV95_noTu10','CNV99_noTu10']:
 v=bp[bp.celltype.eq(group)&bp.n_cells.ge(20)]
 for patient,g in v.groupby('patient'):
  if not {'CCT','LCT'}.issubset(set(g.site)):continue
  a=g[g.site.eq('CCT')].iloc[0];b=g[g.site.eq('LCT')].iloc[0]
  for metric in metrics:pair.append(dict(group=group,patient=patient,metric=metric,primary=a[metric],liver=b[metric],delta=b[metric]-a[metric],n_shared_states=np.nan))
for patient,v in bp[bp.celltype.str.match(r'^Tu[0-9]{2}_')&~bp.celltype.str.startswith('Tu10')&bp.n_cells.ge(20)].groupby('patient'):
 shared=set(v[v.site.eq('CCT')].celltype)&set(v[v.site.eq('LCT')].celltype)
 if len(shared)<2:continue
 for metric in metrics:
  a=v[v.site.eq('CCT')&v.celltype.isin(shared)][metric].mean();b=v[v.site.eq('LCT')&v.celltype.isin(shared)][metric].mean()
  pair.append(dict(group='shared_author_subtypes',patient=patient,metric=metric,primary=a,liver=b,delta=b-a,n_shared_states=len(shared)))
pa=pd.DataFrame(pair);pa.to_csv(O/'B_paired_effects.tsv',sep='\t',index=False)
summ=[]
for (group,metric),g in pa.groupby(['group','metric']):
 d=g.delta.to_numpy();d=d[np.isfinite(d)];nz=d[np.abs(d)>1e-12]
 summ.append(dict(group=group,metric=metric,n=len(d),positive=int((d>0).sum()),negative=int((d<0).sum()),mean_delta=float(np.mean(d)),median_delta=float(np.median(d)),P=stats.binomtest(int((nz>0).sum()),len(nz),.5).pvalue if len(nz)>=2 else np.nan,family='B_three_contrasts_'+group,evidence='small_sample_exploratory_site_contrast'))
summ=adjust(pd.DataFrame(summ),['family'],'P','FDR');summ.to_csv(O/'B_paired_summary.tsv',sep='\t',index=False)

(O/'B_SUBTYPE_QA_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'implementation_correction':'Only TuNN_ author subtype rows; aggregate Tumor_ profiles excluded from shared-subtype averaging','primary_contrasts_unchanged':True},indent=2))
print('B_SUBTYPE_QA_COMPLETE')
