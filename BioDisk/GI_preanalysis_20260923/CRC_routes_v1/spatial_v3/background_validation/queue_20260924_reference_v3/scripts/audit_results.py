from common import *
from scipy import stats
records=[]
for name in ['V03_MIXTURE','V03_HEPREF_SENS']:
 D=Q/name;sec=pd.read_csv(D/'section_effects.tsv',sep='\t');pat=pd.read_csv(D/'patient_effects.tsv',sep='\t');summ=pd.read_csv(D/'summary.tsv',sep='\t');reb=sec.groupby(['patient','gene','mask']).effect.mean().reset_index();m=pat.merge(reb,on=['patient','gene','mask'],suffixes=('_saved','_rebuilt'));assert len(m)==len(pat) and np.allclose(m.effect_saved,m.effect_rebuilt)
 for _,r in summ.iterrows():
  v=pat[(pat.gene==r.gene)&(pat['mask']==r['mask'])];x=v.effect.to_numpy();assert len(x)==r.n_patients;assert np.isclose(x.mean(),r.mean_effect);assert int((x>0).sum())==r.positive
  if len(x)>=10:
   p=stats.ttest_1samp(x,0).pvalue;se=stats.sem(x);t=stats.t.ppf(.975,len(x)-1);assert np.isclose(p,r.P) and np.allclose([x.mean()-t*se,x.mean()+t*se],[r.CI_low,r.CI_high])
  else:assert pd.isna(r.P) and pd.isna(r.FDR)
 for mask,z in summ.groupby('mask'):
  pp=np.nan_to_num(z.P.to_numpy(),nan=1);order=np.argsort(pp);adj=np.minimum.accumulate((pp[order]*4/np.arange(1,len(pp)+1))[::-1])[::-1];qq=np.empty(len(pp));qq[order]=np.minimum(adj,1);valid=z.P.notna().to_numpy();assert np.allclose(qq[valid],z.FDR.to_numpy()[valid])
 records.append({'analysis':name,'section_to_patient_means':'passed','direction_n_effect_CI_P_FDR':'passed','small_n_no_inference':'passed'})
exitstatus=json.loads((Q/'QUEUE_EXIT.json').read_text());assert exitstatus['exit_code']==0
for stage in ['NORMAL_REF','V02_SOURCE','V03_MIXTURE','V03_HEPREF_SENS','V04_SPATIAL','V04_DIAGNOSTIC','V04_BROAD_EP_SENS']:
 e=json.loads((Q/stage/'EXIT.json').read_text());assert e['exit_code']==0 and (Q/stage/'RESULT.json').exists()
dump(Q/'NUMERICAL_AUDIT.json',{'status':'passed','checks':records,'queue_exit':exitstatus,'interpretation':'Model execution success is distinct from region eligibility and biological association.'})
print('Numerical and completion audit passed')
