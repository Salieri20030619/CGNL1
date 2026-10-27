from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,datetime,hashlib
import numpy as np,pandas as pd
from scipy import stats
O=Path.cwd()/'mechanism_investigation_v2';p=O.parent/'results/bulk_patient_deltas.tsv';panels=['epithelial','stroma','myeloid','adipocyte','endothelial','mesothelial']
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'posthoc':True,'prior_known':'17pairedCGNL1PMminusPTresponse andpriorfullproxyregression; sourcecorrelationsnotyetcomputed','aim':'compareall6fixedtargetexcludedsourceproxieswithoutassumingmesothelialorigin','primary':'Spearmancorrelation ofwithinpatientCGNL1changeversuswithinpatientsourceproxychange; all17patients','sensitivity':'partialrankcorrelation afterlinearregressionofeachrankvectoragainstinterceptandpreservation_delta;df=n-3','families':'6primarySpearman and6partialrank tests correctedseparately','robustness':'all17leave-one-donor-out correlations reported asrange, notselected','limits':['RNAproxyassociationnotcellfractionorcausalfunction','same17patientsnotnewvalidationcohort','correlatedlineagemarkersmaylimitdiscrimination','primaryPusesstandardasymptoticrankapproximation; partialPusespartialcorrelationtapproximation']};sp=O/'BULK_SOURCE_CONTEXT_SPECIFICATION.json';assert not sp.exists();sp.write_text(json.dumps(spec,indent=2),encoding='utf-8')
d=pd.read_csv(p,sep='\t');assert len(d)==17 and d.patient.is_unique and not d[['CGNL1','preservation_delta']+panels].isna().any().any();Z=np.c_[np.ones(len(d)),d.preservation_delta.to_numpy()];x=stats.rankdata(d.CGNL1);xr=x-Z@np.linalg.lstsq(Z,x,rcond=None)[0];out=[]
for panel in panels:
 r=stats.spearmanr(d.CGNL1,d[panel]);loo=[stats.spearmanr(d.drop(i).CGNL1,d.drop(i)[panel]).statistic for i in d.index];y=stats.rankdata(d[panel]);yr=y-Z@np.linalg.lstsq(Z,y,rcond=None)[0];rho=np.corrcoef(xr,yr)[0,1];df=len(d)-3;t=rho*np.sqrt(df/(1-rho*rho));P=2*stats.t.sf(abs(t),df)
 out.extend([{'panel':panel,'model':'Spearman','n':17,'rho':r.statistic,'P':r.pvalue,'LOO_rho_min':min(loo),'LOO_rho_max':max(loo)},{'panel':panel,'model':'partial_rank_preservation','n':17,'rho':rho,'P':P,'df':df}])
r=pd.DataFrame(out);r['FDR']=np.nan
for model,ix in r.groupby('model').groups.items():r.loc[ix,'FDR']=stats.false_discovery_control(r.loc[ix,'P'])
r.to_csv(O/'results/GSE314812_source_context_correlations.tsv',sep='\t',index=False,encoding='utf-8');print(r.to_string(index=False))
(O/'results/BULK_SOURCE_CONTEXT_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'spec_sha256':hashlib.sha256(sp.read_bytes()).hexdigest(),'n_donors':17,'tests':len(r),'causal_claim':False},indent=2),encoding='utf-8')
