from pathlib import Path
import json,csv,itertools
import numpy as np
from scipy.stats import rankdata
b=Path(__file__).resolve().parents[1]
pdo=list(csv.DictReader((b/'results/independent_pdo_326/patient_deltas.tsv').open(),delimiter='\t'));pdx=list(csv.DictReader((b/'results/secondary_pdx_329/model_deltas.tsv').open(),delimiter='\t'));out=[];detail=[]
for method in ['median_ratio','human_CPM']:
 a={r['case']:float(r['delta']) for r in pdo if r['gene']=='CGNL1' and r['method']==method};z={r['case']:float(r['delta']) for r in pdx if r['gene']=='CGNL1' and r['method']==method and r['phase']=='PDX_72h'};keys=sorted(set(a)&set(z));d=np.array([z[k]-a[k] for k in keys]);rank=rankdata(abs(d));obs=abs(sum(rank*np.sign(d)));p=sum(abs(sum(rank*np.array(s)))>=obs-1e-12 for s in itertools.product([-1,1],repeat=len(d)))/2**len(d)
 out.append({'method':method,'n':len(keys),'mean_interaction':float(d.mean()),'median_interaction':float(np.median(d)),'positive':int(sum(d>0)),'P_exploratory':p,'PDO_matched_mean':float(np.mean([a[k] for k in keys])),'PDX_matched_mean':float(np.mean([z[k] for k in keys]))});detail.extend({'case':k,'method':method,'PDO':a[k],'PDX':z[k],'interaction':z[k]-a[k]} for k in keys)
r={'results':out,'patient_interactions':detail,'status':'exploratory after unblinding marginal results; not independent mainline validation'};(b/'results/CONTEXT_INTERACTION_331.json').write_text(json.dumps(r,indent=2),encoding='utf-8');print(json.dumps(out))
