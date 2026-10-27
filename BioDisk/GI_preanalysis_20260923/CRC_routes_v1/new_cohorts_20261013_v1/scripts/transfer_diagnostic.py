from pathlib import Path
import pandas as pd,numpy as np,json,datetime
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/new_cohorts_20261013_v1');O=R/'results'
p=pd.read_csv(O/'N3_transfer_predictions.tsv',sep='\t');rows=[]
for key,g in p.groupby(['train','test','outcome']):
 eb=g.baseline-g.observed;ea=g.augmented-g.observed
 mb=float((eb**2).mean());ma=float((ea**2).mean());bb=float(eb.mean());ba=float(ea.mean())
 vb=float(eb.var(ddof=0));va=float(ea.var(ddof=0))
 assert abs(mb-(bb*bb+vb))<1e-10 and abs(ma-(ba*ba+va))<1e-10
 rows.append(dict(zip(['train','test','outcome'],key),n=len(g),MSE_baseline=mb,MSE_augmented=ma,improvement=1-ma/mb,bias_baseline=bb,bias_augmented=ba,centered_MSE_baseline=vb,centered_MSE_augmented=va,centered_improvement=1-va/vb,bias_squared_reduction=bb*bb-ba*ba,centered_error_reduction=vb-va))
pd.DataFrame(rows).to_csv(O/'N3_transfer_error_decomposition.tsv',sep='\t',index=False)
(R/'TRANSFER_DECOMPOSITION_NOTE.json').write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),role='post-result descriptive diagnostic of already frozen transfer predictions',identity='MSE = squared mean error + mean squared centered error',heldout_refitting=False,calibrated_prediction_claim=False,no_new_P_values=True),indent=2))
print(pd.DataFrame(rows).to_string(index=False))
