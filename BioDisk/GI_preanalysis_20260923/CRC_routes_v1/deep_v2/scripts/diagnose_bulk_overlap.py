import pathlib,pandas as pd,numpy as np,json
from sklearn.linear_model import LinearRegression
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2');R=O/'results/bulk_validation';m=pd.read_csv(R/'sample_values.tsv',sep='\t');y=(m.site=='LM').astype(float).values;P=pd.get_dummies(m.patient.astype(str),drop_first=True,dtype=float).values;C=np.column_stack([np.ones(len(m)),P]);yr=y-C@np.linalg.lstsq(C,y,rcond=None)[0];out=[]
for name,cols in [('hepatic',['hepatic']),('composition',['stromal','epithelial']),('combined',['hepatic','stromal','epithelial'])]:
 A=m[cols].values;Ar=A-C@np.linalg.lstsq(C,A,rcond=None)[0];pred=Ar@np.linalg.lstsq(Ar,yr,rcond=None)[0];r2=1-np.sum((yr-pred)**2)/np.sum(yr**2);out.append({'model':name,'within_patient_site_R2_from_proxies':float(r2),'site_VIF':float(1/(1-r2)),'remaining_site_variance_fraction':float(1-r2)})
pd.DataFrame(out).to_csv(R/'site_overlap_diagnostics.tsv',sep='\t',index=False);print(pd.DataFrame(out).to_string(index=False));print(m.groupby('site').hepatic.agg(['min','max','median']).to_string())
