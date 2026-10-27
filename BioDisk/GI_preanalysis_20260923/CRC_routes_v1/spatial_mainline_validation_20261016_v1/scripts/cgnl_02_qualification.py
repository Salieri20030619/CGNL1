from pathlib import Path
from datetime import datetime,timezone
import json,hashlib
import pandas as pd,numpy as np
from scipy.spatial import cKDTree
R=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');O=Path.cwd()
p=R/'spatial_v3/results/screen/spot_target_coordinates.tsv.gz';d=pd.read_csv(p,sep='\t')
print('SPOT_COLUMNS',list(d)); print('REGIONS',d.Layer3.value_counts().to_dict())
jobs=json.loads((R/'spatial_v3/metadata/LIVER_DOWNLOAD_PLAN.json').read_text());counts=[]
for j in jobs:
 v=d[d['sample'].eq(j['sample'])];good=v.QC_pass.astype(str).str.lower().isin(['true','1']);v=v[good]
 t=v[v.Layer3.eq('Liver macrometastasis tumor')];dist=t.distance_to_hepatic_um.to_numpy();hep=v[v.Layer3.eq('Hepatic lobule')]
 row={**j,'QC_spots':len(v),'macro_spots':len(t),'hepatic_spots':len(hep),'macro_50_200':int(((dist>=50)&(dist<=200)).sum()),'macro_200_500':int(((dist>200)&(dist<500)).sum()),'macro_500_1500':int(((dist>=500)&(dist<=1500)).sum()),'macro_200_400':int(((dist>200)&(dist<=400)).sum()),'macro_600_1500':int(((dist>=600)&(dist<=1500)).sum()),'recurrence':'Recurrence' in j['sample_type']}
 counts.append(row)
q=pd.DataFrame(counts);q.to_csv(O/'metadata/interface_geometry_qualification.tsv',sep='\t',index=False)
print(q.drop(columns=['url','sample_type']).to_string(index=False))
print('Primary eligible baseline patients',q[(~q.recurrence)&(q.macro_50_200>=20)&(q.macro_500_1500>=20)].patient.unique().tolist())
# Examine input schemas and panels without calculating any new target effect.
print('PANELS',json.loads((R/'v3_pilot_20261006/plan_v3.json').read_text())['panels'])
metadata={'utc':datetime.now(timezone.utc).isoformat(),'source':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'rows':len(d),'qualification_uses_target_expression':False}
(O/'metadata/GEOMETRY_QUALIFICATION.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
