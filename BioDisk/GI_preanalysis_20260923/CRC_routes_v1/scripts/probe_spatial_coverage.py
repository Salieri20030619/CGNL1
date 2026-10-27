import pathlib,json,numpy as np
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');a=np.load(O/'results/route2/region_aggregates.npz',allow_pickle=True);sy=a['genes'];print('GENES',len(sy),sy[:12]);p=json.loads((O/'metadata/ANALYSIS_PLAN.json').read_text())
for k,gs in p['programs'].items():print(k,'found',[g for g in gs if g in sy],'missing',[g for g in gs if g not in sy])
