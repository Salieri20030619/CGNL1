from pathlib import Path
import pandas as pd,json,sys
sys.stdout.reconfigure(encoding='utf-8')
B=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3\background_validation');r=Path(json.loads((B/'LATEST_RUN.json').read_text())['local_run']);s=pd.read_csv(r/'V01/summary.tsv',sep='\t');print(s[s.gene=='CGNL1'].to_json(orient='records',indent=2))
a=pd.read_csv(B.parent/'results/screen/patient_effects.tsv',sep='\t');a=a[(a.scope=='baseline')&(a.gene=='CGNL1')&(a.method=='pseudobulk')&(a.contrast=='tumor_near300_vs_far600')];p=pd.read_csv(r/'V01/patient_effects.tsv',sep='\t');p=p[(p.gene=='CGNL1')&(p.method=='pseudobulk')&(p.rule=='5_of_6')];m=a.merge(p,on='patient',suffixes=('_initial','_strict'));print('MATCHED',json.dumps({'n':len(m),'initial_mean':m.effect_initial.mean(),'strict_mean':m.effect_strict.mean(),'excluded_patients':sorted(set(a.patient)-set(p.patient))}));m.to_csv(r/'V01/matched_initial_strict.tsv',sep='\t',index=False)
