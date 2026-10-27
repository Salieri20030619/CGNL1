import pathlib,pandas as pd,numpy as np
B=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');O=B/'deep_v2';r=pd.read_csv(O/'results/bulk_validation/candidate_validation_all_models.tsv',sep='\t');d=pd.read_csv(B/'results/route1/cross_subtype_sensitivity_screen.tsv',sep='\t').set_index('gene');out=[]
for g,z in r.groupby('candidate'):
 m=z[z.model=='main'].iloc[0]
 if m.FDR_fixed126<.1:
  v={'gene':g,'discovery':d.loc[g,'median_delta_log2CPM'],'same_direction':np.sign(d.loc[g,'median_delta_log2CPM'])==np.sign(m.logFC)}
  for _,x in z.iterrows():v[x.model+'_FC']=x.logFC;v[x.model+'_q126']=x.FDR_fixed126;v[x.model+'_qgenome']=x['adj.P.Val']
  out.append(v)
a=pd.DataFrame(out);a.to_csv(O/'results/bulk_validation/cross_cohort_screen.tsv',sep='\t',index=False);print(a.to_string(index=False))
