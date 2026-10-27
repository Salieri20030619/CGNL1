from common import *
v=pd.read_csv(O/'D4_patient_profiles.tsv',sep='\t')
rows=[];elig=[]
for definition,g in v.groupby('definition'):
 for gene in S['genes']:
  f=fitted(g,S['models']['M1'],gene)
  if f:rows.append(dict(definition=definition,gene=gene,family='postprimary_D4_components7_'+definition,**f[0]))
  else:elig.append(dict(definition=definition,gene=gene,status='not_estimable'))
t=pd.DataFrame(rows);t['FDR']=np.nan
for _,ix in t.groupby('family').groups.items():t.loc[ix,'FDR']=bh(t.loc[ix,'P'],7)
t.to_csv(O/'D4_postprimary_all7_components.tsv',sep='\t',index=False)
print(t[['definition','gene','beta','P','FDR']].to_json(orient='records'))
