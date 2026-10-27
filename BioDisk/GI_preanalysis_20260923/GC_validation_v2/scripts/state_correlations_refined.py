import pathlib,json,pandas as pd,numpy as np
from scipy import stats
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/GC_validation_v2');O=R/'04_results'/'independent_scRNA_refined'
x=pd.read_csv(O/'marker_expression_patient_compartment.tsv',sep='\t');x=x[(x.celltype=='Macrophage_enriched')&(x.n_cells>=20)];w=x.pivot(index='patient',columns='gene',values='log2CPM');partners=['APOE','LIPA','PLIN2','GPNMB','MRC1','SPP1'];res=[]
(R/'00_docs'/'state_correlation_plan_refined.json').write_text(json.dumps({'anchors':['CD36','TREM2'],'partners':partners,'family_size':12,'interpretation':'exploratory within-compartment patient correlations; not causal, no candidate selection after P values'},indent=2))
for a in ['CD36','TREM2']:
 for b in partners:
  z=w[[a,b]].dropna();rho,p=stats.spearmanr(z[a],z[b]);res.append({'anchor':a,'partner':b,'patients':len(z),'rho':rho,'P':p})
r=pd.DataFrame(res);r['FDR_BH']=stats.false_discovery_control(r.P.values,method='bh');r.to_csv(O/'exploratory_state_correlations.tsv',sep='\t',index=False);print(r.to_string(index=False))
