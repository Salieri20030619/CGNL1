from pathlib import Path
import json,pandas as pd,numpy as np,datetime
O=Path.cwd()/'mechanism_investigation_v2';out=[];cells=[]
for sample in ['PM1','PM2','PM3','PM4']:
 d=pd.read_csv(O/'processed/PMGC2026_scRNA'/(sample+'_cell_source_QC.tsv.gz'),sep='\t');meso=d.group.eq('strict_mesothelial');specific=d.UPK3B+d.LRRN4;support=d.WT1+d.CALB2+d.MSLN;both=(d.UPK3B>0)&(d.LRRN4>0);strong=(specific>=2)&(support>=2);epi=(d.EPCAM>0)&((d[['KRT8','KRT18','KRT19']]>0).sum(axis=1)>=2)
 for scheme in ['primary','sensitivity']:
  keep=d[scheme+'_QC'].astype(bool)
  for label,rule in [('original',meso),('marker_UMI2',meso&strong),('both_specific',meso&both),('UPK3B_only',meso&(d.UPK3B>0)&(d.LRRN4==0)),('LRRN4_only',meso&(d.LRRN4>0)&(d.UPK3B==0))]:
   q=d.loc[rule&keep];n=len(q);tot=float(q.total_UMI.sum());target=float(q.CGNL1.sum());out.append({'sample':sample,'scheme':scheme,'rule':label,'cells':n,'CGNL1_UMI':target,'CGNL1_positive_cells':int((q.CGNL1>0).sum()),'CGNL1_ln1pCP10k':float(np.log1p(target/tot*1e4)) if tot else None,'epithelial_marker_coexpression_cells':int(epi.loc[q.index].sum()),'coverage20':n>=20})
 cells.append(d[meso].assign(marker_UMI2=strong[meso],both_specific=both[meso],epithelial_marker_coexpression=epi[meso]))
r=pd.DataFrame(out);r.to_csv(O/'results/PMGC2026_mesothelial_marker_audit_155.tsv',sep='\t',index=False,encoding='utf-8');pd.concat(cells).to_csv(O/'results/PMGC2026_mesothelial_cells_155.tsv.gz',sep='\t',index=False,encoding='utf-8');record={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'posthoc_identity_audit':True,'same_rules_as_prior_tissue_audit151b':True,'interpretation':'epithelial marker coexpression is not proof of malignancy; mesothelium can express keratins. CGNL1 outcome not used in marker rules. Primary assignments/results retained.','no_P':True};(O/'results/PMGC2026_MARKER_AUDIT_155.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(r.to_string(index=False))
