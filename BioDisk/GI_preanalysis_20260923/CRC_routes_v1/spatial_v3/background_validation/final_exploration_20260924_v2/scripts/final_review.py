from common import *
D=Q/'F03_REVIEW';D.mkdir(exist_ok=True);p=pd.read_csv(Q/'F02_SPATIAL/patient_effects.tsv',sep='\t');s=pd.read_csv(Q/'F02_SPATIAL/summary.tsv',sep='\t')
for _,r in s.iterrows():
 z=p[(p.dataset==r.dataset)&(p.gene==r.gene)&np.isclose(p.q,r.q)&np.isclose(p.pseudocount,r.pseudocount)];assert len(z)==r.n_patients and np.isclose(z.effect.mean(),r['mean']) and int((z.effect>0).sum())==r.positive
 if len(z)<10:assert pd.isna(r.P) and pd.isna(r.FDR)
# Confirm reproduction of previous external q25 endpoint, independent of the new sensitivity settings.
old=pd.read_csv(pathlib.Path(C['prior_run'])/'V08_DISTANCE_RANK/patient_effects.tsv',sep='\t');new=p[(p.dataset=='external_author_processed')&np.isclose(p.q,.25)&np.isclose(p.pseudocount,1)].copy();new.patient=new.patient.astype(int);merged=new.merge(old[old.gene.isin(C['final_targets'])],on=['patient','gene'],suffixes=('_new','_old'));assert len(merged)==6 and np.allclose(merged.effect_new,merged.effect_old)
dump(D/'RESULT.json',{'status':'numerical_review_passed','previous_external_endpoint_reproduced':True,'no_new_genes_selected':True,'scope':'Final exploration complete; scientific synthesis and primary literature review remain separate.'})
