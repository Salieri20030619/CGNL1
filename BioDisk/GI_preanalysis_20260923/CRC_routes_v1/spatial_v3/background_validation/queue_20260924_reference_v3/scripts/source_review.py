from common import *
D=Q/'SOURCE_REVIEW';D.mkdir(exist_ok=True);old=pathlib.Path(C['previous_queue'])/'V02_SOURCE'
x=pd.read_csv(old/'source_expression.tsv',sep='\t');x=x[x.cells>=20];don=x.groupby(['dataset','donor_label','role','label','gene'],as_index=False).agg(mean_CP10K=('mean_CP10K','mean'),detected_pct=('detected_pct','mean'),sample_profiles=('sample','nunique'))
don.to_csv(D/'donor_source_expression.tsv',sep='\t',index=False)
s=don.groupby(['dataset','role','label','gene'],as_index=False).agg(donor_labels=('donor_label','nunique'),median_CP10K=('mean_CP10K','median'),median_detected_pct=('detected_pct','median'),donors_detected=('detected_pct',lambda v:int((v>0).sum())))
s.to_csv(D/'source_summary.tsv',sep='\t',index=False)
cg=s[(s.gene=='CGNL1')&s.role.str.startswith('liver')];cg.to_csv(D/'CGNL1_liver_sources.tsv',sep='\t',index=False)
dump(D/'RESULT.json',{'status':'descriptive_source_review_computed','CGNL1_liver_sources':cg.to_dict('records'),'scope':'>=20 high-confidence lineage cells/profile; donor-equal sample averaging; donor labels provisional; RNA-CNV support not DNA; no unpaired tumor-versus-normal P-values or across-dataset causal comparison'})
