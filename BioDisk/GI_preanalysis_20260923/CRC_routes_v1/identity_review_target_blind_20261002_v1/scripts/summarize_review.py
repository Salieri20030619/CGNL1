from pathlib import Path
import json,hashlib,datetime,tarfile
import numpy as np,pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1];C=json.loads((R/'SPECIFICATION.json').read_text());ex=set(C['excluded_target_panel']);ex10=ex|set(C['extra_exclusion_sensitivity']);summary={};qa=[]
assert json.loads((R/'PROCESS_EXIT.json').read_text())['exit_code']==0
for ref in ['multiome','normal_liver']:
 t=pd.read_csv(R/'references'/(ref+'_LODO.tsv'),sep='\t');w=t.endo_cells.to_numpy();rec=float(np.nansum(t.endo_recall*w)/w.sum());summary[ref+'_reference_validation']={'donors':len(t),'weighted_accuracy':float(np.average(t.accuracy,weights=t.n)),'median_balanced_accuracy':float(t.balanced_accuracy.median()),'weighted_endothelial_recall':rec}
for ds in ['GSE245552','GSE285991']:
 d=R/'results'/ds;c=json.loads((d/'COMPLETE.json').read_text());assert c['exit_code']==0
 for f,blocked in [('primary_HVG.tsv',ex),('extra_exclusion_HVG.tsv',ex10),('multiome_features.tsv',ex10),('normal_liver_features.tsv',ex10)]:
  genes=set(pd.read_csv(d/f,sep='\t').gene);assert not genes&blocked;qa.append({'dataset':ds,'file':f,'features':len(genes),'excluded_overlap':0,'sha256':hashlib.sha256((d/f).read_bytes()).hexdigest()})
 o=pd.read_csv(d/'cell_identity_evidence.tsv.gz',sep='\t',dtype={'cluster05':str,'cluster10':str});assert len(o)==c['QC_cells'] and o.cell_id.is_unique
 cnt=pd.crosstab(o.old_group,o.reference_agreement);group=[]
 fig,ax=plt.subplots(1,2,figsize=(12,5))
 for a,selected,color,title in [(ax[0],'one_marker_removed','#d94141','Previously removed one-marker cells'),(ax[1],'two_plus_retained','#176b9b','Previously retained two-marker cells')]:
  a.scatter(o.UMAP1,o.UMAP2,s=.4,color='#dddddd',rasterized=True)
  g=o[o.old_group==selected];a.scatter(g.UMAP1,g.UMAP2,s=6,color=color,alpha=.8,rasterized=True)
  a.set_title(title+' (n='+str(len(g))+')');a.set_xticks([]);a.set_yticks([])
 fig.suptitle(ds+' | target-excluded UMAP, candidate overlays');fig.tight_layout();fig.savefig(d/'candidate_overlay.png',dpi=170);plt.close(fig)

 for label,v in o.groupby('old_group'):
  group.append(dict(group=label,cells=len(v),both_immune=int(v.reference_agreement.eq('both_immune').sum()),both_endothelial=int(v.reference_agreement.eq('both_endothelial').sum()),discordant=int(v.reference_agreement.eq('discordant').sum()),both_immune_pct=float(v.reference_agreement.eq('both_immune').mean()*100),both_endothelial_pct=float(v.reference_agreement.eq('both_endothelial').mean()*100)))
 pd.DataFrame(group).to_csv(d/'review_group_summary.tsv',sep='\t',index=False)
 coverage=[]
 for donor,g in o.groupby('donor'):
  both=g.reference_agreement.eq('both_endothelial');coverage.append(dict(donor=donor,QC_cells=len(g),two_refs_endo_cells=int(both.sum()),two_refs_endo_from_old_single=int((both&g.old_group.eq('one_marker_removed')).sum()),two_refs_endo_from_other_QC=int((both&g.old_group.eq('other_QC_cells')).sum()),disease_ref_endo_cells=int(g.multiome_prediction.eq('endothelial').sum()),normal_ref_endo_cells=int(g.normal_liver_prediction.eq('endothelial').sum()),reference_disagreement=int(g.reference_agreement.eq('discordant').sum())))
 pd.DataFrame(coverage).to_csv(d/'patient_reference_endo_coverage.tsv',sep='\t',index=False)
 markers=pd.read_csv(d/'cluster_markers_descriptive.tsv',sep='\t',dtype={'cluster':str});multi=pd.read_csv(d/'multiome_cluster_reference.tsv',sep='\t',dtype={'cluster':str});normal=pd.read_csv(d/'normal_liver_cluster_reference.tsv',sep='\t',dtype={'cluster':str});review=[]
 for cl,g in o.groupby('cluster05'):
  a=multi[multi.cluster==cl].iloc[0];b=normal[normal.cluster==cl].iloc[0];review.append(dict(cluster=cl,cells=len(g),donors=g.donor.nunique(),largest_donor_fraction=float(g.donor.value_counts(normalize=True).max()),old_single_marker=int(g.old_group.eq('one_marker_removed').sum()),old_two_plus=int(g.old_group.eq('two_plus_retained').sum()),disease_reference=a.prediction,normal_reference=b.prediction,disease_score=a.score,normal_score=b.score,top_markers=';'.join(markers[markers.cluster==cl].head(15).gene)))
 pd.DataFrame(review).to_csv(d/'cluster_review_packet.tsv',sep='\t',index=False);summary[ds]={'QC_cells':len(o),'ARI_excluding_original_identity_markers':c['extra_marker_exclusion_ARI05'],'groups':group}
summary['scope']='identity review only; no CGNL1 or junction target association tested; proportions descriptive, not biological replicate tests';summary['updated']=datetime.datetime.now().astimezone().isoformat();(R/'REVIEW_SUMMARY.json').write_text(json.dumps(summary,indent=2));(R/'FEATURE_EXCLUSION_QA.json').write_text(json.dumps(qa,indent=2))
files=[f for sub in ['scripts','results','logs'] for f in (R/sub).rglob('*') if f.is_file() and f.suffix!='.h5ad']+[f for f in (R/'references').iterdir() if f.suffix=='.tsv']+[R/'SPECIFICATION.json',R/'COMPLETE.json',R/'PROCESS_EXIT.json',R/'ENVIRONMENT.json',R/'REVIEW_SUMMARY.json',R/'FEATURE_EXCLUSION_QA.json']
manifest=[{'path':str(f.relative_to(R)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files];(R/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,indent=2));files.append(R/'DELIVERY_MANIFEST.json')
with tarfile.open(R/'identity_review_delivery.tar.gz','w:gz') as t:
 for f in files:t.add(f,arcname=str(f.relative_to(R)))
print(json.dumps(summary),flush=True)
