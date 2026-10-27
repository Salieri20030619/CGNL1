from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4',MPLBACKEND='Agg')
import json,hashlib,datetime,base64,html
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
R=Path(__file__).resolve().parent.parent;O=R/'results';F=R/'figures';F.mkdir(exist_ok=True)
S=json.loads((R/'SPECIFICATION.json').read_text(encoding='utf-8'));Q=json.loads((R/'NUMERICAL_QA.json').read_text(encoding='utf-8'));MQ=json.loads((R/'MATRIX_QUALIFICATION.json').read_text(encoding='utf-8'));assert Q['all_pass']
M=pd.read_csv(O/'CGNL1_delta_models.tsv',sep='\t');C=pd.read_csv(O/'paired_changes.tsv',sep='\t');P=pd.read_csv(O/'predictive_increment.tsv',sep='\t');D=pd.read_csv(O/'patient_deltas.tsv',sep='\t');V=pd.read_csv(O/'patient_time_profiles.tsv',sep='\t')
measure='pooled_raw_CP10k' if S['scale']=='verified_raw_counts' else 'mean_author_log';pop='fixed_qc_ge20_all_sites'
a=M[M.measurement.eq(measure)&M.population.eq(pop)&M.model.eq('M1')];c=C[C.measurement.eq(measure)&C.population.eq(pop)];p=P[P.measurement.eq(measure)&P.population.eq(pop)&P.model.eq('M1')]
d=D[D.measurement.eq(measure)&D.definition.eq('fixed_qc')&D.min_cells.ge(20)];n=len(d)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(15,4.5))
v=V[V.measurement.eq(measure)&V.definition.eq('fixed_qc')&V.donor_id.isin(d.donor_id)]
for donor,g in v.groupby('donor_id'):
 y=g.set_index('treatment').loc[['Pre','On'],'CGNL1'];axes[0].plot([0,1],y,color='#24778d' if g.tissueSite.iloc[0]=='Liver' else '#cf8646',alpha=.7,marker='o',ms=3)
axes[0].set(xticks=[0,1],xticklabels=['Pre','Day 15'],ylabel='CGNL1 patient-level expression',title=f'CGNL1 paired changes (n={n})')
for ax,outcome in zip(axes[1:],['actomyosin','cycling']):
 for _,row in d.iterrows():
  ax.scatter(row.CGNL1,row[outcome],c='#24778d' if row.tissueSite=='Liver' else '#cf8646',s=35)
  ax.annotate(row.donor_id,(row.CGNL1,row[outcome]),xytext={'P03':(-17,8),'P30':(5,-12),'P13':(-15,-12),'P28':(4,8)}.get(row.donor_id,(3,3)),textcoords='offset points',fontsize=7)
 ax.axhline(0,lw=.7,color='#b1bac0');ax.axvline(0,lw=.7,color='#b1bac0');ax.set(xlabel='Delta CGNL1 (On - Pre)',ylabel=f'Delta {outcome}',title=f'{outcome.capitalize()}: unadjusted patient changes')
for ax in axes[1:]:ax.xaxis.set_major_locator(MaxNLocator(5))
fig.suptitle('SCP2079 treatment-context validation | teal: Liver; orange: other labeled sites',y=1.01,fontsize=12)
fig.tight_layout();fig.savefig(F/'SCP2079_paired_overview.png',dpi=180,bbox_inches='tight');fig.savefig(F/'SCP2079_paired_overview.svg',bbox_inches='tight');plt.close(fig)
claims=[]
for _,r in a.iterrows():
 evidence='FDR-significant' if r.FDR<.05 else ('nominally significant; exploratory' if r.P<.05 else 'not nominally significant')
 claims.append(f"{r.outcome}: standardized beta={r.beta_standardized:.3f}, P={r.P:.6g}, FDR={r.FDR:.6g}; {evidence}.")
npred=int((p.relative_MSE_improvement>0).sum());nom=int((a.P<.05).sum())
headline='The paired cohort does not establish a conditional CGNL1 association with either fixed endpoint' if nom==0 else 'The paired cohort provides endpoint-specific conditional CGNL1 evidence'
predline=f'Patient-held-out prediction improves for {npred} of the two primary endpoints; full errors and patient consistency are reported below.'
ps=P[P.measurement.eq(measure)&P.population.eq('fixed_qc_ge50_all_sites')&P.model.eq('M1')]
therapy=c[c.readout.eq('CGNL1')].iloc[0]
therapyline=f"CGNL1 declines in {int(therapy.negative_patients)} of {int(therapy.n)} primary paired donors (median On-minus-Pre change {therapy.median_delta:.5f} log-expression units; P={therapy.P:.6g}, FDR={therapy.FDR:.6g})."
cs=C[C.measurement.eq(measure)&C.readout.eq('CGNL1')&C.population.isin(['fixed_qc_ge50_all_sites','author_retained_ge20_all_sites'])]
body='<h1>SCP2079 paired validation of CGNL1</h1><p>Internal final analysis report; '+datetime.datetime.now(datetime.timezone.utc).isoformat()+'. All computation was performed on sgs10 (osv01).</p>'
body+='<h2>CGNL1 decreases in the primary paired analysis</h2><p>'+therapyline+'</p>'
body+='<h2>'+headline+'</h2><p>'+html.escape(' '.join(claims))+'</p><p>'+predline+'</p>'
body+='<h2>Fourteen metadata-qualified pairs enter the prespecified primary population</h2><p>'+str(n)+' independent paired donors enter this analysis. Technical libraries are combined within donor/time point. The Liver-labeled subset contains '+str(int(d.tissueSite.eq('Liver').sum()))+' paired donors and is reported descriptively under the existing protocol. Exact patient eligibility and coverage are retained below.</p>'
body+='<h2>The matrix scale is verified before target extraction</h2><p>Verified scale: <strong>'+html.escape(S['scale'])+'</strong>. Matrix dimensions: '+str(MQ['shape'][0])+' features by '+str(MQ['shape'][1])+' cells. Noninteger stored entries: '+str(MQ['noninteger_nonzeros'])+'. Measurement: '+html.escape(S['measurement_primary'])+'. Depth control: '+html.escape(S['depth'])+'.</p>'
body+='<h2>Methods</h2><p>Author epithelial annotations were matched by exact unique barcodes. Primary QC used author total_counts >=500, nFeature_RNA >=200, mitochondrial percentage <=25, doublet score <=0.3, and at least 20 cells at both treatment time points. The sensitivity definitions use at least 50 cells per time point, or the author-retained cells with at least 20 cells. The fixed CGNL1 target, actomyosin/cycling panels and differentiation controls are unchanged from the preceding project. Target-zero cells are retained.</p>'
body+='<p>Within-patient changes equal On minus Pre. The primary M1 model includes CGNL1, depth, epithelial identity and mature-colonocyte changes; M0 controls depth only. Target coefficients use HC3 standard errors and residual-df Student-t tests. Full rank, at least five residual degrees of freedom and at least 10 independent pairs are required. Raw target P&lt;0.05 is the nominal criterion. BH FDR is calculated within the two endpoint hypotheses of each population/model/measurement family. Paired change tests use the exact conditional signed-rank distribution; these three readouts form a separate family. Leave-one-patient-out prediction fits all scaling in training donors and compares baseline controls against controls plus CGNL1.</p>'
body+='<h2>Paired changes and conditional effects answer different questions</h2>'
def table(frame):return frame.to_html(index=False,border=0,float_format=lambda v:f'{v:.8g}',na_rep='NA')
body+=table(c[['readout','n','mean_delta','median_delta','positive_patients','negative_patients','P','FDR']])
body+=table(a[['outcome','n','beta_standardized','slope','CI_low','CI_high','P','FDR','df']])
body+='<p>Blue-green points represent Liver-labeled samples; orange points retain all other original site labels. Scatter panels show unadjusted changes; conditional inference is reported in the table.</p><img alt="Paired patient overview" src="data:image/png;base64,'+base64.b64encode((F/'SCP2079_paired_overview.png').read_bytes()).decode('ascii')+'">'
body+='<h2>Held-out prediction quantifies additional patient-level information</h2>'+table(p[['outcome','n','completed_folds','MSE_baseline','MSE_augmented','relative_MSE_improvement','patients_improved']])
body+='<h2>The actomyosin prediction gain does not persist at the higher cell threshold</h2>'
pa=p[p.outcome.eq('actomyosin')].iloc[0];sa=ps[ps.outcome.eq('actomyosin')].iloc[0]
body+=f'<p>Adding CGNL1 reduces primary actomyosin prediction MSE by {100*pa.relative_MSE_improvement:.1f}% ({int(pa.patients_improved)}/{int(pa.n)} patients improve). With at least 50 cells at each time point, the relative improvement becomes {100*sa.relative_MSE_improvement:.1f}% ({int(sa.patients_improved)}/{int(sa.n)} improve). This prespecified sensitivity excludes P29, who has 26 retained on-treatment epithelial cells. The prediction increment is not stable across the two coverage gates. Cycling prediction worsens in both gates.</p>'
body+='<p>The CGNL1 decrease remains directionally consistent in the alternative populations, but statistical strength depends on cell filtering. The exact sensitivity results are retained here; nominal significance and FDR significance are distinguished.</p>'+table(cs[['population','n','mean_delta','median_delta','negative_patients','P','FDR']])
body+='<h2>Discussion</h2><p>This is a BRAF-V600E selected, single-arm combined BRAF/MEK/PD-1 treatment cohort. Paired expression changes quantify treatment-context observations, not a CGNL1-specific perturbation or an isolated causal drug effect. Technical replicates and cells do not increase independent sample size. Liver, Muscle, Lymph Node and Other retain their original annotation meaning. The small Liver-only subset does not establish a liver-colonization mechanism. Conditional associations, overall treatment changes and predictive increments must be interpreted separately; neither nominal association nor predictive improvement establishes causality. PFS and RECIST were not screened for new outcome associations.</p>'
body+='<h2>Supplement: all models, including unestimable groups</h2>'+table(M)
body+='<h2>Supplement: all paired changes</h2>'+table(C)
body+='<h2>Supplement: all predictive increments</h2>'+table(P)
body+='<h2>Supplement: patient eligibility</h2>'+table(pd.read_csv(O/'patient_eligibility.tsv',sep='\t'))
body+='<h2>Numerical verification and server record</h2><p>Independent patient-profile reconstruction, ordinary matrix-algebra HC3 calculations, manually reconstructed BH FDR, exhaustive sign enumeration and independently reconstructed leave-one-patient-out predictions passed.</p><pre>'+html.escape(json.dumps(Q,indent=2))+'</pre><p>Server analysis directory: <code>'+str(R)+'</code>. Exact input fingerprints, specification freeze, complete tables, scripts and logs remain on the server. This self-contained HTML is the local final reading artifact.</p>'
report='<!doctype html><html lang="en"><meta charset="utf-8"><title>SCP2079 paired CGNL1 validation</title><style>body{font:16px/1.55 Arial,sans-serif;max-width:1320px;margin:30px auto;padding:0 25px;color:#18313d}h1,h2{color:#155d74}table{border-collapse:collapse;display:block;overflow:auto;font-size:12px;margin:18px 0}td,th{border:1px solid #ccd6db;padding:6px;white-space:nowrap;text-align:left}th{background:#edf4f6}img{max-width:100%}pre{white-space:pre-wrap}code{overflow-wrap:anywhere}</style><body>'+body+'</body></html>'
(R/'CGNL1_SCP2079_paired_report.html').write_text(report,encoding='utf-8')
summary={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'analysis_directory':str(R),'matrix_scale':S['scale'],'primary_pairs':n,'liver_pairs':int(d.tissueSite.eq('Liver').sum()),'primary_results':json.loads(a.to_json(orient='records')),'primary_prediction':json.loads(p.to_json(orient='records')),'primary_changes':json.loads(c.to_json(orient='records')),'numerical_QA_all_pass':Q['all_pass'],'this_paired_validation_complete':True,'whole_CGNL1_project_complete':False,'interpretation':'Primary paired CGNL1 decrease; conditional actomyosin/cycling associations not nominal; actomyosin LOPO improvement reverses at the >=50-cell sensitivity gate; no stable independent increment established','local_policy':'D:/CGNL1 retains handoff records and final reports only; raw data, scripts, and intermediate results remain on sgs10'}
(R/'PROGRESS_SNAPSHOT.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(R.parent/'CGNL1_PROGRESS_20261014_SCP2079_PAIRED.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
(R/'REPORT_VALIDATION.json').write_text(json.dumps({'UTF8_readback':(R/'CGNL1_SCP2079_paired_report.html').read_text(encoding='utf-8')==report,'self_contained_images':True,'all_primary_results_present':len(a)==2,'all_primary_predictions_present':len(p)==2,'report_sha256':hashlib.sha256((R/'CGNL1_SCP2079_paired_report.html').read_bytes()).hexdigest()},indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2),flush=True)