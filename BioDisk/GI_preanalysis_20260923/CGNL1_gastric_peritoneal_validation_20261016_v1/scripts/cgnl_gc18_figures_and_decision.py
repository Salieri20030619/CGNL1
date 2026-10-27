from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4')
import json,datetime
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path.cwd();m=pd.read_csv(O/'processed/original31_sample_values.tsv',sep='\t');d=pd.read_csv(O/'results/bulk_patient_deltas.tsv',sep='\t');p=pd.read_csv(O/'results/identity_source_profiles.tsv',sep='\t');summary=pd.read_csv(O/'results/identity_source_summary.tsv',sep='\t');c=pd.read_csv(O/'results/identity_source_contrasts.tsv',sep='\t')
group=m.groupby('group_from_filename').CGNL1.agg(['count','mean','median','min','max']).reset_index();group.to_csv(O/'results/original31_group_summary.tsv',sep='\t',index=False);print('GROUPS',group.to_string(index=False))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,3,figsize=(14,5),layout='constrained');ax=axes[0];rng=np.random.default_rng(20261016);groups=['without_PC','with_PC'];colors=['#507e9b','#c47a39']
for i,g in enumerate(groups):
 vals=m[m.group_from_filename==g].CGNL1.to_numpy();ax.scatter(i+rng.uniform(-.12,.12,len(vals)),vals,c=colors[i],s=28,alpha=.85);ax.plot([i-.18,i+.18],[np.mean(vals)]*2,c='black',lw=2)
ax.set_xticks([0,1],['Without PC\nn=19','With PC\nn=12']);ax.set_ylabel('CGNL1, vendor log2 scaled-RPM');ax.set_title('Published cohort: actual CHP files');ax.text(.03,.02,'Difference +2.026; FC 4.07\nWelch P=0.00609\nFilename groups; folder-label conflict flagged',transform=ax.transAxes,fontsize=9)
ax=axes[1]
for _,r in d.iterrows():ax.plot([0,1],[r.primary,r.metastasis],color='#1b8290' if r.CGNL1>0 else '#b16b53',alpha=.65,lw=1,marker='o',ms=3)
ax.set_xticks([0,1],['Primary','Peritoneal lesion']);ax.set_ylabel('CGNL1, TMM log2CPM');ax.set_title('Separate cohort: 17 patient pairs');ax.text(.03,.02,'11/17 increases; median +0.518\nPaired signed-rank P=0.07968\nPreservation-adjusted P=0.02181',transform=ax.transAxes,fontsize=9)
ax=axes[2];z=p[p.n_cells>=20];pair=z[z.group.isin(['strict_mesothelial','conservative_epithelial'])].pivot(index='patient',columns='group',values='CGNL1').dropna()
for _,r in pair.iterrows():ax.plot([0,1],[r.conservative_epithelial,r.strict_mesothelial],color='#1b8290',alpha=.6,marker='o',ms=3)
ax.set_xticks([0,1],['Epithelial-like','Strict mesothelial-like']);ax.set_ylabel('CGNL1, pooled ln(1+CP10k)');ax.set_title('Within-patient source contrast');ax.text(.03,.94,'10/10 higher in mesothelial-like\nP=0.00195; FDR=0.00195',transform=ax.transAxes,va='top',fontsize=9)
fig.suptitle('Gastric peritoneal dissemination: complementary evidence, distinct estimands',fontsize=14);fig.savefig(O/'figures/01_core_evidence.png',dpi=170);fig.savefig(O/'figures/01_core_evidence.svg');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained');ax=axes[0];order=['conservative_epithelial','strict_mesothelial','strict_mesothelial_low_cancer_proxy','macrophage'];labels=['Epithelial-like','Strict mesothelial-like','Strict mesothelial-like\nlow cancer proxy','Macrophage-enriched']
for i,g in enumerate(order):
 vals=z[z.group==g].CGNL1.to_numpy();ax.scatter(i+rng.uniform(-.12,.12,len(vals)),vals,s=23,c='#1b8290',alpha=.65);ax.plot([i-.18,i+.18],[np.median(vals)]*2,c='black',lw=2);ax.text(i,-.08,'n='+str(len(vals)),ha='center',transform=ax.get_xaxis_transform(),fontsize=9)
ax.set_xticks(range(4),labels,rotation=15,ha='right');ax.set_ylabel('Patient-level pooled ln(1+CP10k)');ax.set_title('Target-excluded source definitions; >=20 cells')
ax=axes[1];s=pd.read_csv(O/'results/mesothelial_state_models.tsv',sep='\t');q=s[s.model=='HC3_depth_identity'];ax.barh(q.panel,-np.log10(q.P),color='#a2bac4');ax.axvline(-np.log10(.05),color='#a3673f',ls='--',label='Raw P=0.05');ax.set_xlabel('-log10(P), CGNL1 effect adjusted for depth and identity');ax.set_title('No established independent state mechanism');ax.legend(frameon=False)
fig.savefig(O/'figures/02_source_and_limits.png',dpi=170);fig.savefig(O/'figures/02_source_and_limits.svg');plt.close(fig)
decision={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'candidate_worth_next_validation_not_established_mechanism','question_zh':'胃癌腹膜播散相关的CGNL1组织信号，是否对应可区分于细胞数量变化的宿主间皮/基质状态，并与腹膜播散易感或病灶形成有关？','go':['actualpublishedcohortdirectionreproduced','separatepairedlesioncohortpositiveeffectwithnominalpreservationadjustedassociation','targetexcludedmesothelialsourceenrichmentreplicatedacrosspatients'],'not_established':['cancer-cell-intrinsicdriver','mesothelialbarrierlossormechanicalstate','independentjunction/mesenchymal/cyclingassociation','mesothelialCGNL1differencebyverifiedpatientPCstatus','clinicalpredictionvalue'], 'first_next_gate':'patient-levelPC/cytologymappingandindependenttissue/mesothelialsourcevalidation, rather than generic pathway expansion','material_limitations':['original31folder/filenameconflict','CHP effect/P differs from publishedfigure','original31nominalgenomewideFDRnotpassed','17pairedprimarysignedranknotnominal','publicGSE239676targetabsent','GSE228598patientclinicalmappingmissing'],'whole_project_complete':False,'this_actual_data_validation_complete':True,'user_review_pending':True};(O/'DECISION.json').write_text(json.dumps(decision,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(decision,ensure_ascii=False))
