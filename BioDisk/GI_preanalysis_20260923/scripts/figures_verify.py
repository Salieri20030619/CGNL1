import pathlib,json,gzip,hashlib,pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import t as student_t
R=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');O=R/'04_results';F=O/'figures';G=O/'GC_bulk'
def fig(n):plt.tight_layout();plt.savefig(F/(n+'.png'),dpi=180,bbox_inches='tight');plt.savefig(F/(n+'.pdf'),bbox_inches='tight');plt.close()
a=pd.read_csv(G/'diagnostic_panel_statistics.tsv',sep='\t');b=pd.read_csv(G/'composition_sensitivity_panels.tsv',sep='\t');z=a.merge(b,on='panel',suffixes=('','_composition'))
y=np.arange(len(z));plt.figure(figsize=(8,5));plt.errorbar(z.adjusted_mean_log2CPM_difference,y-.12,xerr=z.SE*student_t.ppf(.975,15),fmt='o',label='Patient + preservation');plt.errorbar(z.adjusted_difference,y+.12,xerr=z.SE_composition*student_t.ppf(.975,13),fmt='s',label='Plus composition scores');plt.yticks(y,z.panel);plt.axvline(0,color='grey',lw=1);plt.xlabel('Adjusted mean log2 CPM difference (PM - primary); 95% t CI');plt.legend(fontsize=8);plt.title('GC: 17 paired patients; exploratory gene panels');fig('GC_03_panel_sensitivity')
m=pd.read_csv(G/'paired_sample_qc.tsv',sep='\t');s=pd.read_csv(G/'diagnostic_panel_sample_scores.tsv',sep='\t');j=m.merge(s,on='sample');plt.figure(figsize=(5,5))
for pat,d in j.groupby('patient'):
 vals=[d.loc[d.site==site,'Lipid_uptake'].iloc[0] for site in ['GCPT','GCPM']];plt.plot([0,1],vals,'o-',alpha=.7,color='#E45756' if vals[1]>vals[0] else '#4C78A8')
plt.xticks([0,1],['Primary','Peritoneal metastasis']);plt.ylabel('Mean log2 CPM, lipid uptake marker panel');plt.title('GC: lipid uptake signal; 13/17 pairs increase');fig('GC_04_paired_lipid_uptake')
p=pd.read_csv(G/'PCA.tsv',sep='\t');plt.figure(figsize=(6,5))
for (site,pre),d in p.groupby(['site','preservation']):plt.scatter(d.PC1,d.PC2,label=f'{site} / {pre}',marker='s' if pre=='FFPE' else 'o')
plt.xlabel('PC1');plt.ylabel('PC2');plt.legend(fontsize=7);plt.title('GC paired RNA-seq: site and preservation');fig('GC_05_PCA')
# Single validation pass for completed outputs, no resampling.
checks=[]
def check(n,c,detail=''):checks.append({'check':n,'pass':bool(c),'detail':detail})
c=pd.read_csv(O/'CRC'/'descriptive_gene_effects.tsv',sep='\t');check('CRC_unique_gene_ID',not c.gene_id.duplicated().any());check('CRC_no_inferential_P',not any(k.lower() in ['p','p.value','fdr'] for k in c.columns));check('GC_17_independent_pairs',m.patient.nunique()==17 and len(m)==34);check('GC_each_patient_one_primary_one_PM',pd.crosstab(m.patient,m.site).eq(1).all().all());check('GC05_recurrence_not_independent_patient',not m['sample'].str.contains('GC05b').any());de=pd.read_csv(G/'paired_preservation_adjusted_DE.tsv',sep='\t');check('GC_gene_P_FDR_valid',de['P.Value'].between(0,1).all() and de['adj.P.Val'].between(0,1).all());check('GC_gene_consistency_totals',((de.positive_pairs+de.negative_pairs)<=17).all());check('GC_panel_FDR_family_size',len(a)==8);check('GC_model_complete',(G/'COMPLETE.txt').exists());check('preanalysis_complete',(O/'COMPLETE.json').exists());check('all_Figures_nonempty',all(p.stat().st_size>1000 for p in F.glob('*.png')))
pd.DataFrame(checks).to_csv(R/'logs'/'validation_checks.tsv',sep='\t',index=False)
# Record all stored inputs without treating .part files as completed data.
manifest=[]
for directory in ['01_metadata','02_raw']:
 for p in sorted((R/directory).glob('*')):
  if not p.is_file():continue
  if p.name.endswith('.part'):manifest.append({'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':'','status':'partial_not_used'});continue
  h=hashlib.sha256()
  with p.open('rb') as f:
   while True:
    chunk=f.read(8*1024*1024)
    if not chunk:break
    h.update(chunk)
  status='stored_nonempty' if p.stat().st_size else 'empty_invalid_not_used'
  manifest.append({'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':h.hexdigest(),'status':status})
pd.DataFrame(manifest).to_csv(R/'01_metadata'/'file_integrity_manifest.tsv',sep='\t',index=False)
print('VALIDATION',checks);print('STORED_VALID_BYTES',sum(x['bytes'] for x in manifest if x['status']=='stored_nonempty'))
