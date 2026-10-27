from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,time,datetime,hashlib,itertools
import numpy as np,pandas as pd
from scipy import stats
O=Path.cwd()/'mechanism_investigation_v2';done=O/'metadata/DRA_QUANTIFICATION_COMPLETE.json'
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':'prospectivestatisticalimplementationbeforeanyDRAquantificationoutcomes','intent_reference':'DRA_ANALYSIS_INTENT.json','target_gene':'ENSG00000128849','measure':'sumtranscriptTPMpergene,log2(TPM+1);estimatedNumReadsreportedseparately','primary':'combinedcytokinesminuscontrolwithinsameGCdonor','primary_test':'exactpairedsignedrank,3donors; oneprespecifiedprimaryendpointnoFDR','secondary':'TGFminuscontrolandIL1minuscontrol;BHfamilyoftwo','sensitivity':'pairedtP/95percentmeanCI,assumesapproximatelynormalpairedlogdifferences,assumptioncannotbeestablishedatn3','interpretation':'allresultsremainexploratory;minimumattainableexacttwo-sidedPfor3nonzeroequal-directionpairsis0.25,so lackofnominalsignificancenotproof ofnodirectionalresponse','replicate_definition':'patientsGC40GC41GC42,not12librariesornumberreads','no_test_when':'allzeroresponse,less3qualifiedpairs,missingmeasurement;reportstatusanddirectiononly','additional_fixed_genes':['CGN','UPK3B','LRRN4','WT1','CALB2','MSLN','CAMSAP3','MYH9','MYH10','TNC','FN1','COL1A1','ACTA2','TAGLN'],'context_genes_role':'descriptiveonly,noscreeningorselectingnominalP','normalization_sensitivity':'log2(estimatedgeneNumReads/sumgeneNumReads*1e6+1)','limits':'cytokineperturbationisnotCGNL1perturbation; noread-levelbiologicalP; nofoldclaimfromexponentiatingpseudocountdifference'}
sp=O/'DRA_TARGET_STATISTICAL_SPECIFICATION.json';assert sp.exists();print('DRA statistical specification frozen',flush=True)
def live():
 for p in Path('/proc').iterdir():
  if p.name.isdigit():
   try:
    if any(a.endswith(b'/cgnl_mech46b_quantify_when_ready.py') for a in (p/'cmdline').read_bytes().split(b'\0')):return True
   except (FileNotFoundError,PermissionError):pass
 return False
while not done.exists() and live():time.sleep(900)
assert done.exists(),'Quantifier exited without completion receipt';receipt=json.loads(done.read_text());assert receipt['all_pass'],'Incomplete quantification; review missing/failing runs before statistics'
plan=json.loads((O/'metadata/DRA_DOWNLOAD_PLAN.json').read_text());runs={t['run']:t for t in plan['tasks']};rows=[];context=[]
for run,t in runs.items():
 p=O/'processed/DRA014736_quant'/run;g=pd.read_csv(p/'gene_quant.tsv.gz',sep='\t');meta=json.loads((p/'aux_info/meta_info.json').read_text());geneid=g.gene_id.str.split('.').str[0];v=g[geneid=='ENSG00000128849'];assert len(v)==1 and v.gene_symbol.iloc[0]=='CGNL1';r=v.iloc[0];total=float(g.NumReads.sum());assert total>0
 rows.append({'run':run,'donor':t['donor'],'condition':t['condition'],'TPM':float(r.TPM),'estimated_reads':float(r.NumReads),'total_estimated_gene_reads':total,'log2TPM1':float(np.log2(r.TPM+1)),'log2CPM1':float(np.log2(r.NumReads/total*1e6+1)),'mapping_percent':meta.get('percent_mapped'),'gene_quant_sha256':hashlib.sha256((p/'gene_quant.tsv.gz').read_bytes()).hexdigest()})
 for sym in spec['additional_fixed_genes']:
  vv=g[g.gene_symbol==sym];context.append({'run':run,'donor':t['donor'],'condition':t['condition'],'gene':sym,'gene_rows':len(vv),'TPM':float(vv.TPM.sum()) if len(vv)==1 else np.nan})
x=pd.DataFrame(rows);assert len(x)==12 and x.donor.nunique()==3 and not x[['donor','condition']].duplicated().any();out=[];deltas=[]
for measure in ['log2TPM1','log2CPM1']:
 v=x.pivot(index='donor',columns='condition',values=measure)
 for cond in ['combined','TGF','IL1']:
  d=(v[cond]-v.control).dropna();r={'measure':measure,'contrast':cond+'_minus_control','n_donors':len(d),'mean_difference':float(d.mean()),'median_difference':float(d.median()),'positive':int((d>0).sum()),'negative':int((d<0).sum()),'P':np.nan,'FDR_secondary':np.nan,'paired_t_P':np.nan,'mean_CI_low':np.nan,'mean_CI_high':np.nan,'status':'descriptive_only'}
  if len(d)==3 and (d!=0).any():
   nonzero=d.to_numpy()[d.to_numpy()!=0];ranks=stats.rankdata(abs(nonzero));observed=abs(np.dot(np.sign(nonzero),ranks));exact=float(np.mean([abs(np.dot(sign,ranks))>=observed-1e-12 for sign in itertools.product([-1,1],repeat=len(nonzero))]));r.update(P=exact,status='exploratory_three_donors',nonzero_pairs=len(nonzero))
   se=stats.sem(d)
   if se>0:
    ci=stats.t.ppf(.975,2)*se;r.update(paired_t_P=float(stats.ttest_1samp(d,0).pvalue),mean_CI_low=float(d.mean()-ci),mean_CI_high=float(d.mean()+ci))
  out.append(r);deltas.extend({'donor':donor,'measure':measure,'contrast':cond+'_minus_control','difference':float(val)} for donor,val in d.items())
r=pd.DataFrame(out)
for measure in ['log2TPM1','log2CPM1']:
 ix=r.index[(r.measure==measure)&r.contrast.ne('combined_minus_control')];pvals=r.loc[ix,'P'].fillna(1);q=stats.false_discovery_control(pvals.to_numpy());r.loc[ix,'FDR_secondary']=q;r.loc[ix[r.loc[ix,'P'].isna()],'FDR_secondary']=np.nan
for name,df in [('target_samples',x),('target_donor_deltas',pd.DataFrame(deltas)),('target_contrasts',r),('context_gene_samples',pd.DataFrame(context))]:df.to_csv(O/'results'/('DRA014736_'+name+'.tsv'),sep='\t',index=False,encoding='utf-8')
(O/'results/DRA_TARGET_ANALYSIS_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'spec_sha256':hashlib.sha256(sp.read_bytes()).hexdigest(),'independent_donors':3,'libraries':12,'all_inputs_quantified':True,'mechanism_established':False,'independent_numeric_QA_pending':True},indent=2),encoding='utf-8');print(r.to_string(index=False),flush=True)
