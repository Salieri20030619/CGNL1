from pathlib import Path
import json,csv,math,hashlib,statistics,datetime,re
O=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');W=O/'C04_search_pilot';V=W/'run_v2_all_spectra'
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n');assert json.loads(f.read_text())==z
spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),timing='engine running,only outputheader seen,no finishedsearchresults inspected',unit='topPSMper scan-chargequery;scan-duplicate counts reported',FDR='tiesgroupedbyEvalue;min1((D+1)/T);reversecumulativeminimum;alltargetdecoyrows',phospho_scope='add separate phosphoPSM stratum using modifications mass79.966;claimphosphoidentification onlyif globalq<=.01 andphosphoq<=.01',ambiguous_target_decoy='anyDECOY protein => decoy',contaminants='includedinsearch/FDR but flagged notCGNL1support',site='AScoreProreportedwithoutprobabilityconversion orsiteFDRclaim',no_biological_effect_test=True)
f=W/'POSTPROCESS_SPEC_BEFORE_RESULTS.json'
if not f.exists():save(f,spec)
if __import__('sys').argv[-1]=='freeze':print('Postprocess frozen');raise SystemExit
assert json.loads((V/'SEARCH_EXIT.json').read_text())['exit']==0
csv.field_size_limit(20000000)
with (V/'pilot_hcd.txt').open() as h:
 version=h.readline();r=csv.DictReader(h,delimiter='\t');rows=list(r)
assert all(int(x['num'])==1 for x in rows)
for i,x in enumerate(rows):
 x['_id']=i;x['_e']=float(x['e-value']);assert math.isfinite(x['_e']) and x['_e']>=0
 x['_decoy']='DECOY_' in x['protein'];x['_contam']='CON__' in x['protein'];x['_phospho']='79.966' in x['modifications']
def qvals(rr,key):
 groups={}
 for x in rr:groups.setdefault(x['_e'],[]).append(x)
 target=decoy=0;ranked=[]
 for e,g in sorted(groups.items()):
  decoy+=sum(x['_decoy'] for x in g);target+=sum(not x['_decoy'] for x in g);ranked.append((g,min(1,(decoy+1)/max(target,1))))
 minimum=1.
 for g,fdr in reversed(ranked):
  minimum=min(minimum,fdr)
  for x in g:x[key]=minimum
qvals(rows,'q_global');phos=[x for x in rows if x['_phospho']];qvals(phos,'q_phospho')
accepted=[x for x in rows if not x['_decoy'] and x['q_global']<=.01]
acceptedphos=[x for x in phos if not x['_decoy'] and x['q_global']<=.01 and x['q_phospho']<=.01]
accessions=[x.split('|')[1] for x in (W/'mouse_UP000000589.fasta').read_text().splitlines() if x.startswith('>') and ' GN=Cgnl1 ' in x]
target=[x for x in rows if not x['_decoy'] and any(a in x['protein'] for a in accessions)]
errors=[]
for x in acceptedphos:
 delta=float(x['exp_neutral_mass'])-float(x['calc_neutral_mass']);calc=float(x['calc_neutral_mass']);errors.append(min(abs(delta-k*1.0033548378) for k in [0,1,2])/calc*1e6)
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),total_top_PSMs=len(rows),unique_scan_charge=len(set((x['scan'],x['charge']) for x in rows)),unique_scans=len(set(x['scan'] for x in rows)),target_PSMs=sum(not x['_decoy'] for x in rows),decoy_PSMs=sum(x['_decoy'] for x in rows),phospho_PSMs=len(phos),accepted_global_target_PSMs=len(accepted),accepted_phospho_target_PSMs=len(acceptedphos),accepted_phospho_contaminant_PSMs=sum(x['_contam'] for x in acceptedphos),accepted_phospho_median_abs_corrected_ppm=statistics.median(errors) if errors else None,CGNL1_candidate_PSMs=len(target),CGNL1_accepted_phospho_PSMs=sum(x in acceptedphos for x in target),CGNL1_accessions=accessions,site_localization_qualified=False,protein_FDR_not_established=True,biological_comparison_not_run=True)
save(V/'IDENTIFICATION_SUMMARY.json',result);save(V/'CGNL1_CANDIDATE_PSMS.json',target);save(V/'ALL_PSM_QVALUES.json',rows)
print(json.dumps(result));print('CGNL1candidates',json.dumps(target)[:6000])
