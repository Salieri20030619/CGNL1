from pathlib import Path
import json,csv,hashlib,datetime,math,collections
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'H02_localization_reaudit_20261024'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n');assert read(f)==z
source=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576/H02/selected_585/Phospho (STY)Sites.txt');assert sha(source)=='1b6987f36b7b8f64ef4a3c5ffc5d08c8c1a8f83e1e19c793011741d05d2c8b67'
save(D/'SPEC_BEFORE_COUNTS.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=sha(source),scope='S283/S297 in the17patients of frozenC03_H02 contrast',rule='Count experiment-levelMaxQuant Localizationprob>=0.75 for each of4measurements;missing andbelowthreshold separate;do not infer directPSM for MBR',no_new_effect_test=True,do_not_select_patients_by_result=True))
with (O/'wave1/patient_contrasts.csv').open() as h:primary=[x for x in csv.DictReader(h) if x['test']=='C03_H02']
print('primary columns',list(primary[0]));cases={x['case'] for x in primary};assert len(cases)==17
rows={x['Amino acid']+x['Position']:x for x in read(D/'SOURCE_TARGET_ROWS.json')};pairs={x['case']:x for x in read(P/'results/H02_recovered_590/PAIRING.json')['eligible']}
out=[]
def num(s):
 if s in ['',None,'NA','NaN']:return None
 x=float(s);return x if math.isfinite(x) else None
for case in sorted(cases):
 for typ in ['T','P']:
  exp=pairs[case]['raw']+'_'+typ
  for site in ['S283','S297']:
   row=rows[site];key='Localization prob '+exp;assert key in row
   prob=num(row[key]);intensity=num(row['Intensity '+exp]);parts=[num(row.get('Intensity '+exp+'___'+str(i))) for i in [1,2,3]]
   out.append(dict(case=case,experiment=exp,site=site,localization_prob=prob,prob_status='missing' if prob is None else 'pass_0.75' if prob>=.75 else 'below_0.75',intensity=intensity,multiplicity_intensities=parts))
counts={site:dict(collections.Counter(x['prob_status'] for x in out if x['site']==site)) for site in ['S283','S297']}
passing=[c for c in sorted(cases) if all(x['prob_status']=='pass_0.75' for x in out if x['case']==c)]
mult={site:dict(positive_observations=sum((x['intensity'] or 0)>0 for x in out if x['site']==site),any_multiply_phosphorylated_signal=sum(any((v or 0)>0 for v in x['multiplicity_intensities'][1:]) for x in out if x['site']==site),mono_positive=sum((x['multiplicity_intensities'][0] or 0)>0 for x in out if x['site']==site)) for site in ['S283','S297']}
z=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256=sha(source),primary_patients=17,observations=len(out),by_site=counts,all_four_localization_pass_patients=passing,all_four_localization_pass_n=len(passing),multiplicity_counts=mult,new_tests=0,correction='H02 original file does contain experiment-level localization probability columns;previousblanketabsence wording incorrect',boundary='OriginalH02analysis usedglobalprob>=0.75 andpositiveintensity;experimentprobQA doesnot retroactivelyreplace17patientprimarynull or establish allPSMs localized')
save(D/'OBSERVATION_AUDIT.json',out);save(D/'RESULTS.json',z);print(json.dumps(z))
