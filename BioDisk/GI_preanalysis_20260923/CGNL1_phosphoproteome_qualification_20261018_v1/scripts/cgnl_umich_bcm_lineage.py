from pathlib import Path
import json,csv,hashlib,statistics,datetime,collections
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'material_search_restart_20261024';A=D/'lineage';A.mkdir(exist_ok=True)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n');assert json.loads(f.read_text())==z
uf=next((D/'UMich').glob('*.tsv'));mf=P/'metadata/PDC_PDAC_BIOSPECIMENS_394.json'
save(A/'SPEC_BEFORE_COMPARISON.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),purpose='source-lineage QA only,notnewclinicalmodel',inputs={str(f):sha(f) for f in [uf,mf]},mapping='Exact case/type from existingPDCbiospecimen;only unique source column;no numeric matching',sites='all jointly present CGNL1 single-site labels',metrics=['missingmask','BCM-minus-UMich offsets','within-patientS198-S297contrast agreement'],no_new_P=True))
with uf.open() as h:
 r=csv.DictReader(h,delimiter='\t');uh=r.fieldnames;u={}
 for row in r:
  if row['Gene']=='ENSG00000128849.10':u[row['Index'].rsplit('_',1)[1]]=row
maps=collections.defaultdict(set)
for r in json.loads(mf.read_text())['data']['biospecimenPerStudy']:
 for col in [r['aliquot_submitter_id'],r['sample_submitter_id']]:
  if col in uh:maps[(r['case_submitter_id'],r['sample_type'])].add(col)
records=[];mapping=[];skipped=[];bcms={}
for tissue,kind in [('Tumor','Primary Tumor'),('Normal','Solid Tissue Normal')]:
 f=D/'BCM'/('PDAC_phospho_site_abundance_log2_reference_intensity_normalized_'+tissue+'.txt')
 with f.open() as h:
  r=csv.DictReader(h,delimiter='\t');bh=r.fieldnames;b={}
  for row in r:
   if row['idx'].startswith('ENSG00000128849.'):b[row['idx'].split('|')[2]]=row
 bcms[tissue]=b
 for case in bh[1:]:
  cols=maps.get((case,kind),set())
  if len(cols)!=1:skipped.append(dict(case=case,tissue=tissue,source_columns=sorted(cols)));continue
  col=next(iter(cols));mapping.append(dict(case=case,tissue=tissue,source_column=col))
  for site in sorted(set(b)&set(u)):
   bv=b[site][case];uv=u[site][col];bm=bv in ['NA','NaN','nan',''];um=uv in ['NA','NaN','nan','']
   records.append(dict(case=case,tissue=tissue,site=site,source_column=col,bcm_missing=bm,umich_missing=um,bcm=None if bm else float(bv),umich=None if um else float(uv),difference=None if bm or um else float(bv)-float(uv)))
valid=[x['difference'] for x in records if x['difference'] is not None];offset=statistics.median(valid)
summary=[]
for site in sorted(set(x['site'] for x in records)):
 rr=[x for x in records if x['site']==site];dd=[x['difference'] for x in rr if x['difference'] is not None]
 summary.append(dict(site=site,comparisons=len(rr),missing_mismatch=sum(x['bcm_missing']!=x['umich_missing'] for x in rr),joint_numeric=len(dd),median_offset=statistics.median(dd) if dd else None,max_global_offset_residual=max(abs(x-offset) for x in dd) if dd else None))
lookup={(x['case'],x['tissue'],x['site']):x for x in records};contrasts=[]
for case in sorted(set(x['case'] for x in mapping)):
 rr=[lookup.get((case,t,s)) for t,s in [('Tumor','S198'),('Normal','S198'),('Tumor','S297'),('Normal','S297')]]
 if all(x and x['difference'] is not None for x in rr):
  calc=lambda k:rr[0][k]-rr[1][k]-rr[2][k]+rr[3][k]
  contrasts.append(dict(case=case,bcm=calc('bcm'),umich=calc('umich'),residual=calc('bcm')-calc('umich')))
z=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),unique_mapped_samples=len(mapping),skipped=skipped,shared_site_count=len(summary),comparisons=len(records),joint_numeric=len(valid),missing_mask_mismatches=sum(x['bcm_missing']!=x['umich_missing'] for x in records),global_offset=offset,max_global_offset_residual=max(abs(x-offset) for x in valid),per_site=summary,S198_S297_complete_patients=len(contrasts),S198_S297_max_residual=max(abs(x['residual']) for x in contrasts),boundary='Comparison establishes exported numerical relation only;not direct PSM observation or independent biologicalvalidation',new_tests=0)
save(A/'RESULTS.json',z);save(A/'OBSERVATIONS.json',records);save(A/'SAMPLE_MAPPING.json',mapping);save(A/'CONTRASTS.json',contrasts)
print(json.dumps({k:v for k,v in z.items() if k not in ['per_site','skipped']}));print('skipped',json.dumps(skipped)[:1500])
