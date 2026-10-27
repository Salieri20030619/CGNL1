from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import json,time,datetime,hashlib
import pandas as pd,numpy as np
O=Path.cwd()/'mechanism_investigation_v2';receipt=O/'metadata/ROUTE_HTTP_COMPLETE.json'
def live():
 for p in Path('/proc').iterdir():
  if p.name.isdigit():
   try:
    if any(a.endswith(b'/cgnl_mech65_route_http.py') for a in (p/'cmdline').read_bytes().split(b'\0')):return True
   except (FileNotFoundError,PermissionError):pass
 return False
mapping=json.loads((O/'metadata/GRCh37_CGNL1_expanded.json').read_text());assert mapping['assembly_name']=='GRCh37';tx=[t['id'] for t in mapping['Transcript']];pc=[t['id'] for t in mapping['Transcript'] if t['biotype']=='protein_coding'];assert len(tx)==3 and len(pc)==1
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input':'GSE237876作者SalmontranscriptTPM andestimatedNumReads,hg19','target_transcripts':tx,'coding_transcript_sensitivity':pc,'mapping':'officialEnsemblGRCh37symbollookup; exactauthorannotationreleaseunspecified','require':'all3targettranscriptsmustappearuniquelyineachsample; absenttranscriptisnotzero','patient_site_aggregation':'meanlog2(sumtargetTPM+1) acrossavailabletissue sampleswithinpatient/site; coding-onlysensitivity','missingness':'reportexpectedandactualsamplesforeachpatient/site; primaryresultsrestrictedtocompletesites; partialsitesdescriptiveappendixonly','comparisons':'PM/liver/lung/ovary minusprimarywithindonor; smallroutegroupsdescriptive,noinferentialP','sample_count':'66tissues14patients,5expectedPMdonors; no pooledtissue-leveltest'}
(O/'ROUTE_TARGET_SPECIFICATION.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
while not receipt.exists() and live():time.sleep(900)
assert receipt.exists(),'NoHTTPcompletionreceipt; inspectjobexit'
plan=json.loads((O/'metadata/ROUTE_INPUT_PLAN_43.json').read_text());rows=[];coverage=[]
for t in plan:
 title=t['title'];site='primary' if 'primary' in title else ('peritoneal' if 'perit' in title else ('liver' if 'liver' in title else ('lung' if 'lung' in title else ('ovary' if 'ovary' in title else 'unclassified'))));p=O/'raw/GSE237876'/t['files'][0].rsplit('/',1)[1];c={'sample':t['accession'],'patient':t['patient'],'site':site,'title':title,'file_present':p.exists()};coverage.append(c)
 if not p.exists():continue
 a=pd.read_csv(p,sep='\t');assert a.Name.is_unique and np.isfinite(a[['TPM','NumReads']]).all().all() and (a[['TPM','NumReads']]>=0).all().all();present=set(tx)&set(a.Name);c['target_transcripts_covered']=len(present)
 if len(present)!=3:continue
 target=a[a.Name.isin(tx)];coding=a[a.Name.isin(pc)];rows.append(dict(c,TPM=float(target.TPM.sum()),estimated_reads=float(target.NumReads.sum()),coding_TPM=float(coding.TPM.sum()),total_TPM=float(a.TPM.sum()),total_estimated_reads=float(a.NumReads.sum()),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
s=pd.DataFrame(rows);cv=pd.DataFrame(coverage);s['log2TPM1']=np.log2(s.TPM+1);s['coding_log2TPM1']=np.log2(s.coding_TPM+1);agg=s.groupby(['patient','site']).agg(n_samples=('sample','size'),mean_log2TPM1=('log2TPM1','mean'),mean_coding_log2TPM1=('coding_log2TPM1','mean')).reset_index();expected=cv.groupby(['patient','site']).size().rename('expected_samples').reset_index();agg=expected.merge(agg,how='left',on=['patient','site']);agg['complete_site']=agg.n_samples.eq(agg.expected_samples);deltas=[]
for patient,g in agg.groupby('patient'):
 primary=g[g.site=='primary']
 if len(primary)!=1:continue
 p=primary.iloc[0]
 for _,r in g[g.site!='primary'].iterrows():deltas.append({'patient':patient,'site':r.site,'complete_pair':bool(p.complete_site and r.complete_site),'primary_expected':int(p.expected_samples),'primary_available':float(p.n_samples),'metastatic_expected':int(r.expected_samples),'metastatic_available':float(r.n_samples),'delta_log2TPM1':r.mean_log2TPM1-p.mean_log2TPM1,'delta_coding_log2TPM1':r.mean_coding_log2TPM1-p.mean_coding_log2TPM1})
d=pd.DataFrame(deltas)
for n,df in [('sample_values',s),('coverage',cv),('patient_site_values',agg),('paired_deltas',d)]:df.to_csv(O/'results'/('GSE237876_'+n+'.tsv'),sep='\t',index=False,encoding='utf-8')
complete=d[d.complete_pair];summary=complete.groupby('site').agg(n_patients=('patient','nunique'),mean_difference=('delta_log2TPM1','mean'),median_difference=('delta_log2TPM1','median'),positive=('delta_log2TPM1',lambda x:int((x>0).sum())),negative=('delta_log2TPM1',lambda x:int((x<0).sum())));summary.to_csv(O/'results/GSE237876_complete_route_summary.tsv',sep='\t',encoding='utf-8');print(summary.to_string());(O/'results/ROUTE_ANALYSIS_COMPLETE.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'download':json.loads(receipt.read_text()),'qualified_samples':len(s),'expected_samples':66,'complete_pairs':len(complete),'no_inferential_P':True,'mapping_release_limit_retained':True},indent=2),encoding='utf-8')
