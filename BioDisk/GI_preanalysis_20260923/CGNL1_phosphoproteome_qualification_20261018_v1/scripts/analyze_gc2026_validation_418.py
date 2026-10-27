from pathlib import Path
import json,csv,datetime
import numpy as np
from scipy.stats import rankdata
b=Path(__file__).resolve().parents[1];o=b/'results/GC2026_validation_418';o.mkdir(exist_ok=True);ph=json.loads((b/'metadata/GC_TMT_phospho_ratio_GencodeV42_Tumor_NAT.cct.target_416.json').read_text());pr=json.loads((b/'metadata/GC_TMT_proteomics_ratio_gene_MD_GencodeV42_Tumor_NAT.cct.target_416.json').read_text());assert len(pr['rows'])==1;protein=dict(zip(pr['header'],pr['rows'][0]));mp=json.loads((b/'metadata/GC2026_PAIRING_417.json').read_text());assert not mp['identical_target_vectors'];sites=set(json.loads((b/'metadata/GC2026_VALIDATION_SPEC_410.json').read_text())['candidate_sites']);seq=''.join((b/'raw/CGNL1_Q0VF96_380.fasta').read_text().splitlines()[1:]);out=[];detail=[]
def exact(d):
 d=np.array(d);d=d[d!=0];r=np.rint(rankdata(abs(d))*2).astype(int);total=int(sum(r));dist=[0]*(total+1);dist[0]=1
 for z in r:
  for j in range(total,z-1,-1):dist[j]+=dist[j-z]
 obs=int(sum(r[d>0]));tail=min(obs,total-obs)
 return min(1.,2*sum(dist[:tail+1])/2**len(d)) if len(d) else None

for row in ph['rows']:
 parts=row[0].split('|');site=parts[3];motif=parts[4]
 if site not in sites:continue
 pos=int(site[1:]);assert seq[pos-8:pos+7]==motif
 vals=dict(zip(ph['header'],row));ds=[];raw=[];pd=[]
 for case,samples in mp['case_samples'].items():
  if not all(k in samples for k in ['Tumor','NAT']):continue
  t=samples['Tumor'];n=samples['NAT'];v=[vals[t],vals[n],protein[t],protein[n]]
  if any(x in ['NA','nan','NaN',''] for x in v):continue
  a=list(map(float,v));sd=a[0]-a[1];pp=a[2]-a[3];dd=sd-pp;ds.append(dd);raw.append(sd);pd.append(pp);detail.append({'site':site,'case':case,'site_T':a[0],'site_N':a[1],'protein_T':a[2],'protein_N':a[3],'raw_site_delta':sd,'protein_delta':pp,'adjusted_delta':dd})
 out.append({'site':site,'n_patients':len(ds),'mean_adjusted_delta':float(np.mean(ds)) if ds else None,'median_adjusted_delta':float(np.median(ds)) if ds else None,'up':sum(x>0 for x in ds),'down':sum(x<0 for x in ds),'mean_raw_site_delta':float(np.mean(raw)) if raw else None,'mean_protein_delta':float(np.mean(pd)) if pd else None,'P':exact(ds) if len(ds)>=10 else None,'FDR_within_GC2026':None})
fam=sorted([r for r in out if r['P'] is not None],key=lambda r:r['P']);q=1
for j in range(len(fam),0,-1):r=fam[j-1];q=min(q,r['P']*len(fam)/j);r['FDR_within_GC2026']=q
external=json.loads((b/'results/GASTRIC_validation_415/ANALYSIS_COMPLETE_415.json').read_text())['cumulative_external_family'];external=[{k:r[k] for k in ['cohort','feature','P']} for r in external];external.extend({'cohort':'CPTAC_GC2026','feature':r['site'],'P':r['P']} for r in fam);external.sort(key=lambda r:r['P']);q=1
for j in range(len(external),0,-1):r=external[j-1];q=min(q,r['P']*len(external)/j);r['FDR_cumulative_external']=q
with (o/'patient_values.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(detail[0]),delimiter='\t');w.writeheader();w.writerows(detail)
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'results':out,'cumulative_external_family':external,'source_pairs':30,'QA_pending':True,'PPI_target_qualification':'No CGNL1-identified quantitative crosslink row found in supplied linkage matrix; not evidence of absence of interactions'};(o/'ANALYSIS_COMPLETE_418.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res))
