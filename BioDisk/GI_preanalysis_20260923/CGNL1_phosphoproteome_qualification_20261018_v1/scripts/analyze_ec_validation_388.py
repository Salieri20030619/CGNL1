from pathlib import Path
import json,csv,math,datetime
import numpy as np
from scipy.stats import rankdata
b=Path(__file__).resolve().parents[1];o=b/'results/EC_phosphosite_validation_388';o.mkdir(exist_ok=True)
a=json.loads((b/'metadata/ESCC_TARGET_EXTRACTION_382.json').read_text());spec=json.loads((b/'metadata/ESCC_ANALYSIS_SPEC_381.json').read_text());allowed=set(spec['endpoints']);prot=a['tables']['b. 9,300 proteins'];phos=a['tables']['c. 61,471 phosphosites'];assert len(prot['target_rows'])==1 and prot['target_rows'][0][0]=='Q0VF96';p=dict(zip(prot['header'][2:],prot['target_rows'][0][2:]));groups={}
for r in phos['target_rows']:
 assert r[1]=='Q0VF96'
 groups.setdefault(tuple(r[4:]),[]).append(r[3])
def exact(d):
 d=np.array(d);d=d[d!=0];r=np.rint(rankdata(abs(d))*2).astype(int);total=int(sum(r));dist=[0]*(total+1);dist[0]=1
 for z in r:
  for j in range(total,z-1,-1):dist[j]+=dist[j-z]
 obs=int(sum(r[d>0]));tail=min(obs,total-obs)
 return min(1.,2*sum(dist[:tail+1])/2**len(d)) if len(d) else None

rows=[];out=[]
for vec,aliases in groups.items():
 matched=sorted(allowed & set(aliases))
 if not matched:continue
 vals=dict(zip(phos['header'][4:],vec));cases=sorted({k.split('_')[0] for k in vals},key=int);ds=[];raw=[];pd=[];included=[]
 for k in cases:
  vt=vals.get(k+'_T');vn=vals.get(k+'_N');pt=p.get(k+'T');pn=p.get(k+'N')
  if not all(isinstance(x,(int,float)) and math.isfinite(x) for x in [vt,vn,pt,pn]):continue
  if vt<=0 or vn<=0 or vt==0.01 or vn==0.01:continue
  sd=math.log2(vt/vn);pr=pt-pn;delta=sd-pr;ds.append(delta);raw.append(sd);pd.append(pr);included.append(k);rows.append({'feature':'/'.join(aliases),'matched_targets':'/'.join(matched),'patient':k,'phospho_T':vt,'phospho_N':vn,'protein_log2_T':pt,'protein_log2_N':pn,'raw_site_delta':sd,'protein_delta':pr,'adjusted_delta':delta})
 out.append({'feature_labels':aliases,'matched_discovery_sites':matched,'n_pairs':len(ds),'mean_adjusted_delta':float(np.mean(ds)) if ds else None,'median_adjusted_delta':float(np.median(ds)) if ds else None,'up':sum(v>0 for v in ds),'down':sum(v<0 for v in ds),'mean_raw_site_delta':float(np.mean(raw)) if raw else None,'mean_total_protein_delta':float(np.mean(pd)) if pd else None,'P':exact(ds) if len(ds)>=10 else None,'FDR':None,'patients':included})
family=sorted([r for r in out if r['P'] is not None],key=lambda r:r['P']);q=1
for j in range(len(family),0,-1):
 r=family[j-1];q=min(q,r['P']*len(family)/j);r['FDR']=q
with (o/'patient_site_values.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_cohort_pairs':31,'matched_measured_features':len(out),'tested_family_size_unique_features':len(family),'results':out,'numerical_QA_pending':True,'cross_cancer_not_PDac_internal_replication':True};(o/'ANALYSIS_COMPLETE_388.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps([{k:v for k,v in r.items() if k!='patients'} for r in out]))
