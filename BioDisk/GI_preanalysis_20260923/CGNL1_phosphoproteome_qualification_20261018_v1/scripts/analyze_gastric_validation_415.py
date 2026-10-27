from pathlib import Path
import json,csv,datetime
import numpy as np
from scipy.stats import rankdata
b=Path(__file__).resolve().parents[1];o=b/'results/GASTRIC_validation_415';o.mkdir(exist_ok=True);p=json.loads((b/'metadata/Gastric_Cancer_Korea_Phosphoproteome.phosphosite.itraq.tsv.CGNL1_412.json').read_text());q=json.loads((b/'metadata/Gastric_Cancer_Korea_Proteome.itraq.tsv.CGNL1_412.json').read_text());prot=dict(zip(q['header'],q['rows'][0]));mapping=json.loads((b/'metadata/GASTRIC_PAIRING_413.json').read_text());assert not mapping['invalid_pairs'];eligible=set(json.loads((b/'metadata/GASTRIC_VALIDATION_SPEC_408.json').read_text())['primary_validation_candidates']);out=[];details=[];tech=[]
def exact(d):
 d=np.array(d);d=d[d!=0];r=np.rint(rankdata(abs(d))*2).astype(int);total=int(sum(r));dist=[0]*(total+1);dist[0]=1
 for z in r:
  for j in range(total,z-1,-1):dist[j]+=dist[j-z]
 obs=int(sum(r[d>0]));tail=min(obs,total-obs)
 return min(1.,2*sum(dist[:tail+1])/2**len(d)) if len(d) else None

for row in p['rows']:
 site=row[0].split(':')[-1].upper()
 if site not in eligible:continue
 v=dict(zip(p['header'],row));patients={}
 for m in mapping['ratio_columns']:
  col=m['column'];pv=prot.get('Unshared '+col,'');sv=v.get(col,'')
  if pv in ['', 'NA','nan','NaN'] or sv in ['', 'NA','nan','NaN']:continue
  sd=float(sv);pd=float(pv);d=sd-pd;patients.setdefault(m['patient'],[]).append((d,sd,pd));tech.append({'site':site,'patient':m['patient'],'channel':col,'site_delta':sd,'protein_delta':pd,'adjusted_delta':d})
 ds=[];ss=[];pp=[]
 for case,rows in patients.items():
  a=np.array(rows).mean(axis=0);ds.append(float(a[0]));ss.append(float(a[1]));pp.append(float(a[2]));details.append({'site':site,'patient':case,'n_technical_channels':len(rows),'adjusted_delta':float(a[0]),'site_delta':float(a[1]),'protein_delta':float(a[2])})
 out.append({'site':site,'n_patients':len(ds),'n_technical_channels':sum(map(len,patients.values())),'mean_adjusted_delta':float(np.mean(ds)) if ds else None,'median_adjusted_delta':float(np.median(ds)) if ds else None,'up':sum(v>0 for v in ds),'down':sum(v<0 for v in ds),'mean_raw_site_delta':float(np.mean(ss)) if ss else None,'mean_protein_delta':float(np.mean(pp)) if pp else None,'P':exact(ds) if len(ds)>=10 else None,'FDR_within_GC':None})
fam=sorted([r for r in out if r['P'] is not None],key=lambda r:r['P']);adj=1
for j in range(len(fam),0,-1):r=fam[j-1];adj=min(adj,r['P']*len(fam)/j);r['FDR_within_GC']=adj
external=[{'cohort':'EC','feature':'/'.join(r['feature_labels']),'P':r['P']} for r in json.loads((b/'results/EC_phosphosite_validation_388/ANALYSIS_COMPLETE_388.json').read_text())['results'] if r['P'] is not None];external.extend({'cohort':'Korean_GC','feature':r['site'],'P':r['P']} for r in fam);external.sort(key=lambda r:r['P']);adj=1
for j in range(len(external),0,-1):r=external[j-1];adj=min(adj,r['P']*len(external)/j);r['FDR_cumulative_external']=adj
for name,rows in [('patient_values.tsv',details),('technical_channel_values.tsv',tech)]:
 with (o/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'results':out,'cumulative_external_family':external,'QA_pending':True,'source_case_population':80};(o/'ANALYSIS_COMPLETE_415.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res))
