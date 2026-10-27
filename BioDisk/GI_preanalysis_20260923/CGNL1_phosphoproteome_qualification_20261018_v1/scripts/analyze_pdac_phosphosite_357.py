from pathlib import Path
import csv,json,datetime
import numpy as np
from scipy.stats import rankdata
b=Path(__file__).resolve().parents[1];o=b/'results/PDAC_phosphosite_357';o.mkdir(exist_ok=True)
def load(kind,condition):
 p=next((b/'raw').glob('PDAC_'+kind+'*_'+condition+'.txt'))
 with p.open(encoding='utf-8') as f:
  rd=csv.reader(f,delimiter='\t');h=next(rd);assert len(set(h))==len(h);rows={}
  for r in rd:
   if r[0].split('|')[0].split('.')[0]=='ENSG00000128849':
    assert r[0] not in rows;rows[r[0]]={k:float(v) if v not in ['NA','NaN','nan',''] else np.nan for k,v in zip(h[1:],r[1:])}
 return h[1:],rows
pt,prt=load('proteomics','Tumor');pn,prn=load('proteomics','Normal');st,srt=load('phospho','Tumor');sn,srn=load('phospho','Normal');assert len(prt)==len(prn)==1
cases=sorted(set(pt)&set(pn)&set(st)&set(sn));assert len(cases)>0;protein_id=next(iter(prt));assert protein_id==next(iter(prn));sites=sorted(set(srt)|set(srn));detail=[];out=[]
def exact(d):
 d=np.array(d);d=d[d!=0];r=np.rint(rankdata(abs(d))*2).astype(int);total=int(sum(r));dist=[0]*(total+1);dist[0]=1
 for z in r:
  for j in range(total,z-1,-1):dist[j]+=dist[j-z]
 obs=int(sum(r[d>0]));tail=min(obs,total-obs)
 return min(1.,2*sum(dist[:tail+1])/2**len(d)) if len(d) else None

for site in sites:
 ds=[];raw=[];prot=[]
 for k in cases:
  vals=[srt.get(site,{}).get(k,np.nan),srn.get(site,{}).get(k,np.nan),prt[protein_id][k],prn[protein_id][k]]
  if all(np.isfinite(vals)):
   sd=vals[0]-vals[1];pd=vals[2]-vals[3];d=sd-pd;ds.append(d);raw.append(sd);prot.append(pd);detail.append({'case':k,'site':site,'site_T':vals[0],'site_N':vals[1],'protein_T':vals[2],'protein_N':vals[3],'site_delta':sd,'protein_delta':pd,'adjusted_delta':d})
 out.append({'site':site,'n_pairs':len(ds),'mean_adjusted_delta':float(np.mean(ds)) if ds else None,'median_adjusted_delta':float(np.median(ds)) if ds else None,'up':sum(d>0 for d in ds),'down':sum(d<0 for d in ds),'mean_unadjusted_site_delta':float(np.mean(raw)) if raw else None,'mean_total_protein_delta_same_pairs':float(np.mean(prot)) if prot else None,'P':exact(ds) if len(ds)>=10 else None,'FDR':None})
family=sorted([r for r in out if r['P'] is not None],key=lambda r:r['P']);q=1
for j in range(len(family),0,-1):
 r=family[j-1];q=min(q,r['P']*len(family)/j);r['FDR']=q
with (o/'patient_site_values.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(detail[0]),delimiter='\t');w.writeheader();w.writerows(detail)
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'matched_cases':len(cases),'cases':cases,'total_target_sites':len(sites),'inferential_site_family_size':len(family),'contrasts':out,'QA_pending':True,'source_mask_provenance_pending':True,'independent_validation_pending':True};(o/'ANALYSIS_COMPLETE_357.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in res.items() if k not in ['contrasts','cases']}));print(json.dumps([{k:v for k,v in r.items() if k not in ['mean_total_protein_delta_same_pairs','mean_unadjusted_site_delta']} for r in family]))
