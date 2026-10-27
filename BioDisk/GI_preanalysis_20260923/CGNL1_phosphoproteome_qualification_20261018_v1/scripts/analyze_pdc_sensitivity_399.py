from pathlib import Path
import json,csv,datetime
import numpy as np
from scipy.stats import rankdata
b=Path(__file__).resolve().parents[1];o=b/'results/PDC_processing_sensitivity_399';o.mkdir(exist_ok=True)
a=json.loads((b/'results/PDAC_phosphosite_357/ANALYSIS_COMPLETE_357.json').read_text());allowed={r['site'].split('|')[2] for r in a['contrasts'] if r['P'] is not None};mapping=json.loads((b/'metadata/PDC_ACTIVE_SAMPLE_QUALIFICATION_395.json').read_text())['active_map'];ph=json.loads((b/'metadata/PDC_PHOSPHO_TARGET_393.json').read_text());pr=json.loads((b/'metadata/PDC_PROTEIN_TARGET_397.json').read_text());assert len(pr['rows'])==1;protein=dict(zip(pr['header'],pr['rows'][0]));out=[];detail=[]
def exact(d):
 d=np.array(d);d=d[d!=0];r=np.rint(rankdata(abs(d))*2).astype(int);total=int(sum(r));dist=[0]*(total+1);dist[0]=1
 for z in r:
  for j in range(total,z-1,-1):dist[j]+=dist[j-z]
 obs=int(sum(r[d>0]));tail=min(obs,total-obs)
 return min(1.,2*sum(dist[:tail+1])/2**len(d)) if len(d) else None

def values(row,aliquots,suffix):
 return [float(row[k+suffix]) for k in aliquots if row.get(k+suffix,'') not in ['', 'NA','NaN','nan']]
for method,suffix in [('unshared',' Unshared Log Ratio'),('all_peptides',' Log Ratio')]:
 for sr in ph['rows']:
  site=sr[0].split(':')[-1].upper()
  if site not in allowed:continue
  row=dict(zip(ph['header'],sr));deltas=[];raw=[];pd=[]
  for case in a['cases']:
   c=mapping[case];ti=c.get('Primary Tumor',[]);ni=c.get('Solid Tissue Normal',[]);st=values(row,ti,' Log Ratio');sn=values(row,ni,' Log Ratio');pt=values(protein,ti,suffix);pn=values(protein,ni,suffix)
   if not all([st,sn,pt,pn]):continue
   sd=float(np.mean(st)-np.mean(sn));pp=float(np.mean(pt)-np.mean(pn));dd=sd-pp;deltas.append(dd);raw.append(sd);pd.append(pp);detail.append({'method':method,'site':site,'case':case,'site_delta':sd,'protein_delta':pp,'adjusted_delta':dd,'n_site_T':len(st),'n_site_N':len(sn),'n_protein_T':len(pt),'n_protein_N':len(pn)})
  out.append({'method':method,'site':site,'source_id':sr[0],'n':len(deltas),'mean_adjusted':float(np.mean(deltas)) if deltas else None,'up':sum(v>0 for v in deltas),'down':sum(v<0 for v in deltas),'P':exact(deltas) if len(deltas)>=10 else None,'FDR':None,'mean_raw_site':float(np.mean(raw)) if raw else None,'mean_protein':float(np.mean(pd)) if pd else None})
for method in ['unshared','all_peptides']:
 fam=sorted([r for r in out if r['method']==method and r['P'] is not None],key=lambda r:r['P']);q=1
 for j in range(len(fam),0,-1):r=fam[j-1];q=min(q,r['P']*len(fam)/j);r['FDR']=q
with (o/'patient_values.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(detail[0]),delimiter='\t');w.writeheader();w.writerows(detail)
res={'results':out,'QA_pending':True,'not_independent_biological_validation':True};(o/'ANALYSIS_COMPLETE_399.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(out))
