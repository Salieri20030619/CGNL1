from pathlib import Path
import json,gzip,csv,datetime
import numpy as np
from scipy.stats import rankdata
b=Path(__file__).resolve().parents[1];o=b/'results/secondary_pdx_329';o.mkdir(exist_ok=True)
def exact(d):
 d=np.array(d);d=d[d!=0];r=np.rint(rankdata(abs(d))*2).astype(int);total=int(sum(r));dist=[0]*(total+1);dist[0]=1
 for z in r:
  for j in range(total,z-1,-1):dist[j]+=dist[j-z]
 obs=int(sum(r[d>0]));tail=min(obs,total-obs)
 return min(1.,2*sum(dist[:tail+1])/2**len(d)) if len(d) else None

groups=json.loads((b/'metadata/GSE204805_PAIR_QUALIFICATION_316.json').read_text())['groups'];link=json.loads((b/'metadata/DONOR_LINKAGE_324.json').read_text());summaries=[];allrows=[]
for phase in ['PDX_72h','PDX_6w']:
 keys=sorted(link[phase]['case_to_donor']);rs=[r for r in groups[phase] if r['QC_label']=='QC ok' and r['patient_key'] in keys]
 with gzip.open(b/'raw/GSE204805_merged_hs_mm.tsv.gz','rt',encoding='utf-8') as f:
  rd=csv.reader(f,delimiter='\t');head=next(rd);ix=[head.index(r['sample_code']) for r in rs];genes=[];rows=[]
  for row in rd:
   if row[0].startswith('H_'):genes.append(row[0]);rows.append([int(row[i]) for i in ix])
 x=np.array(rows);ref=np.all(x>0,axis=1)&(np.array(genes)!='H_CGNL1');geo=np.exp(np.log(x[ref]).mean(axis=1));sf=np.median(x[ref]/geo[:,None],axis=0);sf/=np.exp(np.log(sf).mean())
 np.savez_compressed(o/(phase+'_input.npz'),counts=x,genes=np.array(genes),samples=np.array([r['sample_code'] for r in rs]),sizefactors=sf)
 (o/(phase+'_samples.json')).write_text(json.dumps(rs,indent=2),encoding='utf-8')
 for method,arr in [('median_ratio',np.log2(x/sf+1)),('human_CPM',np.log2(x/x.sum(axis=0)*1e6+1))]:
  for gene in ['CGNL1','CGN','CDH1','EPCAM','KRT8','VIM','EGFR','DUSP6']:
   gi=genes.index('H_'+gene);ds=[];donors={}
   for k in keys:
    nt=[i for i,r in enumerate(rs) if r['patient_key']==k and r['treatment']=='NT'];ct=[i for i,r in enumerate(rs) if r['patient_key']==k and r['treatment']=='CET'];delta=float(arr[gi,ct].mean()-arr[gi,nt].mean());ds.append(delta);donor=link[phase]['case_to_donor'][k]
    if donor:donors.setdefault(donor,[]).append(delta)
    allrows.append({'phase':phase,'method':method,'gene':gene,'case':k,'donor':donor or 'unmapped','delta':delta})
   dd=np.array([np.mean(v) for v in donors.values()]);summaries.append({'phase':phase,'method':method,'gene':gene,'all_models_n':len(ds),'all_models_mean':float(np.mean(ds)),'all_models_up':sum(z>0 for z in ds),'known_donors_n':len(dd),'known_donors_mean':float(dd.mean()),'known_donors_up':int(sum(dd>0)),'P':exact(dd) if gene=='CGNL1' else None})
main=[r for r in summaries if r['gene']=='CGNL1' and r['method']=='median_ratio'];ordered=sorted(main,key=lambda r:r['P']);q=1
for rank in range(2,0,-1):
 r=ordered[rank-1];q=min(q,r['P']*2/rank);r['FDR_two_secondary_phases']=q
with (o/'model_deltas.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t');w.writeheader();w.writerows(allrows)
(o/'ANALYSIS_COMPLETE_329.json').write_text(json.dumps({'contrasts':summaries,'QA_pending':True},indent=2),encoding='utf-8');print(json.dumps([r for r in summaries if r['gene'] in ['CGNL1','DUSP6']]))
