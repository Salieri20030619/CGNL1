from pathlib import Path
import json,gzip,csv,datetime
import numpy as np
from scipy.stats import rankdata
b=Path(__file__).resolve().parents[1];o=b/'results/independent_pdo_326';o.mkdir(exist_ok=True)
rs=json.loads((b/'metadata/GSE204805_PAIR_QUALIFICATION_316.json').read_text())['groups']['PDO_72h']
rs=[r for r in rs if r['QC_label']=='QC ok']
keys=sorted({r['patient_key'] for r in rs});keys=[k for k in keys if {r['treatment'] for r in rs if r['patient_key']==k}=={'CET','NT'}]
rs=[r for r in rs if r['patient_key'] in keys]; assert len(keys)==16 and len(rs)==62
with gzip.open(b/'raw/GSE204805_merged_hs_mm.tsv.gz','rt',encoding='utf-8') as f:
 rd=csv.reader(f,delimiter='\t');head=next(rd);assert len(set(head))==len(head)
 ix=[head.index(r['sample_code']) for r in rs]; genes=[];rows=[];mouse=np.zeros(len(rs),dtype=np.int64)
 for row in rd:
  assert len(row)==len(head)
  vals=np.array([int(row[i]) for i in ix],dtype=np.int64);assert np.all(vals>=0)
  if row[0].startswith('H_'):genes.append(row[0]);rows.append(vals)
  elif row[0].startswith('M_'):mouse+=vals
  else:raise ValueError(row[0])
x=np.array(rows);assert len(set(genes))==len(genes)
ref=np.all(x>0,axis=1)&(np.array(genes)!='H_CGNL1');assert sum(ref)>1000
geo=np.exp(np.log(x[ref]).mean(axis=1));sf=np.median(x[ref]/geo[:,None],axis=0);sf/=np.exp(np.log(sf).mean())
y=np.log2(x/sf+1);cpm=np.log2(x/x.sum(axis=0)*1e6+1)
np.savez_compressed(o/'input_and_normalization.npz',counts=x,genes=np.array(genes),samples=np.array([r['sample_code'] for r in rs]),sizefactors=sf,mouse=mouse)
(o/'samples.json').write_text(json.dumps(rs,indent=2),encoding='utf-8')
def exact(d):
 d=np.array(d);d=d[d!=0];r=np.rint(rankdata(abs(d))*2).astype(int);total=int(sum(r));dist=[0]*(total+1);dist[0]=1
 for z in r:
  for j in range(total,z-1,-1):dist[j]+=dist[j-z]
 obs=int(sum(r[d>0]));tail=min(obs,total-obs)
 return min(1.,2*sum(dist[:tail+1])/2**len(d)) if len(d) else None
out=[];delta_rows=[]
for method,arr in [('median_ratio',y),('human_CPM',cpm)]:
 for gene in ['CGNL1','CGN','CDH1','EPCAM','KRT8','VIM','EGFR','DUSP6']:
  gi=genes.index('H_'+gene);ds=[]
  for k in keys:
   nt=[i for i,r in enumerate(rs) if r['patient_key']==k and r['treatment']=='NT'];ct=[i for i,r in enumerate(rs) if r['patient_key']==k and r['treatment']=='CET']
   a=float(arr[gi,nt].mean());z=float(arr[gi,ct].mean());ds.append(z-a)
   delta_rows.append({'method':method,'gene':gene,'case':k,'NT':a,'CET':z,'delta':z-a,'NT_replicates':len(nt),'CET_replicates':len(ct)})
  out.append({'method':method,'gene':gene,'n':len(ds),'mean_delta':float(np.mean(ds)),'median_delta':float(np.median(ds)),'up':sum(d>0 for d in ds),'down':sum(d<0 for d in ds),'P':exact(ds) if gene=='CGNL1' and x[gi].sum()>0 else None,'role':'primary' if method=='median_ratio' and gene=='CGNL1' else 'sensitivity' if gene=='CGNL1' else 'descriptive'})
with (o/'patient_deltas.tsv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(delta_rows[0]),delimiter='\t');w.writeheader();w.writerows(delta_rows)
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'samples':len(rs),'independent_donors':16,'human_genes':len(genes),'reference_genes':int(sum(ref)),'human_fraction_range':[float(v) for v in (min(x.sum(axis=0)/(x.sum(axis=0)+mouse)),max(x.sum(axis=0)/(x.sum(axis=0)+mouse)))],'target_count_range':[int(x[genes.index('H_CGNL1')].min()),int(x[genes.index('H_CGNL1')].max())],'contrasts':out,'independent_QA_pending':True}
(o/'ANALYSIS_COMPLETE_326.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res))
