from common import *
import subprocess,zipfile
D=Q/'V05_GATE';D.mkdir(exist_ok=True);m=Q/'V03_MIXTURE/summary.tsv';e=Q/'V04_SPATIAL/effects.tsv';ok=False;reasons=[]
if m.exists() and e.exists() and e.stat().st_size>5:
 x=pd.read_csv(m,sep='\t');x=x[(x.gene=='CGNL1')&(x['mask']=='strict5')];y=pd.read_csv(e,sep='\t');y=y[y.gene=='CGNL1']
 if len(x):
  r=x.iloc[0];ok=r.n_patients>=5 and r.mean_effect>0 and r.positive>=np.ceil(.7*r.n_patients) and len(y)>=2 and (y.effect>0).all()
  if r.n_patients>=10:ok=ok and pd.notna(r.FDR) and r.FDR<.05
else:reasons.append('No adequate strict residual or external effect table')
if not ok:dump(D/'RESULT.json',{'status':'skipped_evidence_gate','reason':'Insufficient consistent source-aware spatial evidence for automatic mechanism expansion','detail':reasons,'interpretation':'This is a prioritization gate, not proof of no biological function.'});sys.exit(0)
j=json.loads((S/'metadata/mosaic_liver_zenodo.json').read_text());f=next(f for f in j['files'] if f.get('key',f.get('name'))=='BulkRNAseq.zip');url=f['links'].get('content',f['links'].get('self'));p=D/'BulkRNAseq.zip'
with (D/'download.log').open('w') as log:r=subprocess.run(['curl','-L','--fail','--retry','2','--retry-all-errors','--max-time','1800','-C','-','-o',str(p)+'.part',url],stdout=log,stderr=subprocess.STDOUT)
assert r.returncode==0,'perturbation data download failed';pathlib.Path(str(p)+'.part').replace(p)
with zipfile.ZipFile(p) as z:
 assert z.testzip() is None;pd.DataFrame([{'file':f.filename,'bytes':f.file_size} for f in z.infolist()]).to_csv(D/'perturbation_inventory.tsv',sep='\t',index=False)
dump(D/'RESULT.json',{'status':'perturbation_assets_ready_pending_design_review','scope':'Treatment/replicate design must be checked before any response comparison; not a validated CGNL1 mechanism.'})
