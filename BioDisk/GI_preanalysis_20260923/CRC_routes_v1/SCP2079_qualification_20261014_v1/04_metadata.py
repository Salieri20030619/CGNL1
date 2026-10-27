from pathlib import Path
import requests,json,hashlib,concurrent.futures,pandas as pd,io,datetime
ROOT=Path(__file__).resolve().parent
D=ROOT/'metadata';D.mkdir(exist_ok=True)
annotations=json.loads((ROOT/'sources'/'annotations.txt').read_text(encoding='utf-8'))
selected=[a for a in annotations if a['scope']=='study']
base='https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP2079/annotations/'
def fetch(a):
 name=a['name'];p=D/(name+'.tsv')
 url=base+name+'/cell_values'
 params={'annotation_type':a['type'],'annotation_scope':'study','cluster':'Epithelial Cells'}
 r=requests.get(url,params=params,timeout=90);r.raise_for_status()
 p.write_bytes(r.content)
 df=pd.read_csv(io.BytesIO(r.content),sep='\t',dtype=str,keep_default_na=False)
 assert df.shape[1]==2 and df['NAME'].is_unique
 print(name,len(df),flush=True)
 return dict(name=name,url=r.url,status=r.status_code,rows=len(df),sha256=hashlib.sha256(r.content).hexdigest()),df.set_index('NAME')
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: results=list(ex.map(fetch,selected))
idx=results[0][1].index
assert all(set(df.index)==set(idx) for _,df in results)
meta=pd.concat([df.reindex(idx) for _,df in results],axis=1)
assert len(meta)==95421 and meta.index.is_unique
meta.to_csv(ROOT/'metadata_all_cells.tsv',sep='\t',encoding='utf-8')
(ROOT/'METADATA_FETCH.json').write_text(json.dumps([d for d,_ in results],indent=2),encoding='utf-8')
print('COLUMNS',list(meta.columns),flush=True)
for key in ['donor_id','patientID','biosample_id','specimenID','sampleID','treatment','tissueSite']:
 print(key,meta[key].nunique(),flush=True)
assert (meta.donor_id==meta.patientID).all()
assert (meta.biosample_id==meta.specimenID).all()
assert (meta.specimenID.str[:3]==meta.donor_id).all()
assert ((meta.treatment=='Pre')==(meta.specimenID.str.endswith('D01'))).all()
summary=meta.groupby(['donor_id','treatment','tissueSite']).size().rename('n_cells').reset_index()
summary.to_csv(ROOT/'patient_time_site_coverage.tsv',sep='\t',index=False,encoding='utf-8')
print(summary.to_string(index=False),flush=True)
for col in ['total_counts','nCount_RNA','nFeature_RNA','pct_counts_mt','doublet_score']:
 v=pd.to_numeric(meta[col],errors='raise');print(col,v.describe().to_dict(),flush=True)