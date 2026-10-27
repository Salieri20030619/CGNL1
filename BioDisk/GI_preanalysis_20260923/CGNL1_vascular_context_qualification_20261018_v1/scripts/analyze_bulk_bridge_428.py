from pathlib import Path
import json,gzip,csv,datetime,hashlib
import openpyxl,numpy as np
b=Path(__file__).resolve().parents[1]
meta={};cur=None
for line in gzip.decompress((b/'raw/GSE151165_family.soft.gz').read_bytes()).decode('utf-8').splitlines():
 if line.startswith('^SAMPLE = '):cur={'gsm':line.split(' = ',1)[1]}
 elif cur is not None and line.startswith('!Sample_title = '):
  title=line.split(' = ',1)[1];sid,group=title.split('_');cur.update(sample=sid,group=group);meta[sid]=cur
 elif cur is not None and line.startswith('!Sample_characteristics_ch1 = '):
  k,v=line.split(' = ',1)[1].split(': ',1);cur[k.lower()]=v
w=openpyxl.load_workbook(b/'raw/GSE151165_RNA_seq.raw_read_count.xlsx',read_only=True,data_only=True)
ws=w['Sheet1'];it=ws.iter_rows(values_only=True);header=next(it);cols=[i for i,v in enumerate(header) if v in meta]
assert len(cols)==30
samples=[header[i] for i in cols];genes=[];values=[]
for row in it:
 if not row[0]:
  assert all(row[i] is None for i in cols);continue
 genes.append(str(row[0]));values.append([row[i] for i in cols])
a=np.array(values,dtype=float);assert np.isfinite(a).all() and (a>=0).all() and (a==np.floor(a)).all()
assert len(set(genes))==len(genes)
markers=['PECAM1','CDH5','CLDN5','ESAM','VWF'];targets=['CGNL1']+markers
assert all(genes.count(g)==1 for g in targets)
libs=a.sum(axis=0);assert (libs>0).all();log=np.log2(a*1e6/libs+1);idx=[genes.index(g) for g in targets];delta=log[idx[0]]-log[idx[1:]].mean(axis=0)
rows=[]
for j,s in enumerate(samples):
 r=dict(meta[s]);r['library_sum']=int(libs[j]);r['CGNL1_relative_to_endothelial']=float(delta[j])
 for g,i in zip(targets,idx):r[g+'_count']=int(a[i,j]);r[g+'_log2CPM1']=float(log[i,j])
 rows.append(r)
outdir=b/'results/bulk_bridge_428';outdir.mkdir(exist_ok=True)
with (outdir/'sample_values.tsv').open('w',encoding='utf-8',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
summary={}
for group in ['DT','RT','DN','RN']:
 inds=[j for j,s in enumerate(samples) if meta[s]['group']==group]
 summary[group]={'sample_count':len(inds),'target_mean_log2CPM1':float(log[idx[0],inds].mean()),'target_relative_mean':float(delta[inds].mean()),'target_counts_range':[int(a[idx[0],inds].min()),int(a[idx[0],inds].max())]}
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'descriptive_complete_inference_withheld_pending_donor_identity','gene_count':len(genes),'raw_values_checked':int(a.size),'groups':summary,'replacement_minus_desmoplastic_tumor_relative_effect':summary['RT']['target_relative_mean']-summary['DT']['target_relative_mean'],'replacement_minus_desmoplastic_tumor_unadjusted_effect':summary['RT']['target_mean_log2CPM1']-summary['DT']['target_mean_log2CPM1'],'primary_p':None,'reason_no_p':'GEO labels unique samples and demographics but no explicit individual identifier; no source-confirmed patient linkage yet. Do not infer pairing from consecutive sample identifiers.','scope':'Bulk expression context only, not endothelial-specific localization, causality or treatment response','source_sha256':hashlib.sha256((b/'raw/GSE151165_RNA_seq.raw_read_count.xlsx').read_bytes()).hexdigest()}
(outdir/'DESCRIPTIVE_COMPLETE_428.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res))
