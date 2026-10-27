from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
import numpy as np,pandas as pd,json,gzip,hashlib,datetime,csv
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1');R=B/'new_cohorts_20261013_v1';O=R/'results'
S=json.loads((R/'SPECIFICATION.json').read_text())
amend=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),metadata_basis='GSE200997 annotation contains sample,Condition,Location,MSI_Status,CMS predictions but no cell lineage. GSE132465 contains Patient,Class,Sample,Cell_type,Cell_subtype.',
 target_associations_seen=False,
 change='GSE200997 uses frozen target-excluded lineage rules as primary. GSE132465 primary remains author epithelial+QC, and frozen marker intersection sensitivity. Missing author labels are not inferred from CMS predictions.',
 strict_subset='For GSE200997, marker margin>=0.8 (vs frozen0.4 main) as pre-result identity sensitivity; both use >=2 markers. Not a new malignant-cell call.')
p=R/'ANNOTATION_AMENDMENT.json'
if not p.exists():p.write_text(json.dumps(amend,indent=2));(R/'ANNOTATION_AMENDMENT_FREEZE.json').write_text(json.dumps(dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest())))
needed=sorted(set(['CGNL1']+S['ordinary_markers']+sum(S['panels'].values(),[])+sum(S['outcomes'].values(),[])+sum(S['lineage_markers'].values(),[])))
profiles=[];coverage=[];qa={};sourcegenes={}
manifest=json.loads((R/'DOWNLOAD_MANIFEST.json').read_text())
for ds,ext,sep in [('GSE200997','csv',','),('GSE132465','txt','\t')]:
 raw=R/'raw'/(ds+'_GEO_processed_CRC_10X_raw_UMI_count_matrix.'+ext+'.gz')
 ann=R/'metadata'/(ds+'_GEO_processed_CRC_10X_cell_annotation.'+ext+'.gz')
 ok=all(any(x['path']==str(p) and x['status']=='verified' for x in manifest) for p in [raw,ann])
 if not ok:qa[ds]={'status':'download_not_verified'};continue
 m=pd.read_csv(ann,sep=sep,index_col=0);assert m.index.is_unique
 print('STREAM_START',ds,flush=True)
 with gzip.open(raw,'rt') as f:
  header=next(csv.reader([f.readline()],delimiter=sep))
  cells=np.array(header[1:]);assert len(cells)>1000 and pd.Index(cells).is_unique
  assert set(cells)==set(m.index),(ds,len(cells),len(m))
  m=m.reindex(cells)
  n=len(cells);umi=np.zeros(n);ng=np.zeros(n,dtype=np.int64);mt=np.zeros(n);small={};names=[]
  for line in f:
   gene,values=line.rstrip('\n\r').split(sep,1);gene=gene.strip('"')
   v=np.fromstring(values,sep=sep,dtype=np.float64)
   assert len(v)==n,(ds,gene,len(v),n)
   assert np.isfinite(v).all() and (v>=0).all() and np.allclose(v,np.round(v),atol=1e-8,rtol=0)
   umi+=v;ng+=(v>0);names.append(gene)
   if gene.startswith('MT-'):mt+=v
   if gene in needed:small[gene]=small.get(gene,np.zeros(n))+v
  mt=mt/np.maximum(umi,1)*100
 sourcegenes[ds]=names
 pd.Series(names,name='gene').to_csv(R/'metadata'/(ds+'_genes.tsv'),sep='\t',index=False)
 assert all(g in small for g in ['CGNL1']+S['outcomes']['actomyosin'])
 good=(umi>=500)&(ng>=200)&(mt<=25)
 labels=list(S['lineage_markers']);scores=[];hits=[]
 for label,gs in S['lineage_markers'].items():
  X=np.array([small[g] for g in gs if g in small]);assert len(X)>=.8*len(gs)
  scores.append(np.log1p(X/np.maximum(umi,1)*1e4).mean(0));hits.append((X>0).sum(0))
 scores=np.array(scores).T;hits=np.array(hits).T;j=labels.index('epithelial');margin=np.sort(scores,axis=1)[:,-1]-np.sort(scores,axis=1)[:,-2]
 marker=good&(scores.argmax(1)==j)&(hits[:,j]>=2)&(margin>=.4)
 if ds=='GSE200997':
  tumor=m.Condition.eq('Tumor').to_numpy();patients=m.samples.str.replace(r'^[BT]_','',regex=True).to_numpy();samples=m.samples.to_numpy()
  defs={'main':marker,'identity_sensitivity':marker&(margin>=.8)};author=np.full(n,False)
 else:
  tumor=m.Class.eq('Tumor').to_numpy();patients=m.Patient.astype(str).to_numpy();samples=m.Sample.astype(str).to_numpy()
  author=m.Cell_type.eq('Epithelial cells').to_numpy()
  defs={'main':good&author,'identity_sensitivity':good&author&marker}
 pd.DataFrame(dict(cell=cells,patient=patients,sample=samples,tumor=tumor,UMI=umi,genes=ng,mito_pct=mt,QC=good,author_epithelial=author,frozen_epithelial=marker,marker_margin=margin)).to_csv(O/(ds+'_cell_QC.tsv.gz'),sep='\t',index=False)
 # Cache target and panels only, not the complete dense expression matrix.
 np.savez_compressed(O/(ds+'_panel_counts.npz'),genes=np.array(sorted(small)),counts=np.array([small[g] for g in sorted(small)]),cells=cells,UMI=umi)
 for definition,epi in defs.items():
  for patient in sorted(set(patients[tumor])):
   ix=np.flatnonzero(tumor&epi&(patients==patient))
   cov=dict(dataset=ds,patient=ds+':'+patient,definition=definition,cells=len(ix),eligible20=len(ix)>=20,eligible50=len(ix)>=50)
   coverage.append(cov)
   if len(ix)<20:continue
   den=float(umi[ix].sum());vals={g:float(np.log1p(v[ix].sum()/den*1e4)) for g,v in small.items()}
   for name,gs in {**S['panels'],**S['outcomes']}.items():
    vv=[vals[g] for g in gs if g in vals];vals[name]=float(np.mean(vv)) if len(vv)>=.8*len(gs) else np.nan;vals[name+'_coverage']=len(vv)
   profiles.append(dict(cov,site='primary_CRC',depth=float(np.log1p(umi[ix]).mean()),total_UMI=den,CGNL1_positive_cells=int((small['CGNL1'][ix]>0).sum()),**vals))
 qa[ds]=dict(status='extracted',matrix_cells=n,annotation_cells=len(m),genes=len(names),raw_nonnegative_integer=True,tumor_patients=len(set(patients[tumor])),eligible_patients=sum(x['dataset']==ds and x['definition']=='main' and x['eligible20'] for x in coverage))
 print('STREAM_COMPLETE',ds,json.dumps(qa[ds]),flush=True)
pd.DataFrame(profiles).to_csv(O/'primary_patient_profiles.tsv',sep='\t',index=False)
pd.DataFrame(coverage).to_csv(O/'primary_cell_coverage.tsv',sep='\t',index=False)
(R/'PRIMARY_EXTRACTION_QA.json').write_text(json.dumps(qa,indent=2))
assert profiles,'No new primary cohort could be extracted'
print('PRIMARY_EXTRACTION_COMPLETE',flush=True)
