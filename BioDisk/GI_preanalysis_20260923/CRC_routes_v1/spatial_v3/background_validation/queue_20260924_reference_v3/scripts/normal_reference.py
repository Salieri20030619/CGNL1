from common import *
D=Q/'NORMAL_REF';D.mkdir(exist_ok=True);R=pathlib.Path(C['normal_reference_input'])
manifest=json.loads((R/'MANIFEST.json').read_text())
for r in manifest:
 p=R/r['file'];assert p.stat().st_size==r['bytes'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256']
anno=pd.read_csv(R/'GSE115469_CellClusterType.txt.gz',sep='\t').set_index('CellName');genes=[];chunks=[];cols=None
for t in pd.read_csv(R/'GSE115469_Data.csv.gz',index_col=0,chunksize=256):
 if cols is None:cols=t.columns;assert cols.is_unique and set(cols)==set(anno.index);anno=anno.reindex(cols)
 assert t.columns.equals(cols);z=t.to_numpy(float);assert np.isfinite(z).all() and z.min()>=0 and z.max()<40
 chunks.append(sparse.csr_matrix(np.expm1(z*np.log(2))));genes.extend(t.index.astype(str))
X=sparse.vstack(chunks,format='csr');del chunks;sy=np.array(genes);assert pd.Index(sy).is_unique
lib=np.asarray(X.sum(axis=0)).ravel();ng=X.getnnz(axis=0);mt=np.asarray(X[np.char.startswith(sy,'MT-')].sum(axis=0)).ravel()/lib*100
good=(ng>=200)&(mt<=25)&(lib>0);Y=X.multiply(1e6/lib).tocsr();author=anno.CellType.to_numpy();sample=anno.Sample.to_numpy()
assert set(sample)=={'P1TLH','P2TLH','P3TLH','P4TLH','P5TLH'}
panels=C['lineage_markers'];scores=[];hits=[]
for gs in panels.values():
 z=Y[np.isin(sy,gs)].toarray()/100;scores.append(np.log1p(z).mean(axis=0));hits.append((z>0).sum(axis=0))
scores=np.array(scores).T;hits=np.array(hits).T;best=scores.argmax(axis=1);names=np.array(list(panels));margin=np.sort(scores,axis=1)[:,-1]-np.sort(scores,axis=1)[:,-2]
hep=good&np.char.startswith(author.astype(str),'Hepatocyte')&(names[best]=='hepatocyte')&(hits[:,list(panels).index('hepatocyte')]>=2)&(margin>=.4)
cell=anno.reset_index();cell['QC_pass']=good;cell['mitochondrial_pct_from_linearized_expression']=mt;cell['genes_detected']=ng;cell['reference_hepatocyte']=hep;cell['lineage_margin']=margin;cell.to_csv(D/'cell_annotation_QC.tsv.gz',sep='\t',index=False)
meta=[];source=[];profiles=[];donors=[]
for donor in sorted(set(sample)):
 k=hep&(sample==donor);meta.append({'sample':donor,'author_hepatocytes':int(((sample==donor)&np.char.startswith(author.astype(str),'Hepatocyte')).sum()),'QC_marker_hepatocytes':int(k.sum()),'eligible':int(k.sum())>=50})
 if k.sum()>=50:donors.append('GSE115469:'+donor);profiles.append(np.asarray(Y[:,k].mean(axis=1)).ravel())
 for label in sorted(set(author)):
  k=good&(sample==donor)&(author==label)
  if k.sum()==0:continue
  for g in C['targets']:
   ix=np.flatnonzero(sy==g)
   if len(ix)==1:
    v=Y[ix[0],k].toarray().ravel();source.append({'donor':donor,'label':label,'gene':g,'cells':int(k.sum()),'mean_CP10K':float(v.mean()/100),'detected_pct':float((v>0).mean()*100)})
pd.DataFrame(meta).to_csv(D/'donor_QC.tsv',sep='\t',index=False);pd.DataFrame(source).to_csv(D/'target_source_expression.tsv',sep='\t',index=False)
assert len(donors)>=3,'fewer than3 qualifying hepatocyte donors'
np.savez_compressed(D/'hepatocyte_profiles.npz',genes=sy.astype(str),donors=np.array(donors),CPM=np.stack(profiles,axis=1))
dump(D/'RESULT.json',{'status':'processed_reference_ready','input_cells':len(cols),'QC_cells':int(good.sum()),'genes':len(sy),'eligible_hepatocyte_donors':donors,'donor_QC':meta,'normalization':C['normal_reference_scale'],'scope':'Author-filtered, normalized values; no raw UMI QC/count model. Hepatocyte donor means equally weighted; >=50 cells per donor. Descriptive normal tissue source.'})
print(json.dumps(meta),flush=True)
