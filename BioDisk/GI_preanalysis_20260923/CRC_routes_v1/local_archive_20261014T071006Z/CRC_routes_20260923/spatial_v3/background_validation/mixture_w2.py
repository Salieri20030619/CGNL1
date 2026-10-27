from common import *
from mixlib import *
import tarfile
D=Q/'V03_MIXTURE';D.mkdir(exist_ok=True);gate=json.loads((Q/'V02_SOURCE/GATE.json').read_text())
if not gate['allowed']:dump(D/'RESULT.json',{'status':'skipped_insufficient_reference','reason':gate});sys.exit(0)
a=pd.read_csv(S/'results/screen/spot_target_coordinates.tsv.gz',sep='\t');a=a[~a.recurrence];masks=pd.read_csv(I/'V01/masks.tsv.gz',sep='\t');masks=masks[masks.rule=='5_of_6'];out=[];qc=[]
for name,v in a.groupby('sample',sort=True):
 v=v.reset_index(drop=True)
 try:
  with tarfile.open(S/'raw/liver'/(name+'.tar.gz')) as t:
   ns=t.getnames();read=lambda suffix:gzip.GzipFile(fileobj=t.extractfile(next(x for x in ns if x.endswith(suffix))));ft=pd.read_csv(read('filtered_feature_bc_matrix/features.tsv.gz'),sep='\t',header=None);bc=pd.read_csv(read('filtered_feature_bc_matrix/barcodes.tsv.gz'),sep='\t',header=None)[0];
   with threadpool_limits(limits=2):X=io.mmread(read('filtered_feature_bc_matrix/matrix.mtx.gz')).tocsr()
  ix=pd.Index(bc).get_indexer(v.barcode);assert (ix>=0).all();X=X[:,ix];lib=np.asarray(X.sum(axis=0)).ravel();assert np.array_equal(lib,v.UMI.to_numpy());sy=ft[1].astype(str).to_numpy();types,frac,pred,diag=fit_reference(X,sy,lib);hepi=list(types).index('hepatocyte');anchor=(v.Layer3=='Hepatic lobule').to_numpy()&v.QC_pass.to_numpy()&(v.distance_to_macro_um.to_numpy()>=600);assert anchor.sum()>=50,'too_few_distant_liver_calibration_spots';hep_fit=float(np.median(frac[anchor,hepi]));assert hep_fit>=.4,'poor_hepatic_anchor_fit';z=v[['barcode','patient','QC_pass','distance_to_hepatic_um']].copy();z=pd.concat([z,pd.DataFrame(frac,columns=['RNAweight_'+x for x in types])],axis=1);mm=masks[masks['sample']==name].set_index('barcode').reindex(v.barcode);eligible=v.QC_pass.to_numpy();macro=(v.Layer3=='Liver macrometastasis tumor').to_numpy();d=v.distance_to_hepatic_um.to_numpy();groups={'strict5':(mm.near.fillna(False).to_numpy(),mm.far.fillna(False).to_numpy()),'original':(macro&eligible&(d<=300+1e-6),macro&eligible&(d>=600-1e-6))}
  for gene,expected in pred.items():
   obs=v[gene+'_counts'].to_numpy();den=np.sum(expected[anchor]*lib[anchor]/1e6);total=obs[anchor].sum()
   if den<=0 or total<10:continue
   scale=total/den;expected=expected*scale;z[gene+'_expected_CPM']=expected;z[gene+'_residual_log2CPM']=np.log2(1+obs/np.maximum(lib,1)*1e6)-np.log2(1+expected)
   for mask,(near,far) in groups.items():
    if near.sum()<20 or far.sum()<20:continue
    def residual(k):return np.log2(1+obs[k].sum()/lib[k].sum()*1e6)-np.log2(1+np.sum(expected[k]*lib[k])/lib[k].sum())
    out.append({'sample':name,'patient':int(v.patient.iloc[0]),'gene':gene,'mask':mask,'effect':float(residual(near)-residual(far)),'near_spots':int(near.sum()),'far_spots':int(far.sum()),'calibration_scale':float(scale)})
  z.to_csv(D/(name+'_mixture.tsv.gz'),sep='\t',index=False);qc.append(dict(sample=name,status='fitted',hepatic_anchor_RNAweight=hep_fit,anchor_spots=int(anchor.sum()),**diag))
 except Exception as err:qc.append({'sample':name,'status':'not_interpretable','reason':str(err)})
 pd.DataFrame(qc).to_csv(D/'sample_QC.tsv',sep='\t',index=False);print(name,qc[-1]['status'],flush=True)
summarize(out,D)
