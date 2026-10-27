from common import *
from scipy import stats
D=Q/'F04_MARKER_SENS';D.mkdir(exist_ok=True);old=pathlib.Path(C['source_run']);samples=json.loads((I/'SAMPLES.json').read_text());panels={'endothelial_junction_no_label_markers':['CLDN5','ESAM','TJP1','OCLN']};qc=[];effects=[];availability=[]
for s in samples:
 if s['sample'] not in ['GSM7290775','GSM8714595','GSM8714597']:continue
 name=s['sample'];X,sy,bc,lib=load_sample(s);ann=pd.read_csv(old/'cells'/(name+'.tsv.gz'),sep='\t').set_index('barcode').reindex(bc);assert ann.label.notna().all();ng=X.getnnz(axis=0);feat=np.c_[np.log1p(lib),np.log1p(ng)];unique=~pd.Index(sy).duplicated(keep=False);idx={g:i for i,g in enumerate(sy) if unique[i]}
 for label in ['endothelial']:
  k=np.flatnonzero(ann.label.to_numpy()==label)
  if len(k)<50:continue
  for g in C['final_targets']+C['partner_genes']:
   if g in idx:
    z=X[idx[g],k].toarray().ravel();availability.append({'sample':name,'dataset':s['dataset'],'donor':s['patient_label_provisional'],'lineage':label,'gene':g,'cells':len(k),'mean_CP10K':float((z/lib[k]*1e4).mean()),'detected_pct':float((z>0).mean()*100)})
  for target in ['CGNL1']:
   if target not in idx:continue
   z=X[idx[target]].toarray().ravel();pos=k[z[k]>0];neg=k[z[k]==0];used=set();pairs=[]
   for ip in sorted(pos,key=lambda j:(z[j],j)):
    cand=np.array([j for j in neg if j not in used],dtype=int)
    if len(cand)==0:break
    delta=np.abs(feat[cand]-feat[ip]);ok=(delta<=.25).all(axis=1)
    if not ok.any():continue
    cand=cand[ok];delta=delta[ok];jn=cand[np.argmin((delta**2).sum(axis=1))];used.add(jn);pairs.append((ip,jn))
   r={'sample':name,'dataset':s['dataset'],'donor':s['patient_label_provisional'],'lineage':label,'target':target,'cells':len(k),'positive_cells':len(pos),'positive_pct':len(pos)/len(k)*100,'matched_pairs':len(pairs),'eligible':len(pairs)>=10};qc.append(r)
   if len(pairs)<10:continue
   a,b=np.array(pairs).T
   for panel,genes in panels.items():
    ii=[idx[g] for g in genes if g in idx and g!=target]
    if len(ii)<len(genes)*.7:continue
    pp=np.r_[a,b];value=np.log1p(X[ii][:,pp].toarray()/lib[pp]*1e4).mean(axis=0);delta=value[:len(a)]-value[len(a):];effects.append({**r,'panel':panel,'effect':float(delta.mean()),'abs_logUMI_match_difference':float(np.abs(feat[a,0]-feat[b,0]).mean()),'abs_logGene_match_difference':float(np.abs(feat[a,1]-feat[b,1]).mean())})
 pd.DataFrame(qc).to_csv(D/'eligibility.tsv',sep='\t',index=False);print(name,'done',flush=True)
pd.DataFrame(availability).to_csv(D/'partner_availability.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(D/'sample_state_effects.tsv',sep='\t',index=False);e=pd.DataFrame(effects);don=e.groupby(['dataset','donor','lineage','target','panel'],as_index=False).effect.mean();don.to_csv(D/'donor_state_effects.tsv',sep='\t',index=False);rows=[]
for (lineage,target,panel),v in don.groupby(['lineage','target','panel']):
 x=v.effect.to_numpy();n=len(x);r={'lineage':lineage,'target':target,'panel':panel,'n_donors':n,'positive':int((x>0).sum()),'mean':float(x.mean()),'median':float(np.median(x)),'P':np.nan,'FDR':np.nan}
 if n>=10:r['P']=float(stats.ttest_1samp(x,0).pvalue)
 rows.append(r)
ss=pd.DataFrame(rows)
for (lineage,target),z in ss.groupby(['lineage','target']):
 q=bh(z.P.to_numpy(),len(panels));q[z.P.isna().to_numpy()]=np.nan;ss.loc[z.index,'FDR']=q
ss.to_csv(D/'summary.tsv',sep='\t',index=False)
av=pd.DataFrame(availability).groupby(['dataset','donor','lineage','gene'],as_index=False).agg(mean_CP10K=('mean_CP10K','mean'),detected_pct=('detected_pct','mean'));av.to_csv(D/'donor_partner_availability.tsv',sep='\t',index=False)
dump(D/'RESULT.json',{'status':'descriptive_source_and_state_completed','state_rows':len(ss),'eligible_donor_counts':ss.groupby(['lineage','target']).n_donors.max().to_dict() if False else ss[['lineage','target','n_donors']].drop_duplicates().to_dict('records'),'scope':'Within-lineage detection-associated state; not molecular binding or proteolytic activity. Existing donor labels provisional.'})
