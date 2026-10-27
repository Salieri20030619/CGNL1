from pathlib import Path
import json,re,math,datetime,zipfile,xml.etree.ElementTree as E,posixpath
import numpy as np
from scipy.stats import rankdata,t as tdist
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');D=B/'results/core_qualification_568';O=B/'results/HCC_exploratory_572';O.mkdir(exist_ok=True);spec=json.loads((B/'metadata/HCC_EXPLORATORY_SPEC571.json').read_text());ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
def col(x):return re.sub(r'\d+$','',x)
def read_matrix(record,sheetname):
 s=next(x for x in record['sheets'] if x['name']==sheetname);headers=s['headers'][1]['cells'];columns={col(k):v for k,v in headers.items() if isinstance(v,str) and re.fullmatch(r'\d+[TN]',v)};labels=list(columns.values());assert len(labels)==len(set(labels));idx={k:i for i,k in enumerate(columns)};matrix=np.full((s['rows']+1,len(labels)),np.nan)
 with zipfile.ZipFile(record['source']) as z:
  relations={x.get('Id'):x.get('Target') for x in E.fromstring(z.read('xl/_rels/workbook.xml.rels'))};sheets=E.fromstring(z.read('xl/workbook.xml')).findall('m:sheets/m:sheet',ns);sh=next(x for x in sheets if x.get('name')==sheetname);target=relations[sh.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')];member=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
  with z.open(member) as f:
   for event,el in E.iterparse(f,events=('end',)):
    if not el.tag.endswith('}row'):continue
    rn=int(el.get('r'));assert rn<matrix.shape[0]
    if rn>=3:
     for c in el.findall('m:c',ns):
      cc=col(c.get('r'));v=c.find('m:v',ns)
      if cc in idx and v is not None and c.get('t') not in ['s','str','e','b']:
       matrix[rn,idx[cc]]=float(v.text)
    el.clear()
 raw=matrix[3:];valid=np.isfinite(matrix)&(matrix>0);logged=np.full(matrix.shape,np.nan);logged[valid]=np.log2(matrix[valid]);med=np.nanmedian(logged[3:],axis=0);assert np.isfinite(med).all()
 qa=dict(samples=len(labels),matrix_rows=s['rows']-2,positive=int((np.isfinite(raw)&(raw>0)).sum()),zero=int((raw==0).sum()),negative=int((raw<0).sum()),missing=int(np.isnan(raw).sum()),infinite=int(np.isinf(raw).sum()))
 return s,labels,logged,med,qa
pr=json.loads((D/'HCC_protein.json').read_text());sr=json.loads((D/'HCC_phospho.json').read_text());ps,pl,pm,pmed,pqa=read_matrix(pr,'1. 6512 proteins');ss,sl,sm,smed,sqa=read_matrix(sr,'1. phosphoproteomic');pi={x:i for i,x in enumerate(pl)};si={x:i for i,x in enumerate(sl)};patients=sorted({x[:-1] for x in sl if x.endswith('T')} & {x[:-1] for x in sl if x.endswith('N')} & {x[:-1] for x in pl if x.endswith('T')} & {x[:-1] for x in pl if x.endswith('N')});proteinrow=ps['target_rows'][0]['row'];fasta=''.join(x.strip() for x in (B/'raw/CGNL1_Q0VF96_380.fasta').read_text().splitlines() if not x.startswith('>'))
peptides=[]
for row in ss['target_rows']:
 v=list(row['cells'].values());start=int(v[2]);sequence=v[3];plain=re.sub(r'\[[^\]]+\]','',sequence);assert fasta[start-1:start-1+len(plain)]==plain;pos=0;sites=[]
 for token in re.finditer(r'([A-Z])(\[[^\]]+\])?',sequence):
  if token.group(2) and 'Phospho' in token.group(2):sites.append(token.group(1)+str(start+pos))
  pos+=1
 if len(sites)==1 and sites[0] in spec['family']:peptides.append(dict(row=row['row'],sequence=sequence,site=sites[0]))
def signed_rank(values):
 d=np.asarray(values,dtype=float);d=d[d!=0]
 if len(d)==0:return None
 ranks=np.rint(rankdata(abs(d))*2).astype(int);total=int(ranks.sum());observed=int(ranks[d>0].sum());small=min(observed,total-observed);counts=[0]*(total+1);counts[0]=1;used=0
 for rank in ranks:
  rank=int(rank)
  for k in range(used,-1,-1):counts[k+rank]+=counts[k]
  used+=rank
 return min(1.,2*sum(counts[:small+1])/2**len(d))
rows=[];summaries=[]
for site in spec['family']:
 pp=[x for x in peptides if x['site']==site];values=[]
 for patient in patients:
  a,b=pi[patient+'T'],pi[patient+'N'];pa,pb=pm[proteinrow,a],pm[proteinrow,b]
  if not (np.isfinite(pa) and np.isfinite(pb)):continue
  dp_raw=pa-pb;dp=dp_raw-(pmed[a]-pmed[b]);a,b=si[patient+'T'],si[patient+'N'];byseq={};raw_byseq={};complete_rows=[]
  for pep in pp:
   ta,nb=sm[pep['row'],a],sm[pep['row'],b]
   if np.isfinite(ta) and np.isfinite(nb):
    delta=ta-nb;byseq.setdefault(pep['sequence'],[]).append(delta-(smed[a]-smed[b]));raw_byseq.setdefault(pep['sequence'],[]).append(delta);complete_rows.append(pep['row'])
  if not byseq:continue
  ds=float(np.median([np.median(x) for x in byseq.values()]));ds_raw=float(np.median([np.median(x) for x in raw_byseq.values()]));row=dict(site=site,patient=patient,delta_site=ds,delta_protein=float(dp),adjusted_delta=ds-dp,uncentered_adjusted_delta=ds_raw-dp_raw,complete_peptide_rows=complete_rows,unique_peptide_sequences=len(byseq));rows.append(row);values.append(row)
 d=np.array([x['adjusted_delta'] for x in values]);raw=np.array([x['uncentered_adjusted_delta'] for x in values]);n=len(d);summary=dict(site=site,n=n,positive=int((d>0).sum()),negative=int((d<0).sum()),mean=float(d.mean()) if n else None,median=float(np.median(d)) if n else None,P=signed_rank(d) if n>=10 else None,FDR=None,mean_delta_site=float(np.mean([x['delta_site'] for x in values])) if n else None,mean_delta_protein=float(np.mean([x['delta_protein'] for x in values])) if n else None,uncentered_mean=float(raw.mean()) if n else None,uncentered_P=signed_rank(raw) if n>=10 else None)
 if n>=10:
  xx=np.array([x['delta_protein'] for x in values]);yy=np.array([x['delta_site'] for x in values]);X=np.column_stack([np.ones(n),xx]);inv=np.linalg.inv(X.T@X);coef=inv@X.T@yy;err=yy-X@coef;hat=np.einsum('ij,jk,ik->i',X,inv,X);cov=inv@(X.T@((err/(1-hat))[:,None]**2*X))@inv;se=np.sqrt(np.diag(cov));summary['regression_sensitivity']=dict(intercept=float(coef[0]),protein_slope=float(coef[1]),intercept_HC3_SE=float(se[0]),intercept_P=float(2*tdist.sf(abs(coef[0]/se[0]),n-2)),protein_delta_min=float(xx.min()),protein_delta_max=float(xx.max()))
 if site=='S283':
  common=[x for x in values if len(x['complete_peptide_rows'])==len(pp)];cd=np.array([x['adjusted_delta'] for x in common]);summary['all_peptide_rows_complete_sensitivity']=dict(n=len(cd),mean=float(cd.mean()) if len(cd) else None,P=signed_rank(cd) if len(cd)>=10 else None)
 summaries.append(summary)
eligible=[x for x in summaries if x['P'] is not None];order=sorted(eligible,key=lambda x:x['P']);last=1.
for rank in range(len(order),0,-1):last=min(last,order[rank-1]['P']*len(order)/rank);order[rank-1]['FDR']=last
(O/'patient_site_values572.json').write_text(json.dumps(rows,indent=2),encoding='utf-8');(O/'NORMALIZATION572.json').write_text(json.dumps(dict(protein_QC=pqa,phospho_QC=sqa,protein_medians=dict(zip(pl,pmed.tolist())),phospho_medians=dict(zip(sl,smed.tolist()))),indent=2),encoding='utf-8');result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope=spec['scope'],family_size=len(eligible),source_qualification_pending=True,mainline_validated=False,results=summaries);(O/'RESULT572.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)
