from pathlib import Path
import csv,json,re,math,statistics,hashlib,datetime,collections
import numpy as np
import openpyxl
from scipy.stats import rankdata,wilcoxon
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
D=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576/H01')
O=B/'results/H01_582';O.mkdir(exist_ok=True)
def save(name,x): (O/name).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
sites=['S112','S198','S202','S261','S283','S296','S297','S349','S414','S486','S514','S708']
save('SPEC.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),family=sites,primary='Paired tumour minus adjacent normal phosphosite log ratio minus paired total CGNL1 log ratio; each omics independently sample-median centered; no imputation',sensitivity='As-supplied log ratios without additional sample centering',n='Independent clinically mapped donor pairs; n<10 descriptive',test='Two-sided exact Wilcoxon signed rank conditional sign distribution; zeros excluded; average tied ranks; BH within evaluable prespecified site family',source_limit='PDC reprocessed matrix; 160 clinical pairs; case 109 absent from both quantitative matrices leaves 159 eligible pairs; correspondence to publication exclusion not independently established; effects in supplied log-ratio units pending log-base source verification'))
w=openpyxl.load_workbook(next(D.glob('*.xlsx')),read_only=True,data_only=True)
clinical=list(w['Table 1. Clinical Information'].values)[1:]
assert len({r[0] for r in clinical})==len(clinical)
assert all(r[7]=='HCC' for r in clinical)
pdc=json.loads((B/'metadata/bulk_576/PDC000198_biospecimens.json').read_text())['biospecimenPerStudy']
mp={r['sample_submitter_id']:r for r in pdc}
pairs=[]
for r in clinical:
    t,n=str(int(r[3])),str(int(r[4]))
    assert mp[t]['case_submitter_id']==mp[n]['case_submitter_id']
    assert mp[t]['sample_type']=='Primary Tumor' and mp[n]['sample_type']=='Solid Tissue Normal'
    pairs.append(dict(case=str(r[0]),tumor='T'+t,normal='P'+n))
assert len({p['tumor'] for p in pairs})==len(pairs)==len({p['normal'] for p in pairs})
save('PAIRING.json',pairs)
seq=''.join(l.strip() for l in (B/'raw/CGNL1_Q0VF96_380.fasta').read_text().splitlines() if not l.startswith('>'))
def val(s):
    if s in ('','NA','NaN','nan','N/A'):return np.nan
    x=float(s)
    return x if math.isfinite(x) else np.nan
def read(p,phos=False):
    with p.open(encoding='utf-8') as f:
        r=csv.DictReader(f,delimiter='\t');heads=r.fieldnames
        cols=[c for c in heads if re.fullmatch(r'[TP]\d+ Log Ratio',c)]
        data=[];ids=[];targets=[]
        for row in r:
            key=row[heads[0]]
            if key in ['Mean','Median','StdDev']:continue
            data.append([val(row[c]) for c in cols]);ids.append(key)
            if (row.get('Gene')=='CGNL1'):
                targets.append((len(data)-1,row))
    a=np.asarray(data,dtype=float);offset=np.nanmedian(a,axis=0)
    return cols,a,offset,targets
pc,pa,po,pt=read(D/'Zhou_Liver_Cancer_Proteome.tmt11.tsv')
sc,sa,so,st=read(D/'Zhou_Liver_Cancer_Phosphoproteome.phosphosite.tmt11.tsv',True)
assert len(pt)==1
pi=pt[0][0];pidx={c[:-10]:i for i,c in enumerate(pc)};sidx={c[:-10]:i for i,c in enumerate(sc)}
excluded=[x for x in pairs if not all(x[k] in pidx and x[k] in sidx for k in ['tumor','normal'])]
save('PAIRING_COVERAGE.json',dict(clinical_pairs=len(pairs),excluded=excluded,reason='Both omics must have tumor and normal columns; no imputation'))
pairs=[x for x in pairs if x not in excluded]
assert len(pairs)==159 and len(excluded)==1 and excluded[0]['case']=='109'
save('ANALYSIS_PAIRING.json',pairs)
save('NORMALIZATION.json',dict(protein=dict(zip(pc,po.tolist())),phosphosite=dict(zip(sc,so.tolist())),protein_rows=len(pa),phosphosite_rows=len(sa)))
sequence=[];mapped={}
for i,row in st:
    sid=row['Phosphosite'].split(':')[-1].upper()
    if not re.fullmatch('[STY]\d+',sid):continue
    checks=[]
    for peptide in row['Peptide'].split(';'):
        starts=[m.start() for m in re.finditer('(?='+re.escape(peptide.upper())+')',seq)]
        coords=[[peptide[j].upper()+str(start+j+1) for j in range(len(peptide)) if peptide[j].islower()] for start in starts]
        checks.append(dict(peptide=peptide,positions=coords,pass_mapping=[sid] in coords))
    ok=all(x['pass_mapping'] for x in checks)
    sequence.append(dict(source=row['Phosphosite'],site=sid,checks=checks,valid=ok))
    if ok: mapped.setdefault(sid,[]).append(i)
save('SEQUENCE_QA.json',sequence)
def exact_signed(x):
    y=[float(v) for v in x if v!=0]
    if not y:return None
    ranks=rankdata(np.abs(y));ranks2=[int(round(2*v)) for v in ranks]
    observed=sum(r for r,v in zip(ranks2,y) if v>0)
    total=sum(ranks2);cut=min(observed,total-observed)
    ways=[0]*(cut+1);ways[0]=1
    for r in ranks2:
        for j in range(cut,r-1,-1):ways[j]+=ways[j-r]
    return min(1.,2*sum(ways)/(2**len(y)))
results=[];values=[]
for site in sites:
    rows=mapped.get(site,[]);v=[];vraw=[];sd=[];pd=[]
    for pair in pairs:
        ti,ni=pidx[pair['tumor']],pidx[pair['normal']]
        ts,ns=sidx[pair['tumor']],sidx[pair['normal']]
        dp=pa[pi,ti]-pa[pi,ni];dpc=dp-po[ti]+po[ni]
        dsites=[sa[j,ts]-sa[j,ns] for j in rows if np.isfinite(sa[j,ts]) and np.isfinite(sa[j,ns])]
        if not dsites or not np.isfinite(dp):continue
        ds=float(np.median(dsites));dsc=ds-so[ts]+so[ns]
        a=float(dsc-dpc);raw=float(ds-dp);v.append(a);vraw.append(raw);sd.append(float(dsc));pd.append(float(dpc))
        values.append(dict(site=site,**pair,site_difference=float(dsc),protein_difference=float(dpc),adjusted=a,unscaled_adjusted=raw,source_rows=len(dsites)))
    r=dict(site=site,n=len(v),up=sum(x>0 for x in v),down=sum(x<0 for x in v),mean=float(np.mean(v)) if v else None,median=float(np.median(v)) if v else None,mean_site=float(np.mean(sd)) if sd else None,mean_protein=float(np.mean(pd)) if pd else None,p=exact_signed(v) if len(v)>=10 else None,sensitivity_mean=float(np.mean(vraw)) if vraw else None,sensitivity_p=exact_signed(vraw) if len(vraw)>=10 else None)
    results.append(r)
valid=sorted([r for r in results if r['p'] is not None],key=lambda r:r['p']);q=1.
for i in range(len(valid)-1,-1,-1):
    q=min(q,valid[i]['p']*len(valid)/(i+1));valid[i]['fdr']=q
save('PATIENT_VALUES.json',values)
save('RESULT.json',dict(cohort='H01',clinical_pairs=len(pairs),tested_family_size=len(valid),results=results,source_files=[dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in D.iterdir() if p.suffix in ['.tsv','.xlsx','.txt']],status='exploratory_complete_pending_published_n_and_scale_provenance_reconciliation'))
print(json.dumps(results))
