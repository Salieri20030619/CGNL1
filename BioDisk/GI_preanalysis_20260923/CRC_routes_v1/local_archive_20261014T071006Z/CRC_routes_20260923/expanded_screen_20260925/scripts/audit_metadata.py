from pathlib import Path
import gzip,json,pandas as pd
A=Path(__file__).resolve().parent.parent;M=A/'metadata';out={}
for g in ['GSE131418','GSE41258','GSE203612','GSE234804','GSE315534']:
 rows=[];s=None
 with gzip.open(M/(g+'.soft.gz'),'rt',errors='replace') as f:
  for l in f:
   if l.startswith('^SAMPLE = '):s={'GSM':l.strip().split(' = ')[1]};rows.append(s)
   if s is None:continue
   if l.startswith('!Sample_title = '):s['title']=l.strip().split(' = ')[1]
   if l.startswith('!Sample_source_name_ch1 = '):s['source']=l.strip().split(' = ')[1]
   if l.startswith('!Sample_characteristics_ch1 = '):
    v=l.strip().split(' = ',1)[1]
    if ': ' in v:k,val=v.split(': ',1);s[k.lower()]=val
 m=pd.DataFrame(rows);m.to_csv(M/(g+'_sample_metadata.tsv'),sep='\t',index=False)
 sel=m.apply(lambda r:any('liver' in str(x).lower() and ('met' in str(x).lower() or str(x).lower()=='liver') for x in r),axis=1);z=m[sel];out[g]={'samples_total':len(m),'columns':list(m.columns),'liver_labeled_samples':len(z),'liver_rows':z.astype(object).where(pd.notna(z),None).to_dict('records')}
 print(g,len(m),'liver_labels',len(z),'fields',list(m.columns),flush=True)
 if len(z):print(z.head(2).to_string(index=False),flush=True)
for acc in ['13651','13652','13655']:
 m=pd.read_csv(M/('E-MTAB-'+acc+'.sdrf.txt'),sep='\t');cols=[c for c in m if any(k in c.lower() for k in ['individual','source name','image','file'])];out['E-MTAB-'+acc]={'rows':len(m),'unique_fields':{c:m[c].dropna().astype(str).unique().tolist() for c in cols if 'fastq' not in c.lower()}};print('E-MTAB-'+acc,{c:m[c].nunique() for c in cols},flush=True)
(A/'EXTENDED_METADATA_AUDIT.json').write_text(json.dumps(out,indent=2))
