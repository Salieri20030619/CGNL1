from pathlib import Path
import csv,json,collections,datetime
import numpy as np
p=Path.cwd();platform={};inside=False;header=None
with (p/'raw/GPL10332_self_297.soft').open(encoding='utf-8',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row:continue
  if row[0]=='!platform_table_begin':inside=True;continue
  if row[0]=='!platform_table_end':break
  if not inside:continue
  if header is None:header=row;continue
  d=dict(zip(header,row));platform[d['ID']]=d
raw=collections.defaultdict(list);h=None;params=None;protocol=None
with (p/'processed/Agilent_GSM4191604_cp1252_to_utf8_307.txt').open(encoding='utf-8',newline='') as f:
 for row in csv.reader(f,delimiter='\t'):
  if not row:continue
  if row[0]=='FEPARAMS':params=row;continue
  if params is not None and protocol is None and row[0]=='DATA':protocol=dict(zip(params,row))
  if row[0]=='FEATURES':h=row;continue
  if h is not None and row[0]=='DATA':d=dict(zip(h,row));raw[d['ProbeName']].append(d)
a=np.load(p/'processed/matrix_input_306.npz',allow_pickle=False);i=list(a['accessions']).index('GSM4191604');values=a['values'][:,i];matched=[];ambiguous=[]
for probe,value in zip(a['probe_ids'],values):
 name=platform[str(probe)]['NAME'];rr=raw[name]
 if len(rr)!=1:ambiguous.append({'ID':str(probe),'name':name,'raw_rows':len(rr)});continue
 v=float(rr[0]['gProcessedSignal']);matched.append((float(value),v))
errors=[abs(x-y) for x,y in matched];relative=[abs(x-y)/max(1,abs(y)) for x,y in matched];target=raw['A_23_P163306'];assert len(target)==1
x={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'issue':'EGFR-Q01','sample':'GSM4191604','raw_protocol':protocol['Protocol_Name'],'raw_channel_count':protocol['Scan_NumChannels'],'green_processed_signal_present':'gProcessedSignal' in h,'red_processed_signal_present':'rProcessedSignal' in h,'matched_probe_values':len(matched),'ambiguous_mappings':ambiguous[:10],'max_absolute_difference':max(errors),'max_relative_difference':max(relative),'all_match_within_1e-6_relative':max(relative)<1e-6,'target_raw_QC':{k:target[0][k] for k in ['ProbeName','gProcessedSignal','gIsFeatNonUnifOL','gIsPosAndSignif','gIsWellAboveBG','gIsSaturated']},'interpretation':'Raw file contains single-channel GE1 green fluorescence data. Quantitative match determines deposited-value channel; Cy5/GE2 sample metadata conflicts with actual raw file.','whole_project_complete':False};(p/'metadata/RAW_CHANNEL_VERIFICATION_308.json').write_text(json.dumps(x,indent=2),encoding='utf-8');print(json.dumps(x))
