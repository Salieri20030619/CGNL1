from pathlib import Path
import csv,json,hashlib,re,math
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1/C02_new_sources_20261024')
f=P/'Set04_Phospho/evidence.txt';receipt=f.with_name(f.name+'.receipt.json')
if not (f.exists() and receipt.exists()):raise SystemExit('Pending complete evidence and CRC receipt; no inference')
rec=json.loads(receipt.read_text());assert rec['crc_pass'] and hashlib.sha256(f.read_bytes()).hexdigest()==rec['sha256']
fasta=P.parents[1]/'raw/CGNL1_Q0VF96_380.fasta';seq=''.join(x.strip() for x in fasta.read_text().splitlines() if not x.startswith('>'))
sites={x['Positions within proteins']:x for x in json.loads((P/'Set04_CGNL1_rows.json').read_text())['rows'] if x['Positions within proteins'] in ['283','486']}
allcg=[];records={}
with f.open() as h:
 for row in csv.DictReader(h,delimiter='\t'):
  if 'Q0VF96' in row['Proteins'].split(';') or 'CGNL1' in row['Gene names'].split(';'):allcg.append(row);records[row['id']]=row
(P/'Set04_CGNL1_EVIDENCE_ROWS.json').write_text(json.dumps(allcg,indent=2),encoding='utf-8')
def parseprob(s):
 aa=[];probs=[];i=0
 while i<len(s):
  assert 'A'<=s[i]<='Z';aa.append(s[i]);i+=1;p=0.
  if i<len(s) and s[i]=='(':
   end=s.index(')',i);p=float(s[i+1:end]);i=end+1
  probs.append(p)
 tokens=re.findall(r'([A-Z])(?:\(([0-9.eE+-]+)\))?',s)
 assert ''.join(a for a,_ in tokens)==''.join(aa)
 assert [float(p) if p else 0. for _,p in tokens]==probs
 return ''.join(aa),probs
results={}
for position,site in sites.items():
 wanted=site['Evidence IDs'].split(';');missing=[x for x in wanted if x not in records];linked=[]
 for ident in wanted:
  if ident not in records:continue
  x=records[ident];peptide=x['Sequence'];starts=[m.start() for m in re.finditer('(?='+re.escape(peptide)+')',seq)];aa,prob=parseprob(x['Phospho (STY) Probabilities']) if x['Phospho (STY) Probabilities'] else ('',[])
  assert not aa or aa==peptide
  local=int(position)-starts[0]-1 if len(starts)==1 else -1
  lp=prob[local] if 0<=local<len(prob) else None
  reciprocal=site['id'] in x['Phospho (STY) site IDs'].split(';')
  direct=x['Type'] in ['MSMS','MULTI-MSMS'] and int(x['MS/MS count'] or 0)>0 and bool(x['MS/MS IDs'])
  mono=int(x['Phospho (STY)'] or 0)==1
  reporters=[float(x[f'Reporter intensity corrected {i}'] or 0) for i in range(1,17)]
  qualifies=reciprocal and direct and mono and lp is not None and lp>=.75 and not x['Reverse'] and not x['Potential contaminant']
  linked.append({'id':ident,'sequence':peptide,'modified_sequence':x['Modified sequence'],'probability_text':x['Phospho (STY) Probabilities'],'canonical_start':starts[0]+1 if len(starts)==1 else None,'local_probability':lp,'reciprocal':reciprocal,'direct':direct,'mono':mono,'qualifies':qualifies,'type':x['Type'],'raw_file':x['Raw file'],'MSMS_IDs':x['MS/MS IDs'],'PEP':x['PEP'],'reporters':reporters})
 mono=[x for x in linked if x['mono'] and x['reciprocal']]
 sums=[sum(x['reporters'][i] for x in mono) for i in range(16)];source=[float(site[f'Reporter intensity corrected {i}___1'] or 0) for i in range(1,17)];errors=[abs(a-b)/max(abs(b),1) for a,b in zip(sums,source)]
 allqualified=bool(mono) and all(x['qualifies'] for x in mono if any(v>0 for v in x['reporters']))
 results[position]={'missing_evidence_IDs':missing,'records':linked,'mono_contributors':len(mono),'mono_sum_relative_error_max':max(errors),'all_mono_contributors_qualified':allqualified,'site_row_mono_reporters':source,'pass':not missing and allqualified and max(errors)<1e-4}
out={'evidence_sha256':rec['sha256'],'canonical_fasta_sha256':hashlib.sha256(fasta.read_bytes()).hexdigest(),'sites':results,'pass':all(x['pass'] for x in results.values()),'probability_parser_two_implementations_agree':True,'no_biological_tests':True}
(P/'Set04_EVIDENCE_GATE.json').write_text(json.dumps(out,indent=2,allow_nan=False),encoding='utf-8')
print(json.dumps({'pass':out['pass'],'sites':{k:{a:b for a,b in v.items() if a not in ['records','site_row_mono_reporters']} for k,v in results.items()}}))
