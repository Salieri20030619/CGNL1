from pathlib import Path
import gzip,json,hashlib,xml.etree.ElementTree as ET
O=Path.cwd()/'mechanism_investigation_v2';p=O/'reference/gencode_v47/gencode.v47.transcripts.fa.gz';h=hashlib.md5(p.read_bytes()).hexdigest();official={l.split()[1].lstrip('*').removeprefix('./'):l.split()[0] for l in (p.parent/'MD5SUMS').read_text().splitlines() if len(l.split())==2};assert h==official[p.name];hits=[];first=None;seen=set();duplicates=[]
with gzip.open(p,'rt',encoding='utf-8') as f:
 for line in f:
  if not line.startswith('>'):continue
  a=line[1:].strip().split('|');assert len(a)>=7
  if first is None:first=a
  if a[0] in seen:duplicates.append(a[0])
  seen.add(a[0])
  if a[1].split('.')[0] in ['ENSG00000128849','ENSG00000143375']:hits.append({'transcript_id':a[0],'gene_id':a[1],'gene_symbol':a[5],'declared_length':a[6],'header':line.strip()})
assert not duplicates and any(r['gene_symbol']=='CGNL1' for r in hits) and any(r['gene_symbol']=='CGN' for r in hits)
out={'reference_MD5_verified':True,'first_header_fields':first,'transcripts':len(seen),'target_paralog_headers':hits,'duplicate_transcript_ids':duplicates,'pipeline_header_gene_symbol_field5_confirmed':True};(O/'metadata/REFERENCE_HEADER_AUDIT_77.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print('REFERENCE',len(seen),'target/paralog',len(hits),'header',first)
root=ET.fromstring((O/'literature/PMC10202830.xml').read_bytes());pars=[]
for e in root.iter('p'):
 t=' '.join(e.itertext())
 if 'CGNL1' in t and any(w in t.lower() for w in ['stiffness','nm2a','nm2b','tortuosity']):pars.append(t)
(O/'literature/CGNL1_MECHANICAL_SPECIFICITY_77.json').write_text(json.dumps(pars,ensure_ascii=False,indent=2),encoding='utf-8');print('FUNCTION_PARAGRAPHS','\n'.join(pars)[:14000])
