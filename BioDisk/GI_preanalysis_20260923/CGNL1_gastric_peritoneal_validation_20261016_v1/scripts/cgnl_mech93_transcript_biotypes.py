from pathlib import Path
import json,datetime
O=Path.cwd()/'mechanism_investigation_v2';a=json.loads((O/'metadata/REFERENCE_HEADER_AUDIT_77.json').read_text());rows=[]
for r in a['target_paralog_headers']:
 h=r['header'].lstrip('>').split('|');rows.append({'transcript_id':h[0],'gene_id':h[1],'gene_symbol':h[5],'length':int(h[6]),'biotype':h[7]})
coding=[r['transcript_id'] for r in rows if r['gene_symbol']=='CGNL1' and r['biotype']=='protein_coding'];noncoding=[r for r in rows if r['gene_symbol']=='CGNL1' and r['biotype']!='protein_coding'];out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reference':'officialMD5verifiedGENCODE47FASTA','transcripts':rows,'CGNL1_coding_transcripts':coding,'CGNL1_other_biotypes':noncoding,'qualifier':'codingRNAabundancedoesnotestablishproteinabundanceorlocalization','status':'prospectiveannotationqualificationbeforeDRAoutcomes'};(O/'metadata/TARGET_TRANSCRIPT_BIOTYPES_93.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
