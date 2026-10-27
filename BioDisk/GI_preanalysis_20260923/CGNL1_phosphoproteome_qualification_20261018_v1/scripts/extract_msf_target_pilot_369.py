from pathlib import Path
import json
from msf_range_vfs_368 import HTTPReader,connect,b
pilot=json.loads((b/'metadata/MSF_PILOT_RANGE_366.json').read_text());size=int(pilot['content_range'].split('/')[-1]);reader=HTTPReader(pilot['url'],size,b/'raw/MSF_pilot_range_cache_368',cap=64*1024*1024);v,c=connect('target369',b/'raw/MSF_pilot_header_366.bin',reader,size)
queries=[]
def query(sql,args=()):
 cur=c.cursor();rows=list(cur.execute(sql,args));cols=[x[0] for x in cur.get_description()] if rows else []
 return [dict(zip(cols,r)) for r in rows]
# APSW description must be captured before statement exhausts
def query(sql,args=()):
 cur=c.cursor();cur.execute(sql,args)
 try:cols=[x[0] for x in cur.get_description()]
 except Exception:cols=[]
 return [dict(zip(cols,r)) for r in cur]
proteins=query("SELECT UniqueSequenceID,Accession,Description,FastaTitleLines,Sequence FROM TargetProteins WHERE Description LIKE '%cingulin%' OR FastaTitleLines LIKE '%CGNL1%'")
print('target/comparator proteins',[(r['Accession'],r['Description']) for r in proteins],flush=True)
psms=[]
for r in proteins:
 sql='SELECT p.WorkflowID,p.PeptideID,p.Sequence,p.ModifiedSequence,p.Modifications,p.ParentProteinCount,p.ParentProteinAccessions,p.MatchConfidence,p.PercolatorqValue,p.PercolatorPEP,p.PrecursorAbundance,p.HasReliableFeature,p.Charge,p.RetentionTime FROM TargetProteinsTargetPsms l JOIN TargetPsms p ON p.WorkflowID=l.TargetPsmsWorkflowID AND p.PeptideID=l.TargetPsmsPeptideID WHERE l.TargetProteinsUniqueSequenceID=?'
 queries.append({'sql':sql,'protein':r['UniqueSequenceID'],'query_plan':list(c.execute('EXPLAIN QUERY PLAN '+sql,(r['UniqueSequenceID'],)))})
 for row in query(sql,(r['UniqueSequenceID'],)):row['protein_record']=r['UniqueSequenceID'];psms.append(row)
meta=query('SELECT CustomDataID,WorkflowID,ProcessingNodeName,Name,CustomDataType,FormatString FROM ProcessingNodeCustomData');mods=query('SELECT * FROM FoundModifications')
out={'source':pilot['url'],'proteins':proteins,'target_and_comparator_psms':psms,'custom_data_metadata':meta,'modifications':mods,'queries':queries,'fetched_bytes':reader.total,'requests':reader.nreq,'range_receipts':reader.receipts,'whole_file_integrity_not_claimed':True,'site_localization_and_normalization_pending':True};(b/'metadata/MSF_PILOT_TARGET_369.json').write_text(json.dumps(out,indent=2),encoding='utf-8');c.close();v.unregister();print('PSMs',len(psms),'bytes',reader.total,'requests',reader.nreq,flush=True)
