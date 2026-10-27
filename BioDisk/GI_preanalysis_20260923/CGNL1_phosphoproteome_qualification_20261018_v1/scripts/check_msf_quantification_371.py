import json
from msf_range_vfs_368 import HTTPReader,connect,b
p=json.loads((b/'metadata/MSF_PILOT_RANGE_366.json').read_text());size=int(p['content_range'].split('/')[-1]);r=HTTPReader(p['url'],size,b/'raw/MSF_pilot_range_cache_368',cap=64*1024*1024);v,c=connect('quant371',b/'raw/MSF_pilot_header_366.bin',r,size)
a=json.loads((b/'metadata/MSF_PILOT_TARGET_369.json').read_text());target=next(z for z in a['proteins'] if 'cingulin-like' in z['Description']);ps=[z for z in a['target_and_comparator_psms'] if z['protein_record']==target['UniqueSequenceID'] and z['PercolatorqValue'] is not None and z['PercolatorqValue']<=.01];links=[]
for z in ps:
 sql='SELECT f.WorkflowId,f.Id,f.Area,f.Intensity,f.ChargeState,f.ApexRT,f.CompVoltageV FROM LcmsFeaturesTargetPsms l JOIN LcmsFeatures f ON f.WorkflowId=l.LcmsFeaturesWorkflowId AND f.Id=l.LcmsFeaturesId WHERE l.TargetPsmsWorkflowID=? AND l.TargetPsmsPeptideID=?';rows=list(c.execute(sql,(z['WorkflowID'],z['PeptideID'])));links.append({'PSM':z['PeptideID'],'features':rows})
print('links',links,flush=True)
counts={t:next(c.execute('SELECT COUNT(*) FROM '+t))[0] for t in ['LcmsFeatures','LcmsFeaturesTargetPsms','TargetPsms']};print('counts',counts,flush=True)
quant={}
if counts['LcmsFeatures']:
 quant['all_feature_area_summary']=list(c.execute('SELECT count(*),sum(Area),min(Area),max(Area) FROM LcmsFeatures WHERE Area>0'))
out={'target_feature_links':links,'counts':counts,'quantification':quant,'new_range_bytes':r.total,'requests':r.nreq,'source':p['url'],'range_receipts':r.receipts};(b/'metadata/MSF_QUANTIFICATION_PILOT_371.json').write_text(json.dumps(out,indent=2),encoding='utf-8');c.close();v.unregister();print('done',r.total,flush=True)
