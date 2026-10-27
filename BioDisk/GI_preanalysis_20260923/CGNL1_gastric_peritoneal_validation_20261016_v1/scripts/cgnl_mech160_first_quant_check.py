from pathlib import Path
import json,datetime,subprocess
O=Path.cwd()/'mechanism_investigation_v2';out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':None,'runs':[]}
p=O/'metadata/DRA_QUANTIFICATION_STATUS.json'
if p.exists():out['status']=json.loads(p.read_text())
for p in (O/'processed/DRA014736_quant').iterdir():
 if not p.is_dir():continue
 r={'run':p.name,'quant_sf':(p/'quant.sf').exists(),'gene_quant':(p/'gene_quant.tsv.gz').exists()};m=p/'aux_info/meta_info.json'
 if m.exists():
  meta=json.loads(m.read_text());r['metrics']={k:meta.get(k) for k in ['salmon_version','num_processed','num_mapped','percent_mapped','library_types','num_decoy_fragments','num_valid_targets']}
 log=O/'logs'/('salmon_'+p.name+'.log')
 if log.exists():r['log_tail']=log.read_text(encoding='utf-8')[-1600:]
 out['runs'].append(r)
ps=subprocess.run(['ps','-u','wangmingyuan','-o','pid,etime,pcpu,pmem,args'],capture_output=True,text=True).stdout;out['quant_processes']=[l for l in ps.splitlines() if str(O/'software/salmon_1_10_3/bin/salmon') in l and ' quant ' in l]
(O/'metadata/FIRST_QUANT_CHECK_160.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out)[:12000])
