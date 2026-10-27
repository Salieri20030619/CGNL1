
from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
import csv,json,math,statistics,datetime,hashlib,select
from scipy.stats import wilcoxon,ttest_rel
p=Path.cwd();m=p/'mechanism_investigation_v2'
launch=json.loads((m/'metadata/RESUME_PIPELINE_LAUNCH_192.json').read_text(encoding='utf-8'));receipt=Path(launch['receipt'])
if not receipt.exists():
 try:
  fd=os.pidfd_open(launch['wrapper_pid']);poll=select.poll();poll.register(fd,select.POLLIN);poll.poll();os.close(fd)
 except ProcessLookupError:pass
assert receipt.exists() and json.loads(receipt.read_text(encoding='utf-8'))['exit_code']==0,'Full pipeline not successfully complete'
def rows(f):
 with f.open(encoding='utf-8') as h:return list(csv.DictReader(h,delimiter='	'))
def write(name,values):
 with (m/'results'/name).open('x',encoding='utf-8',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(values[0]),delimiter='	');w.writeheader();w.writerows(values)
mapping=rows(m/'reference/gencode_v47/transcript_to_gene.tsv');target={r['transcript_id'] for r in mapping if r['gene_id'].split('.')[0]=='ENSG00000128849'};assert len(target)==3
frozen=json.loads((m/'DRA_TARGET_STATISTICAL_SPECIFICATION.json').read_text());needed=set(frozen['additional_fixed_genes']+['CGNL1'])
lookup={r['transcript_id']:r['gene_symbol'] for r in mapping if r['gene_symbol'] in needed}
coding=json.loads((m/'metadata/TARGET_TRANSCRIPT_BIOTYPES_93.json').read_text())['CGNL1_coding_transcripts']
original=rows(m/'results/DRA014736_target_samples.tsv');checks=[];computed={};raw_genes={}
idx=json.loads((m/'metadata/INDEX_COMPLETE.json').read_text())['info']['SeqHash']
for r in original:
 root=m/'processed/DRA014736_quant'/r['run'];f=root/'quant.sf';tr=rows(f);names=[a['Name'] for a in tr];assert len(names)==len(set(names))
 selected={a['Name']:a for a in tr if a['Name'] in target};assert set(selected)==target
 tpm=math.fsum(float(a['TPM']) for a in selected.values());nr=math.fsum(float(a['NumReads']) for a in selected.values());total=math.fsum(float(a['NumReads']) for a in tr);ct=math.fsum(float(selected[x]['TPM']) for x in coding)
 calc={'TPM':tpm,'estimated_reads':nr,'total_estimated_gene_reads':total,'log2TPM1':math.log2(1+tpm),'log2CPM1':math.log2(1+nr/total*1e6)}
 error=max(abs(v-float(r[k])) for k,v in calc.items());assert error<1e-6,(r['run'],error)
 computed[(r['donor'],r['condition'])]=dict(calc,defined_CDS_TPM=ct)
 genes={g:[] for g in needed}
 for a in tr:
  if a['Name'] in lookup:genes[lookup[a['Name']]].append(float(a['TPM']))
 assert all(genes[g] for g in needed)
 raw_genes[(r['donor'],r['condition'])]={g:math.fsum(v) for g,v in genes.items()}
 meta=json.loads((root/'aux_info/meta_info.json').read_text());assert meta['index_seq_hash']==idx and not meta['quant_errors']
 checks.append({'run':r['run'],'max_error':error,'percent_mapped':meta['percent_mapped'],'library_types':meta['library_types'],'raw_quant_sha256':hashlib.sha256(f.read_bytes()).hexdigest()})
assert len(computed)==12
out=rows(m/'results/DRA014736_target_contrasts.tsv')
for r in out:
 key=r['measure'];cond=r['contrast'].split('_minus_')[0];a=[computed[(d,cond)][key] for d in ['GC40','GC41','GC42']];b=[computed[(d,'control')][key] for d in ['GC40','GC41','GC42']];delta=[x-y for x,y in zip(a,b)]
 assert abs(statistics.mean(delta)-float(r['mean_difference']))<1e-10
 assert sum(x<0 for x in delta)==int(r['negative'])
 if any(delta):
  assert all(x!=0 for x in delta),'Zero differences require explicit tie-aware recheck'
  wp=wilcoxon(delta,method='exact').pvalue;tp=ttest_rel(a,b).pvalue
  assert abs(wp-float(r['P']))<1e-12 and abs(tp-float(r['paired_t_P']))<1e-10
for key in ['log2TPM1','log2CPM1']:
 family=sorted([r for r in out if r['measure']==key and r['contrast']!='combined_minus_control'],key=lambda r:float(r['P']))
 expected=[min(1.,2*float(family[0]['P']),float(family[1]['P'])),float(family[1]['P'])]
 assert all(abs(float(r['FDR_secondary'])-q)<1e-12 for r,q in zip(family,expected))
for r in rows(m/'results/DRA014736_context_gene_samples.tsv'):
 assert abs(float(r['TPM'])-raw_genes[(r['donor'],r['condition'])][r['gene']])<1e-6
for r in rows(m/'results/DRA014736_defined_CDS_summary.tsv'):
 v=computed[(r['donor'],r['condition'])]
 assert abs(float(r['all_gene_TPM'])-v['TPM'])<1e-9 and abs(float(r['defined_CDS_TPM'])-v['defined_CDS_TPM'])<1e-9
# Verify the full primary result preserves the previous independently checked checkpoint.
old=rows(m/'results/DRA_primary_checkpoint/primary_contrasts.tsv')
for r in old:
 if r['measure'] in ['log2TPM1','log2CPM1']:
  new=next(x for x in out if x['measure']==r['measure'] and x['contrast']=='combined_minus_control')
  assert abs(float(new['mean_difference'])-float(r['mean_difference']))<1e-12 and float(new['P'])==float(r['exact_signedrank_P'])
qa={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'independent_donors':3,'libraries':12,'method':'Raw transcript CSV reaggregation; scipy Wilcoxon paired comparison; independent two-test BH calculation; context and CDS reaggregation; unchanged primary checkpoint','checks':checks}
(m/'results/DRA_FULL_INDEPENDENT_QA_200.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
spec=json.loads((m/'DRA_IDENTITY_COMPARATOR_SPECIFICATION_200.json').read_text());diag=[]
for donor in ['GC40','GC41','GC42']:
 for cond in ['combined','TGF','IL1']:
  delta={g:math.log2(1+raw_genes[(donor,cond)][g])-math.log2(1+raw_genes[(donor,'control')][g]) for g in needed}
  for panel,genes in spec['panels'].items():
   mean=statistics.mean(delta[g] for g in genes)
   diag.append({'donor':donor,'condition':cond,'panel':panel,'CGNL1_delta_log2TPM1':delta['CGNL1'],'panel_mean_delta_log2TPM1':mean,'target_minus_panel_delta':delta['CGNL1']-mean})
write('DRA_IDENTITY_COMPARATOR_200.tsv',diag)
print(json.dumps(qa,indent=2));print('ALL THREE CONTRASTS',json.dumps(out));print('IDENTITY DIAGNOSTIC',json.dumps(diag))
