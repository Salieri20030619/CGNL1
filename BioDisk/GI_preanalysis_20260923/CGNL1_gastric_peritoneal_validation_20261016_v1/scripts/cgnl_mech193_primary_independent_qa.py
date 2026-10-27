
from pathlib import Path
import csv,gzip,json,math,statistics,datetime,hashlib
from scipy.stats import wilcoxon,ttest_rel,t
p=Path.cwd();m=p/'mechanism_investigation_v2';d=m/'results/DRA_primary_checkpoint'
def rows(f):
 with (gzip.open(f,'rt',encoding='utf-8') if str(f).endswith('.gz') else open(f,encoding='utf-8')) as h:return list(csv.DictReader(h,delimiter='	'))
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
ann=rows(m/'reference/gencode_v47/transcript_to_gene.tsv')
target={r['transcript_id'] for r in ann if r['gene_id'].split('.')[0]=='ENSG00000128849'}
assert len(target)==3
coding=json.loads((m/'metadata/TARGET_TRANSCRIPT_BIOTYPES_93.json').read_text())['CGNL1_coding_transcripts']
original=rows(d/'target_samples.tsv');checks=[];computed={}
for r in original:
 f=m/'processed/DRA014736_quant'/r['run']/'quant.sf'
 tr=rows(f);selected={a['Name']:a for a in tr if a['Name'] in target};assert set(selected)==target
 tpm=math.fsum(float(a['TPM']) for a in selected.values());nr=math.fsum(float(a['NumReads']) for a in selected.values());total=math.fsum(float(a['NumReads']) for a in tr);ct=math.fsum(float(selected[x]['TPM']) for x in coding)
 calc={'TPM':tpm,'estimated_NumReads':nr,'log2TPM1':math.log2(1+tpm),'log2CPM1':math.log2(1+nr/total*1e6),'defined_CDS_TPM':ct,'log2definedCDSTPM1':math.log2(1+ct)}
 error=max(abs(v-float(r[k])) for k,v in calc.items());assert error<1e-9
 computed[(r['donor'],r['condition'])]=calc
 checks.append({'run':r['run'],'max_error':error,'raw_quant_sha256':sha(f)})
for r in rows(d/'primary_contrasts.tsv'):
 key=r['measure'];a=[computed[(donor,'combined')][key] for donor in ['GC40','GC41','GC42']];b=[computed[(donor,'control')][key] for donor in ['GC40','GC41','GC42']];delta=[x-y for x,y in zip(a,b)]
 assert abs(statistics.mean(delta)-float(r['mean_difference']))<1e-10
 if key!='log2definedCDSTPM1':
  w=wilcoxon(a,b,method='exact',alternative='two-sided');tt=ttest_rel(a,b)
  assert abs(w.pvalue-float(r['exact_signedrank_P']))<1e-12 and abs(tt.pvalue-float(r['paired_t_P']))<1e-10
  ci=t.ppf(.975,2)*statistics.stdev(delta)/math.sqrt(3)
  assert abs(statistics.mean(delta)-ci-float(r['mean_CI_low']))<1e-10
 checks.append({'measure':key,'mean_delta':statistics.mean(delta),'negative':sum(x<0 for x in delta)})
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'pass':True,'method':'Independent csv parsing of raw transcript quantification, math.fsum, scipy exact Wilcoxon and paired t; compared all target sample measures and primary summaries','checks':checks,'original_files_unchanged':True,'scope':'completed primary endpoint only; full twelve-library analysis QA pending'}
dest=m/'results/DRA_PRIMARY_INDEPENDENT_QA_193.json';assert not dest.exists();dest.write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))
