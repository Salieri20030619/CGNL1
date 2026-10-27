from pathlib import Path
import sys,json,hashlib,csv,datetime
root=Path('/mnt/BioData/bio/wangmingyuan/CGNL1_vascular_context_20261018')
sys.path.insert(0,str(root/'software/blake3_516'))
import blake3
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
D=root/'technical_probe_513';p=D/'3M-february-2018.malva.txt'
h=blake3.blake3()
with p.open('rb') as f:
    for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
registry=json.loads((D/'simpleaf_chemistries.json').read_text(encoding='utf-8'))['10xv3']
expected=registry['plist_name'];observed=h.hexdigest()
spec=json.loads((B/'metadata/PILOT_QUANTIFICATION_SPEC_496.json').read_text(encoding='utf-8'))
mapping={}
with (Path(spec['reference']).parent/'transcript_to_gene.tsv').open(encoding='utf-8') as f:
    for row in csv.DictReader(f,delimiter='\t'): mapping[row['transcript_id']]=row['gene_symbol']
target_rows=[];mt=0.;total=0.
with (D/'auto_strand_quant/quant.sf').open(encoding='utf-8') as f:
    for row in csv.DictReader(f,delimiter='\t'):
        symbol=mapping.get(row['Name']);reads=float(row['NumReads']);total+=reads
        if symbol and symbol.startswith('MT-'):mt+=reads
        if symbol=='CGNL1':target_rows.append(row)
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),reference_path=str(p),hash_algorithm='BLAKE3',hash_library_version=blake3.__version__,expected_registry_fingerprint=expected,observed_fingerprint=observed,reference_fingerprint_matches=expected==observed,algorithm_source='https://github.com/COMBINE-lab/simpleaf/blob/main/src/utils/prog_utils.rs',prior_SHA256_comparison_not_valid_equality_test=True,technical_R2_probe=dict(sample='first 1000000 R2 reads',CGNL1_reference_transcripts=len(target_rows),CGNL1_estimated_assigned_reads=sum(float(r['NumReads']) for r in target_rows),total_estimated_assigned_reads=total,mitochondrial_estimated_assigned_reads=mt,mitochondrial_assigned_fraction=mt/total if total else None,biological_inference=False))
(B/'metadata/REFERENCE_FINGERPRINT_516.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))

