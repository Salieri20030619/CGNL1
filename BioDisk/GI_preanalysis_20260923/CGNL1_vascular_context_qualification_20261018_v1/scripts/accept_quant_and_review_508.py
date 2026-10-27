from pathlib import Path
import datetime,gzip,hashlib,json,os,shutil,subprocess,sys,traceback
import numpy as np
from scipy import io,sparse
B=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1')
def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(name,data):
    p=B/'metadata'/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8');t.replace(p)
with (B/'metadata/REVIEW508_LAUNCH.json').open('x',encoding='utf-8') as f:
    json.dump(dict(pid=os.getpid(),wrapper_pid=os.getpid(),started=now(),receipt=str(B/'metadata/REVIEW508.exit.json'),log=str(B/'logs/REVIEW508.log')),f,indent=2)
try:
    spec=json.loads((B/'metadata/PILOT_QUANTIFICATION_SPEC_496.json').read_text(encoding='utf-8'))
    original=Path(spec['output'])
    previous=json.loads((B/'metadata/PILOT496.exit.json').read_text(encoding='utf-8'))
    assert previous['exit_code']==1 and 'AssertionError: alevin/quants_mat.mtx' in previous['error']
    meta=json.loads((original/'aux_info/meta_info.json').read_text(encoding='utf-8'))
    qa=json.loads((B/'metadata/FASTQ_QA_COMPLETE_496.json').read_text(encoding='utf-8'))
    assert not meta['quant_errors'] and meta['num_processed']==qa['read_pairs'] and meta['end_time']
    assert 'Finished optimizer' in (B/'logs/PILOT496.log').read_text(encoding='utf-8')
    view=original.with_name('pilot_quantification_view_508')
    view.mkdir(exist_ok=False);(view/'alevin').mkdir()
    source=original/'alevin/quants_mat.mtx.gz';matrix=view/'alevin/quants_mat.mtx'
    with gzip.open(source,'rb') as r,matrix.open('xb') as w: shutil.copyfileobj(r,w)
    for p in (original/'alevin').iterdir():
        if p.name!='quants_mat.mtx': (view/'alevin'/p.name).symlink_to(p)
    for p in original.iterdir():
        if p.name!='alevin': (view/p.name).symlink_to(p,target_is_directory=p.is_dir())
    rows=(view/'alevin/quants_mat_rows.txt').read_text(encoding='utf-8').splitlines()
    cols=(view/'alevin/quants_mat_cols.txt').read_text(encoding='utf-8').splitlines()
    m=sparse.csr_matrix(io.mmread(matrix));m.eliminate_zeros()
    assert m.shape==(len(rows),len(cols)) and len(set(rows))==len(rows) and len(set(cols))==len(cols)
    assert m.nnz>0 and np.isfinite(m.data).all() and (m.data>=0).all()
    acceptance=dict(utc=now(),accepted=True,source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),matrix_sha256=hashlib.sha256(matrix.read_bytes()).hexdigest(),shape=list(m.shape),nnz=int(m.nnz),estimated_UMI_sum=float(m.sum()),compressed_gzip_integrity_passed=True,original_receipt_preserved=True,quantification_rerun=False,biological_validation=False,quantifier_exit_evidence='Original wrapper passed returncode==0 assertion then failed expected uncompressed filename check; original failed wrapper receipt retained')
    save('QUANT_OUTPUT_ACCEPTANCE_508.json',acceptance)
    spec.update(output=str(view),source_quantification_output=str(original),format_correction='Lossless gzip decompression into separate review view')
    save('PILOT_QUANTIFICATION_SPEC_508.json',spec)
    save('PILOT508.exit.json',dict(exit_code=0,operation='accept_existing_quantification_output',finished=now()))
    save('PILOT_QUANT_COMPLETE_508.json',dict(required_files_present=True,acceptance_record=str(B/'metadata/QUANT_OUTPUT_ACCEPTANCE_508.json'),biological_validation=False))
    del m
    cell=(B/'scripts/qualify_pilot_cells_486.py').read_text(encoding='utf-8').replace('SPEC_485','SPEC_508').replace('PILOT485.exit','PILOT508.exit').replace('COMPLETE_485','COMPLETE_508').replace('pilot_cell_qualification_457','pilot_cell_qualification_509').replace('COMPLETE_457','COMPLETE_509')
    cp=B/'scripts/qualify_pilot_cells_509.py';compile(cell,str(cp),'exec')
    with cp.open('x',encoding='utf-8') as f:f.write(cell)
    save('REVIEW508_STATE.json',dict(stage='cell_qualification',started=now()))
    r=subprocess.run([sys.executable,str(cp)])
    save('CELL509.exit.json',dict(exit_code=r.returncode,finished=now()))
    assert r.returncode==0,'Cell qualification failed'
    audit=(B/'scripts/audit_pilot_500.py').read_text(encoding='utf-8').replace('CELL497.exit','CELL509.exit').replace('SPEC_496','SPEC_508').replace('qualification_497','qualification_509').replace('COMPLETE_497','COMPLETE_509').replace('audit_500','audit_510').replace('AUDIT500','AUDIT510').replace('alevin/alevin_meta_info.json','aux_info/alevin_meta_info.json')
    ap=B/'scripts/audit_pilot_510.py';compile(audit,str(ap),'exec')
    with ap.open('x',encoding='utf-8') as f:f.write(audit)
    save('REVIEW508_STATE.json',dict(stage='numeric_audit',started=now()))
    r=subprocess.run([sys.executable,str(ap)]);assert r.returncode==0,'Numeric audit failed'
    status=dict(stage='ready_for_scientific_review',exit_code=0,finished=now(),mainline_validated=False)
except Exception:
    status=dict(stage='failed_needs_review',exit_code=1,finished=now(),error=traceback.format_exc())
save('REVIEW508_STATE.json',status);save('REVIEW508.exit.json',status)
print(json.dumps(status),flush=True)

