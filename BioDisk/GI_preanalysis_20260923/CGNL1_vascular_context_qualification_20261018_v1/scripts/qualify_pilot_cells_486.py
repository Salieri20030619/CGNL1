from pathlib import Path
import json,csv,gzip,hashlib,datetime
import numpy as np
from scipy import sparse,io
b=Path(__file__).resolve().parents[1]
def classify(mat,symbols,spec):
 mat=sparse.csr_matrix(mat,dtype=float);assert np.isfinite(mat.data).all() and (mat.data>=0).all()
 heldout=[i for i,s in enumerate(symbols) if s!='CGNL1'];qc_mat=mat[:,heldout];qc_mat.eliminate_zeros();totals=np.asarray(qc_mat.sum(axis=1)).ravel();det=np.diff(qc_mat.indptr)
 mt=[i for i,s in enumerate(symbols) if s.startswith('MT-')];assert mt,'No mitochondrial features mapped'
 mtcounts=np.asarray(mat[:,mt].sum(axis=1)).ravel();mtfrac=np.divide(mtcounts,totals,out=np.ones_like(totals),where=totals>0)
 norm=sparse.diags(np.divide(10000.,totals,out=np.zeros_like(totals),where=totals>0))@mat;norm.data=np.log1p(norm.data)
 identity=spec['endothelial_identity'];panels={'endothelial':identity['positive_panel'],'epithelial':identity['epithelial_comparators'],'mural':identity['mural_comparators'],**spec['additional_liver_identity_competitors']};scores={};hits={}
 for group,genes in panels.items():
  assert 'CGNL1' not in genes
  indices=[]
  for g in genes:
   found=[i for i,s in enumerate(symbols) if s==g];assert len(found)==1,(g,len(found));indices+=found
  scores[group]=np.asarray(norm[:,indices].mean(axis=1)).ravel();hits[group]=np.asarray((mat[:,indices]>0).sum(axis=1)).ravel()
 q=spec['cell_QC_primary'];qc=(totals>=q['min_total_UMI'])&(det>=q['min_detected_genes'])&(mtfrac<=q['max_mitochondrial_UMI_fraction'])
 candidates=qc&(hits['endothelial']>=2)
 ambiguous=np.zeros(len(totals),dtype=bool)
 for group in panels:
  if group=='endothelial':continue
  candidates &= scores['endothelial']>scores[group];ambiguous|=hits[group]>=2
 return {'totals':totals,'detected':det,'mitochondrial_fraction':mtfrac,'qc':qc,'candidate':candidates,'ambiguous':ambiguous,'qualified':candidates&~ambiguous,'scores':scores,'hits':hits}
def main():
 spec=json.loads((b/'metadata/PILOT_CELL_QUALIFICATION_AMENDMENT_458.json').read_text());qspec=json.loads((b/'metadata/PILOT_QUANTIFICATION_SPEC_485.json').read_text());out=Path(qspec['output']);receipt=json.loads((b/'metadata/PILOT485.exit.json').read_text());assert receipt['exit_code']==0
 assert json.loads((b/'metadata/PILOT_QUANT_COMPLETE_485.json').read_text())['required_files_present']
 cells=(out/'alevin/quants_mat_rows.txt').read_text().splitlines();genes=(out/'alevin/quants_mat_cols.txt').read_text().splitlines();assert len(cells)==len(set(cells)) and len(genes)==len(set(genes))
 mat=sparse.csr_matrix(io.mmread(out/'alevin/quants_mat.mtx'));mat.eliminate_zeros();assert mat.shape==(len(cells),len(genes))
 mapping=Path(qspec['reference']).parent/'transcript_to_gene.tsv';lookup={}
 with mapping.open(encoding='utf-8') as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['gene_id'] in lookup:assert lookup[r['gene_id']]==r['gene_symbol']
   lookup[r['gene_id']]=r['gene_symbol']
 missing=[g for g in genes if g not in lookup];assert not missing,missing[:10];symbols=[lookup[g] for g in genes];res=classify(mat,symbols,spec)
 ti=[i for i,s in enumerate(symbols) if s=='CGNL1'];assert len(ti)==1;target=np.asarray(mat[:,ti[0]].toarray()).ravel()
 rd=b/'results/pilot_cell_qualification_457';rd.mkdir(exist_ok=True)
 with gzip.open(rd/'cell_qualification.tsv.gz','wt',encoding='utf-8',newline='') as f:
  fields=['barcode','QC_estimated_UMI_excluding_CGNL1','QC_detected_genes_excluding_CGNL1','mitochondrial_fraction','QC_pass','endothelial_candidate','competitor_coexpression_flag','ambiguity_filtered_endothelial','CGNL1_estimated_UMI']+[g+'_score' for g in res['scores']];w=csv.writer(f,delimiter='\t');w.writerow(fields)
  for i,c in enumerate(cells):w.writerow([c,res['totals'][i],res['detected'][i],res['mitochondrial_fraction'][i],bool(res['qc'][i]),bool(res['candidate'][i]),bool(res['ambiguous'][i]),bool(res['qualified'][i]),target[i]]+[v[i] for v in res['scores'].values()])
 groups={}
 for name,mask in [('all_called',np.ones(len(cells),dtype=bool)),('QC_pass',res['qc']),('endothelial_candidate',res['candidate']),('ambiguity_filtered_endothelial',res['qualified'])]:
  n=int(mask.sum());groups[name]={'cells':n,'CGNL1_estimated_UMI_sum':float(target[mask].sum()),'CGNL1_positive_cells':int((target[mask]>0).sum()),'CGNL1_detection_fraction':float((target[mask]>0).mean()) if n else None}
 result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'independent_donors':1,'groups':groups,'candidate_cell_count_gate_passed':groups['ambiguity_filtered_endothelial']['cells']>=50,'manual_marker_mapping_review_required':True,'automatic_full_cohort_launch_authorized_by_this_result':False,'biological_P':None,'count_semantics':'Alevin estimated gene-cell UMI abundance, not guaranteed integer raw counts','source_matrix_sha256':hashlib.sha256((out/'alevin/quants_mat.mtx').read_bytes()).hexdigest()}
 (rd/'CELL_QUALIFICATION_COMPLETE_457.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))
if __name__=='__main__':main()
