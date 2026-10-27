from pathlib import Path
import json,csv,openpyxl,math
import numpy as np
from scipy.stats import spearmanr
b=Path(__file__).resolve().parents[1];ann={}
with (b/'raw/PDC_original_384/CPTAC3_Pancreatic_Ductal_Adenocarcinoma_Proteome.tmt11.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['NCBIGeneID'].isdigit() and r['Gene']!='CGNL1':ann[r['Gene']]=r['NCBIGeneID']
w=openpyxl.load_workbook(b/'raw/Fudan_quantification_TableS3_363.xlsx',read_only=True,data_only=True);it=w['Table S3A '].iter_rows(values_only=True);h=next(it);cols=[i for i,k in enumerate(h) if isinstance(k,str) and k.startswith('PDAC_')];samples=[h[i] for i in cols];genes=[];matrix=[]
for r in it:
 if r[0] in ann:
  genes.append(r[0]);matrix.append([float(r[i]) if isinstance(r[i],(int,float)) and r[i]>0 else np.nan for i in cols])
x=np.array(matrix);pilot=json.loads((b/'metadata/FUDAN_RAW_PROTEIN_PILOT_400.json').read_text());out=[]
for p in pilot:
 rows=list(csv.DictReader(Path(p['file']).open(),delimiter='\t'));key=next(k for k in rows[0] if k.startswith('Exp'));rmap={r['Accession']:float(r[key]) for r in rows};raw=np.array([rmap.get(ann[g],np.nan) for g in genes]);scores=[]
 for j,k in enumerate(samples):
  ok=np.isfinite(raw)&(raw>0)&np.isfinite(x[:,j]);n=int(sum(ok))
  if n<1000:continue
  lr=np.log2(x[ok,j]/raw[ok]);med=float(np.median(lr));mad=float(np.median(abs(lr-med)));rho=float(spearmanr(raw[ok],x[ok,j]).statistic);exact=float(np.mean(abs(x[ok,j]-np.round(raw[ok],4))<=5.1e-5));scores.append({'published_sample':k,'n_anchor_genes':n,'median_log2_scale':med,'median_absolute_log2_residual':mad,'spearman':rho,'fraction_equal_at4decimals_without_rescaling':exact})
 named=p['sample'].replace('NAT','N');out.append({'raw_sample_filename_label':p['sample'],'raw_file':p['file'],'target_excluded_from_fingerprint':True,'named_sample_score':next((r for r in scores if r['published_sample']==named),None),'best5_by_scale_residual':sorted(scores,key=lambda r:r['median_absolute_log2_residual'])[:5]})
res={'annotation_source':'PDC original protein table Gene/NCBIGeneID single-valued fields; CGNL1 excluded','common_annotated_proteins':len(genes),'fingerprints':out,'no_relabeling_performed':True,'interpretation':'Fingerprint diagnostic only; does not by itself establish donor identity or authorize replacing source labels'};(b/'metadata/FUDAN_NONTARGET_FINGERPRINT_403.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print(json.dumps(res))
