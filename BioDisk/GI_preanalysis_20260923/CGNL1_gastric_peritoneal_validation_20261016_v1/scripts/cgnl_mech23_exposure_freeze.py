from pathlib import Path
import json,datetime,hashlib,gzip
O=Path.cwd()/'mechanism_investigation_v2';p=O/'EXPOSURE_ANALYSIS_SPECIFICATION.json'
s={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'primary_dataset':'GSE223165','primary_target':'CGNL1 / ENSG00000128849','purpose':'test whether tumor-derived exposure changes expression within primary human mesothelial cells; ovarian exposure is not gastric validation','comparisons':[{'name':'CM_minus_control','independent_unit':'mesothelial donor','expected_n':10,'primary_measure':'paired log2(authorFPKM+1) difference','test':'two-sided exact Wilcoxon signed-rank if no zero/tied absolute differences; otherwise permutation signed-rank','secondary':'paired mean effect and t interval/P; raw-count full-library log2CPM with+1 sensitivity'}, {'name':'mean_ASC_minus_control','independent_unit':'mesothelial donor, average over the same three shared ascites exposures','expected_n':7,'status':'exploratory','test':'paired signed-rank only with valid donor repeats; retain three separate ascites directions descriptively','limitation':'three ascites sources do not multiply independent mesothelial donors, and common ascites limits generalization'}],'multiple_testing':'BH across two target primary contrasts; secondary test and normalization sensitivity separate, not select best P','gene_family':'single prespecified target per contrast; any additional genes or programs require distinct recorded family','qualification':['verify complete gzip and numerical matrix integrity','match actual columns to exact donor and treatment identifiers in GEO descriptions; no matching by expression','no inferred missing patient or pairing','target unique Ensembl row; record gene-name mapping','nonnegative finite author values and correct sample counts'],'limitations':['in vitro exposure and ovarian source','FPKM/RPKM naming inconsistency in GEO; preserve author scale label','do not infer tumor promoting vs protective function from induced RNA alone','n<10 comparisons exploratory; no cell-level inference']}
assert not p.exists();p.write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8');print('FROZEN',hashlib.sha256(p.read_bytes()).hexdigest())
for f in sorted((O/'raw').glob('GSE223165*.gz')):
 try:
  with gzip.open(f,'rt',encoding='utf-8') as h:header=h.readline()
  print(f.name,header[:4000])
 except Exception as e:print(f.name,str(e))
# GEO descriptions may contain exact sequencing sample identifiers.
t=(O/'metadata/GSE223165_samples.soft').read_bytes().decode('utf-8')
for b in t.split('^SAMPLE = ')[1:]:
 ls=b.splitlines();print(ls[0],[l for l in ls if l.startswith(('!Sample_title =','!Sample_description =','!Sample_supplementary_file'))])
