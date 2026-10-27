from pathlib import Path
import json,html,hashlib,zipfile,datetime,shutil
import pandas as pd
R=Path(__file__).resolve().parent
S=json.loads((R/'PROGRESS_SNAPSHOT.json').read_text(encoding='utf-8'))
Q=json.loads((R/'METADATA_QA.json').read_text(encoding='utf-8'))
C=pd.read_csv(R/'patient_pair_eligibility.tsv',sep='\t')
D=json.loads((R/'REQUIRED_DOWNLOADS.json').read_text(encoding='utf-8'))
(R/'raw'/'expression_epi_counts').mkdir(parents=True,exist_ok=True)
shutil.copy2(Path(r'C:\Users\31937\CGNL1_new_cohorts_20261013_v1\delivery\metadata\atlas_loaders\Tian_2023_Nat_Med.py'),R/'sources'/'prior_atlas_loader_Tian_2023_Nat_Med.py')
summary=pd.DataFrame(S['summary'])
files=''.join('<tr><td>'+html.escape(x['bucket_location'])+'</td><td>'+str(x['upload_file_size'])+'</td><td><a href="'+html.escape(x['download_url'],quote=True)+'">Sign-in download</a></td></tr>' for x in D['files'])
report='''<!doctype html><html lang="en"><meta charset="utf-8"><title>CGNL1: SCP2079 qualification</title><style>body{font:16px/1.55 Arial,sans-serif;max-width:1180px;margin:36px auto;padding:0 24px;color:#18313d}table{border-collapse:collapse;font-size:14px;display:block;overflow:auto}td,th{border:1px solid #ccd6db;padding:7px;text-align:left}th{background:#eef3f6}h1,h2{color:#125675}a{color:#126b91}.status{padding:16px;background:#fff3d9;border-left:5px solid #bc8a21}code{overflow-wrap:anywhere}</style><h1>CGNL1: SCP2079 qualification</h1>'''
report+='<p>Internal progress record. '+html.escape(S['updated_utc'])+'. This is a data qualification stage, not a new CGNL1 efficacy result or a client deliverable.</p>'
report+='<div class="status"><strong>Metadata qualification complete; expression matrix qualification pending authenticated download.</strong> No target expression, association, predictive performance, or inferential test has been computed in this round.</div>'
report+='<h2>Fourteen pairs pass the existing metadata QC and coverage gate</h2><p>The public epithelial annotations contain 95,421 unique cells from 27 donor IDs, 46 donor/time specimens and 89 technical library labels. Nineteen donors have at least one cell at both Pre and On. The existing fixed QC retains 88,660 cells; requiring at least 20 retained cells per time point leaves 14 paired donors. Requiring 50 leaves 13. These are metadata-qualified counts, pending full barcode/matrix reconciliation.</p>'
report+='<p>Technical libraries are pooled within donor and time point. Cells, specimen labels, and technical library labels do not replace the biological patient count. Donor/patient IDs agree exactly, biosample/specimen IDs agree, and all 22 exported annotation columns have identical cell sets.</p>'
report+=summary.to_html(index=False,border=0)
report+='<h2>Eight qualifying pairs are explicitly labeled Liver</h2><p>The original tissueSite field contains Liver, Lymph Node, Muscle, and Other. All observed time points of each donor retain the same site label. Only Liver qualifies for the liver subset; ambiguous labels remain unchanged. The old atlas loader explicitly replaced sample_tissue with liver despite its own uncertainty. Its all-liver assignment is not inherited. The liver subset remains descriptive under the preceding validation protocol for fewer than 10 independent pairs.</p>'
report+='<h2>The current portal lists two matrices; the raw-count claim is not yet verified</h2><p>The portal currently lists separate expression_epi_counts and expression_epi_norm directories. The old loader read epi_matrix.mtx.gz and described that older input as normalized. That observation cannot establish the scale of the current files. In the public annotations, total_counts is integer-valued, whereas every nCount_RNA value is fractional. nCount_RNA is therefore not used as raw UMI depth. This metadata observation does not itself establish either matrix scale.</p>'
report+='<p>All fixed panel genes including CGNL1 appear in the public visualization gene catalog. Matrix integrity, actual gene identities, numerical scale, barcode alignment and denominator validity remain to be checked from the actual files.</p>'
report+='<h2>The official download requires an authenticated user</h2><p>The official metadata download returned HTTP 401. Anonymous requests to the publicly listed cloud object locations returned HTTP 403. No credentials were extracted and no access controls were bypassed. Public visualization annotations were obtained through the documented annotation endpoint.</p>'
report+='<p>Sign in to the <a href="https://singlecell.broadinstitute.org/single_cell/study/SCP2079/combined-pd-1-braf-and-mek-inhibition-in-braf-v600e-colorectal-cancer-a-phase-2-trial">official SCP2079 study</a> and obtain the three counts files below. Preserve their directory names so they are not confused with the normalized files.</p><table><tr><th>File</th><th>Expected bytes</th><th>Official link</th></tr>'+files+'</table>'
report+='<p>Local destination: <code>C:\\Users\\31937\\CGNL1_SCP2079_20261014_v1\\raw\\expression_epi_counts\\</code>. The normalized matrix is a conditional reference if the counts file fails scale qualification.</p>'
report+='<h2>The analysis plan is prepared, with measurement frozen only after qualification</h2><p>Retain CGNL1 and the existing actomyosin and cycling panels. The intended main analysis uses within-patient On-minus-Pre changes, with depth, epithelial identity and mature-colonocyte changes as controls. Full rank, adequate residual degrees of freedom, and at least 10 paired donors are required under this protocol. Target-effect raw P&lt;0.05 is the nominal screening rule; BH FDR is reported for the two endpoint hypotheses within each predefined model/population family. Fourteen pairs do not guarantee model estimability. Liver-only results remain descriptive.</p>'
report+='<p>Patient-held-out prediction is conditional on all training folds being estimable. No original cross-cohort coefficients will be transferred without verified compatible units. If only normalized expression can be verified, a separate author-scale aggregation specification must be frozen before target outcomes are viewed. No invented raw counts will be reconstructed.</p>'
report+='<p>The study is a single-arm BRAF/MEK/PD-1 treatment cohort. Paired changes can support treatment-context associations; they do not identify a CGNL1 perturbation effect or establish a general liver-colonization mechanism. PFS/RECIST remain descriptive metadata in this round.</p>'
report+='<h2>Every patient remains visible</h2>'+C.to_html(index=False,border=0)
report+='<h2>Traceability</h2><ul>'
for p in ['patient_pair_eligibility.tsv','patient_metadata.tsv','technical_library_coverage.tsv','qualification_summary.tsv','METADATA_QA.json','CONDITIONAL_ANALYSIS_PLAN.json','REQUIRED_DOWNLOADS.json','PROGRESS_SNAPSHOT.json','METADATA_FETCH.json','GENE_AVAILABILITY.json','PUBLIC_OBJECT_PROBE.json']:
 report+='<li><a href="'+p+'">'+p+'</a></li>'
report+='</ul><p>Primary sources: <a href="https://www.nature.com/articles/s41591-022-02181-8">Tian et al., Nature Medicine (2023)</a>; <a href="https://singlecell.broadinstitute.org/single_cell/api/v1/site/studies/SCP2079">current official file registry</a>; original per-cell annotation exports retained in metadata/. The paper reports 23 whole-biopsy pairs and 19 pairs in the epithelial signature panels; the operational sample count here is recalculated from actual public annotations and the existing project gates.</p></html>'
(R/'SCP2079_qualification_report.html').write_text(report,encoding='utf-8')
# File-level integrity manifest, then a compressed transfer bundle.
manifest=[]
for p in sorted(R.rglob('*')):
 if p.is_file() and p.suffix!='.zip' and p.name not in ['FILE_MANIFEST.json','BUNDLE_RECEIPT.json','SERVER_RECEIPT.json','DELIVERY_QA.json']:
  manifest.append({'path':p.relative_to(R).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(R/'FILE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
qa={'all_manifest_files_reopened_and_hash_verified':all(hashlib.sha256((R/x['path']).read_bytes()).hexdigest()==x['sha256'] for x in manifest),'json_readback':True,'html_links_exist':True,'metadata_rows':95421,'files':len(manifest)}
for p in R.glob('*.json'):json.loads(p.read_text(encoding='utf-8'))
import re
for link in re.findall('href="([^"]+)"',report):
 if not link.startswith('https://'):assert (R/link).is_file(),link
(R/'DELIVERY_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
zp=R/'SCP2079_qualification_bundle.zip'
with zipfile.ZipFile(zp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for item in manifest:z.write(R/item['path'],item['path'])
 z.write(R/'FILE_MANIFEST.json','FILE_MANIFEST.json');z.write(R/'DELIVERY_QA.json','DELIVERY_QA.json')
with zipfile.ZipFile(zp) as z:assert z.testzip() is None
receipt={'path':zp.name,'bytes':zp.stat().st_size,'sha256':hashlib.sha256(zp.read_bytes()).hexdigest(),'files':len(manifest)+2}
(R/'BUNDLE_RECEIPT.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt),flush=True)
print(C[(C.definition=='fixed_qc')&C.paired_ge20].to_string(index=False),flush=True)