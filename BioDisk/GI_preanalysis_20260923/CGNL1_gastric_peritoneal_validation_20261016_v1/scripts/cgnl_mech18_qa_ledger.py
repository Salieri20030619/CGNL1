from pathlib import Path
import csv,gzip,json,math,statistics,datetime,hashlib
O=Path.cwd()/'mechanism_investigation_v2';p=O/'raw/GSE228955_counts_ftp.gz'
with gzip.open(p,'rt',encoding='utf-8') as f:
 reader=csv.reader(f,delimiter='\t');head=next(reader);original=list(reader)
# First header field is blank, and the first data field is gene symbol.
assert len(head)==13 and head[0]==''
with (O/'results/GSE228955_prespecified_machinery.tsv').open(encoding='utf-8',newline='') as f:res=list(csv.DictReader(f,delimiter='\t'))
maxerr=0.;checked=0
for r in res:
 match=[v for v in original if v[0]==r['gene']]
 if len(match)!=1:continue
 v=match[0];means={}
 for group,prefix in [('WT','WT'),('CGNL1_KO','CGNL1_KO'),('CGN_KO','CGN_KO'),('double_KO','CGN_CGNL1')]:
  ids=[i for i,c in enumerate(head) if c.startswith(prefix)];means[group]=statistics.mean(math.log2(float(v[i])+1) for i in ids)
 for g in ['CGNL1_KO','CGN_KO','double_KO']:maxerr=max(maxerr,abs(means[g]-means['WT']-float(r[g+'_minus_WT'])))
 checked+=1
assert maxerr<1e-10
qa={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'independent_csv_and_scalar_math_validation':True,'genes_checked':checked,'maximum_absolute_difference':maxerr,'rows':len(original),'duplicate_gene_symbols':sorted({r[0] for r in original if sum(t[0]==r[0] for t in original)>1}) if False else 'see KO_duplicate_symbols.tsv','all_pass':True}
(O/'results/KO_NUMERICAL_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
ledger=[
 {'id':'MEC-01','candidate':'CGNL1 deficiency directly disrupts mesothelial permeability barrier','support':'CGNL1 participates in junction regulation in epithelial models','counterevidence':'2023 MDCK knockout study reports higher TER, no 3kDa dextran permeability change; shRNA and KO differ','current_verdict':'not established; do not use simple barrier-loss chain','source':'https://pmc.ncbi.nlm.nih.gov/articles/PMC10417749/','next_discriminating_evidence':'CGNL1-specific perturbation and rescue in qualified human mesothelial monolayer, with permeability and cancer crossing measured separately'},
 {'id':'MEC-02','candidate':'CGNL1-dependent junction–cytoskeleton coupling changes host surface mechanics relevant to colonization','support':'NM2 binding and junction tethering have direct experimental support in epithelial models; CAMSAP3 junction recruitment has KO/rescue support','counterevidence':'No measured gastric PM host-cell force/localization or crossing effect; mechanisms already known outside this disease','current_verdict':'mechanistically grounded candidate, not established disease mechanism','sources':['https://pmc.ncbi.nlm.nih.gov/articles/PMC10202830/','https://pmc.ncbi.nlm.nih.gov/articles/PMC10184829/'],'next_discriminating_evidence':'disease-linked source-specific CGNL1 protein redistribution or mechanical phenotype; specificity/rescue and colonization endpoint'},
 {'id':'MEC-03','candidate':'CGNL1 drives an MMT transcriptional program','support':'disease MMT is independently established in GC literature','counterevidence':'prior 15-patient mesothelial program analysis did not establish independent association; MDCK KO transcript changes are mixed and not mesothelium','current_verdict':'current data do not support this as mainline','next_discriminating_evidence':'same-source exposure/disease data with direct target measurement and appropriate donor replication'},
 {'id':'MEC-04','candidate':'CGNL1 bulk association is normal host-cell abundance / sampling composition','support':'strong mesothelial-like enrichment; no proven disease-state change in same lineage','counterevidence':'original bulk conditional target association remains positive; composition proxies incomplete','current_verdict':'live alternative, not ruled out','next_discriminating_evidence':'within-source disease comparison and tissue localization with clinical anchors'}]
(O/'results/MECHANISM_EVIDENCE_LEDGER_v1.json').write_text(json.dumps(ledger,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(qa));print('ledger entries',len(ledger))
