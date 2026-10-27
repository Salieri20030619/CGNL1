import pathlib,re,gzip,json,pandas as pd,numpy as np
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2');R=O/'results/second_validation';R.mkdir(exist_ok=True);s=pd.read_csv(O/'metadata/GSE14297_samples.tsv',sep='\t');rows=[]
for _,r in s.iterrows():
 site='LM' if '-LM-' in r.title else ('PT' if '-Tu-' in r.title else ('NL' if '-NL-' in r.title else 'NC'));match=re.search(r'Patient_No\.(\d+)',r.title);pt=match.group(1) if match else r.title.split('_')[-1]+'_'+site;origin=r.characteristics.split(';')[0].replace('Liver Metastasis of a Colorectal Cancer from ','');strict=site in ['PT','LM'] and ('colon' in origin.lower() or 'coli' in origin.lower() or origin in ['Sigma','Coecum']);rows.append({'sample':r['sample'],'patient':pt,'site':site,'origin':origin,'strict_colon':strict})
m=pd.DataFrame(rows);assert m[m.site=='PT'].shape[0]==18;assert m[(m.site=='PT')&m.strict_colon].shape[0]==13;m.to_csv(R/'design.tsv',sep='\t',index=False)
expr=pd.read_csv(O/'metadata/GSE14297.matrix.gz',sep='\t',comment='!',index_col=0,encoding='utf-8');expr=expr.dropna(how='all');expr.index=expr.index.astype(str);print('MATRIX',expr.shape,'QUANTILES',np.quantile(expr.values,[0,.1,.5,.9,1]));ann=pd.read_csv(O/'metadata/GPL6370_annotation.tsv.gz',sep='\t',dtype=str);mapping={}
with gzip.open(O/'metadata/GRCh38.p14.ncbi.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  c=line.rstrip().split('\t')
  if c[2] not in ['transcript','CDS']:continue
  gene=re.search(r'(?:^|; )gene "([^"]+)"',c[8]);gene=gene.group(1) if gene else re.search(r'gene_id "([^"]+)"',c[8]).group(1)
  for match in re.findall(r'(?:transcript_id|protein_id) "([NX][MPR]_\d+)(?:\.\d+)?"',c[8]):mapping.setdefault(match,set()).add(gene)
info=pd.read_csv(O/'metadata/Homo_sapiens.gene_info.gz',sep='\t',low_memory=False);current=set(info.Symbol.astype(str));aliases={}
for _,v in info.iterrows():
 for alias in str(v.Synonyms).split('|'):
  if alias!='-':aliases.setdefault(alias,set()).add(str(v.Symbol))
rows=[]
for _,r in ann.iterrows():
 if r.ID not in expr.index:continue
 label=str(r.SPOT_ID);accessions=re.findall(r'[NX][MPR]_\d+',label);gs=({label} if label in current else aliases.get(label,set()).copy())
 if re.match(r'^\d+\*',label) or label.startswith('i'):accessions=[];gs=set()
 for ac in accessions:gs.update(mapping.get(ac,set()))
 rows.append({'probe_id':r.ID,'SPOT_ID':r.SPOT_ID,'gene':next(iter(gs)) if len(gs)==1 else '', 'mapping_status':'unique' if len(gs)==1 else ('unmapped' if not gs else 'ambiguous')})
a=pd.DataFrame(rows);a.to_csv(R/'probe_gene_map.tsv',sep='\t',index=False);valid=a[a.mapping_status=='unique'];g=expr.loc[valid.probe_id].copy();g.index=valid.gene.values;g=g.groupby(level=0).mean().dropna();g.to_csv(R/'gene_expression.tsv.gz',sep='\t',index_label='gene');print('GENE_MATRIX',g.shape);print('CGNL1_PROBES',valid[valid.gene=='CGNL1'].to_string(index=False));print(m.groupby(['site','strict_colon']).size().to_string());print('CGNL1_MEANS',[(site,float(g.loc['CGNL1',m[m.site==site]['sample']].mean())) for site in ['PT','LM','NL','NC']]);assert 'CGNL1' in g.index
(O/'results/second_validation/PREPARE_COMPLETE.json').write_text(json.dumps({'exit_code':0,'all_pairs':18,'strict_colon_pairs':13,'genes':len(g),'probe_mapping':'exact current symbol, unique NCBI synonym, or unique RefSeq accession; multiplied or i-prefixed mapping labels conservatively excluded; all-NA filtered probes removed; multi-probe gene mean without outcome selection','scale':'author variance-stabilized and spline-normalized signals; no second log transform'},indent=2))
