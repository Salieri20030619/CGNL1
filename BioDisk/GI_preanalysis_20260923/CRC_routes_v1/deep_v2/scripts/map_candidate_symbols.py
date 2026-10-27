import pathlib,json,pandas as pd,gzip,re,sys
O=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2');a=pd.read_csv(O/'raw/bulk_FPKM.tsv.gz',sep='\t',index_col=0);info=pd.read_csv(O/'metadata/Homo_sapiens.gene_info.gz',sep='\t',low_memory=False);rows=[];symbols=set(info.Symbol.astype(str));alias={}
for _,r in info.iterrows():
 for syn in str(r.Synonyms).split('|'):
  if syn!='-':alias.setdefault(syn,set()).add(str(r.Symbol))
for gene in json.loads((O/'PLAN.json').read_text())['frozen_genes']:
 g=gene.replace('.','-');matches=[]
 if g in a.index:matches=[g];rule='exact or R-name reversal'
 else:matches=[old for old in a.index if alias.get(old)=={g}];rule='unique NCBI synonym'
 rows.append({'candidate':gene,'canonical':g,'bulk_gene':matches[0] if len(matches)==1 else '', 'mapping_rule':rule if len(matches)==1 else 'absent_or_ambiguous','matches':';'.join(matches)})
x=pd.DataFrame(rows);x.to_csv(O/'metadata/candidate_symbol_map.tsv',sep='\t',index=False);print(x[x.mapping_rule!='exact or R-name reversal'].to_string(index=False))
with gzip.open(O/'metadata/GRCh38.p14.ncbi.gtf.gz','rt') as f:
 k=0
 for line in f:
  if line.startswith('#'):continue
  cols=line.rstrip().split('\t')
  if cols[2]=='gene':print('GTF_GENE',line[:300]);k+=1
  if k==2:break
sys.path.insert(0,str(O/'vendor'))
try:import infercnvpy;print('INFER_IMPORT_OK',infercnvpy.__version__)
except Exception as e:print('INFER_IMPORT_ERROR',repr(e))
