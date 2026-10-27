from pathlib import Path
import importlib.util,json,datetime,shutil
O=Path.cwd()/'mechanism_investigation_v2';out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python_packages':{n:bool(importlib.util.find_spec(n)) for n in ['scanpy','anndata','infercnvpy','scipy','h5py']},'R_packages':[],'local_gene_annotations':[]}
for env in Path('/mnt/BioDisk/wangmingyuan/miniconda3/envs').iterdir():
 lib=env/'lib/R/library'
 if lib.exists():
  found=[n for n in ['copykat','infercnv','Seurat','SingleCellExperiment','biomaRt'] if (lib/n).exists()]
  if found:out['R_packages'].append({'env':env.name,'Rscript':str(env/'bin/Rscript'),'packages':found})
for base in [O/'reference',O.parent.parent/'reference',O.parent.parent/'references',Path('/mnt/BioDisk/wangmingyuan/references'),Path('/mnt/BioDisk/wangmingyuan/reference')]:
 if base.is_dir():
  for p in base.rglob('*'):
   if p.is_file() and any(k in p.name.lower() for k in ['gtf','gene_order','geneorder','gencode']) and p.stat().st_size<300000000:out['local_gene_annotations'].append({'file':str(p),'bytes':p.stat().st_size})
(O/'metadata/CNV_CAPABILITY_161.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out))
