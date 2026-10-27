from pathlib import Path
import json,openpyxl,math,datetime
b=Path(__file__).resolve().parents[1];p=b/'raw/41467_2021_25202_MOESM3_ESM.xlsx';w=openpyxl.load_workbook(p,read_only=True,data_only=True);out={};global_min=float('inf');min_count=0;numeric_count=0
for sn,gene_col,first_data in [('b. 9,300 proteins',1,2),('c. 61,471 phosphosites',2,4)]:
 sh=w[sn];it=sh.iter_rows(values_only=True);header=list(next(it));rows=[]
 for row in it:
  if isinstance(row[gene_col],str) and 'CGNL1' in row[gene_col].split(';'):rows.append(list(row))
  if first_data==4:
   for v in row[first_data:]:
    if isinstance(v,(int,float)) and math.isfinite(v):
     numeric_count+=1
     if v<global_min:global_min=v;min_count=1
     elif v==global_min:min_count+=1
 out[sn]={'header':header,'target_rows':rows}
 print(sn,'target_rows',len(rows),flush=True)
res={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_file':str(p),'global_phosphosite_minimum':global_min,'minimum_occurrences':min_count,'phosphosite_numeric_values':numeric_count,'tables':out,'no_effect_contrasts_computed':True};(b/'metadata/ESCC_TARGET_EXTRACTION_382.json').write_text(json.dumps(res,indent=2),encoding='utf-8');print('floor',global_min,'count',min_count,flush=True)
