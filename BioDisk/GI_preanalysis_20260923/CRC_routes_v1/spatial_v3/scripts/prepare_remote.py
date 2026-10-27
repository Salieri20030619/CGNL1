from pathlib import Path
import json,sys,openpyxl
sys.stdout.reconfigure(encoding='utf-8')
O=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3')
p=O/'scripts/check_gene_coverage.py';s=p.read_text(encoding='utf-8-sig').split('\nw=openpyxl')[0].replace(',openpyxl','').replace("O=Path(r'C:\\Users\\31937\\CRC_routes_20260923\\spatial_v3')","O=Path(sys.argv[1])")
(O/'scripts/remote_coverage.py').write_text(s,encoding='utf-8')
w=openpyxl.load_workbook(O/'metadata/BJC2026_gene_panels.xlsx',read_only=True,data_only=True);out={}
for sh in w:
 cells=[str(v) for row in sh.iter_rows(values_only=True) for v in row if v is not None];out[sh.title]={'rows':sh.max_row,'columns':sh.max_column,'CGNL1_matches':[v for v in cells if 'CGNL1' in v.upper() or 'ENSG00000128849' in v.upper()],'first_cells':cells[:10]}
(O/'results/BJC2026_probe_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False))
