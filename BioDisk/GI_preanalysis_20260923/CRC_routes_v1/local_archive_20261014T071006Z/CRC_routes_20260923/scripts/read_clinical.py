import pathlib,subprocess,json,openpyxl
O=pathlib.Path(r'C:\Users\31937\CRC_routes_20260923');M=O/'metadata';p=M/'spatial_table2.xlsx';print('XLSX_MAGIC',p.read_bytes()[:4]);w=openpyxl.load_workbook(p,read_only=True,data_only=True);out={}
for s in w:
 rows=list(s.values);out[s.title]=rows;print('SHEET',s.title,s.max_row,s.max_column);print(rows[:5])
(M/'spatial_table2.json').write_text(json.dumps(out,default=str,ensure_ascii=False,indent=2),encoding='utf-8')
u='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10275599/supplementaryFiles';r=subprocess.run(['curl.exe','-x','http://127.0.0.1:7897','-L','--fail','--max-time','120','-o',str(M/'GSE225857_supplements.zip'),u],capture_output=True);print('SUPP_ZIP',r.returncode)
