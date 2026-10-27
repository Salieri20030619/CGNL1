from pathlib import Path
p=Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2\scripts\prepare_bulk.py');s=p.read_text(encoding='utf-8-sig');s=s.replace("f=pd.read_csv(gzip.GzipFile(fileobj=tar.extractfile(x)),sep='\\t');assert f.iloc[:,0].is_unique;", "f=pd.read_csv(gzip.GzipFile(fileobj=tar.extractfile(x)),sep='\\t');bad=f.iloc[:,0].astype(str).str.match(r'^\\d{1,2}[-/](Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)',case=False);f=f.loc[~bad];assert f.iloc[:,0].is_unique;")
p.write_text(s,encoding='utf-8')
