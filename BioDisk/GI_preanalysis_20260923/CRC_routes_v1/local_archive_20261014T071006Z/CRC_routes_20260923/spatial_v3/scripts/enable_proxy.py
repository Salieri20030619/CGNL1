from pathlib import Path
p=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3\scripts\download_liver.py');s=p.read_text(encoding='utf-8-sig');s=s.replace("['curl','-L'","['curl']+(['-x',sys.argv[2]] if len(sys.argv)>2 else [])+['-L'");p.write_text(s,encoding='utf-8');compile(s,str(p),'exec')
