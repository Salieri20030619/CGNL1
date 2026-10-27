from pathlib import Path
p=Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2\scripts\second_validation.R');s=p.read_text();s=s.replace("sets<-list(","m$strict_colon<-m$strict_colon%in%c('True','TRUE');sets<-list(");p.write_text(s,encoding='ascii')
