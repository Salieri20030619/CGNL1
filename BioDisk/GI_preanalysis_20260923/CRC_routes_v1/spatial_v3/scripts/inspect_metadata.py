from pathlib import Path
import re,json,pandas as pd
O=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3');p=O/'scripts/audit_spatial_inputs.py';s=p.read_text(encoding='utf-8-sig').replace("(O/'metadata/BJC2026.html').read_text()","(O/'metadata/BJC2026.html').read_text(encoding='utf-8')");p.write_text(s,encoding='utf-8')
s=(O/'metadata/BJC2026.html').read_text(encoding='utf-8');links=sorted(set(re.findall(r'https?[^\s<>"\']+\.(?:xlsx|pdf|docx)',s)));print('SUPPLEMENTS',json.dumps(links))
m=pd.read_csv(O.parent/'metadata/spatial_patient_map.tsv',sep='\t');v=m[m.visium&m.sample_type.str.contains('Liver')];print('LIVER',len(v),'sections',v.patient.nunique(),'patients');v.to_csv(O/'results/GSE294385_liver_sections.tsv',sep='\t',index=False)
j=json.loads((O.parent/'metadata/spatial_table2.json').read_text(encoding='utf-8'));print('CLINICAL KEYS',list(j));print('CLINICAL PREVIEW',str(j.get('Tab 3',[]))[:3000])
s=(O/'metadata/iss_panel_2.tmap').read_text(encoding='utf-8');print('TMAPSTART',s[:300]);print('CGNL1_IN_TMAP','CGNL1' in s)
