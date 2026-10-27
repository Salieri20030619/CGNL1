from pathlib import Path
import json,gzip,tarfile,io,xml.etree.ElementTree as ET,requests,concurrent.futures,hashlib
O=Path.cwd()/'mechanism_investigation_v2';s=(O/'metadata/GSE308231_family.soft').read_bytes().decode('utf-8');bad=[l for l in s.splitlines() if l.startswith('!Sample_data_processing') and any(x in l for x in ['PEDV','Intestinal mucosal'])];p=O/'raw/GSE308231_RAW.tar.preserved_prefix';header_rows=[]
# Qualify only a complete feature-file member: no expression measurements from a partial archive.
with tarfile.open(p,mode='r|') as tf:
 for m in tf:
  if m.name.endswith('_features.tsv.gz'):
   compressed=tf.extractfile(m).read();assert len(compressed)==m.size;raw=gzip.decompress(compressed);lines=raw.decode('utf-8').splitlines();header_rows=[l.split('\t') for l in lines];feature={'member':m.name,'compressed_bytes':m.size,'gzip_crc_verified':True,'member_sha256':hashlib.sha256(compressed).hexdigest(),'rows':len(lines),'human_ENSG_ids':sum(r[0].startswith('ENSG') for r in header_rows),'pig_ENSSSCG_ids':sum(r[0].startswith('ENSSSCG') for r in header_rows),'target_rows':[r for r in header_rows if 'CGNL1' in r],'first_rows':header_rows[:4]};break
 else:raise RuntimeError('No complete feature member in preserved prefix')
out={'issue':'GC-Q07','conflicting_processing_lines':bad,'feature_qualification':feature,'interpretation':'Humanfeatureidentifierssupporthumanreferenceannotation; theydonotindependentlyprovegastrictissueprovenance. PEDVprocessingtextisaconfirmedmetadata inconsistency,notproofsamplesarepigs.','expression_analysis_from_partial_archive':False};(O/'metadata/GSE308231_SPECIES_AND_METADATA_106.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(out,ensure_ascii=False)[:6000])
root=ET.fromstring((O/'literature/PMC13181027.xml').read_bytes());links=[e.attrib for e in root.iter() if e.tag in ['media','supplementary-material']];(O/'metadata/GSE314812_SUPPLEMENT_LINKS.json').write_text(json.dumps(links,indent=2),encoding='utf-8');print('SUPPLEMENT_LINKS',links)
u='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13181027/supplementaryFiles';a=requests.get(u,timeout=(15,90));a.raise_for_status();(O/'metadata/GSE314812_article_supplement.zip').write_bytes(a.content);print('Supplementdownload',len(a.content),flush=True)
