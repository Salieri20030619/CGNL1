from pathlib import Path
import importlib.util,subprocess,json,zipfile,hashlib,datetime,traceback
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1")
spec=importlib.util.spec_from_file_location("bulk",B/"scripts/bulk_download_analyze_576.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
# Isolated version: retain original batch metadata and raise the limit only for this file.
src=(B/"scripts/bulk_download_analyze_576.py").read_text().replace("'2147483648'","'107374182400'").replace("timeout=3600","timeout=14400")
exec(compile(src,"bulk_H02_578","exec"),m.__dict__)
p=B/"metadata/H02_RESTART578.json"
st={"id":"H02","stage":"starting","files":[],"errors":[]}
item={"name":"Phospho_MaxQuant_results.zip","url":"https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/03/PXD008373/Phospho_MaxQuant_results.zip","archive":True}
try:
 f=m.download(item,"H02",st,p)
 st.update(stage="verified_archive",file=str(f))
 with zipfile.ZipFile(f) as z:
  members=[dict(name=i.filename,size=i.file_size) for i in z.infolist()]
  st["members"]=members
  dest=f.parent/"Phospho_MaxQuant_results_selected578";dest.mkdir(exist_ok=True)
  selected=[]
  for i in z.infolist():
   if Path(i.filename).name.lower() not in ["phospho (sty)sites.txt","proteingroups.txt","summary.txt","parameters.txt"]:continue
   target=(dest/i.filename).resolve()
   assert dest.resolve() in target.parents
   if not target.exists():z.extract(i,dest)
   selected.append(str(target))
  st["extracted"]=selected
 st["stage"]="complete_verified_selected_tables"
except Exception as e:st.update(stage="failed",error=repr(e),traceback=traceback.format_exc())
st["utc"]=m.now();m.save(p,st)
