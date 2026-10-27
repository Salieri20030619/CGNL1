from pathlib import Path
import subprocess,json,datetime,hashlib,zipfile,os
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1");D=Path("/mnt/BioData/bio/wangmingyuan/CGNL1_vascular_context_20261018/independent_nmp_534");D.mkdir(exist_ok=True)
state=dict(pid=os.getpid(),stage="proxy_download",source="https://zenodo.org/api/records/16939324/files/data.zip/content",expected_md5="dd8885428162d717b8f9486f658d451f",expected_bytes=1029552662)
def save():
 state["utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat();(B/"metadata/NMP534_STATE.json").write_text(json.dumps(state,indent=2),encoding="utf-8")
save();p=D/"data.zip.part"
try:
 assert not (D/"data.zip").exists()
 r=subprocess.run(["curl","-q","--proxy","http://127.0.0.1:7897","--noproxy","","-L","--fail","-sS","--connect-timeout","15","--max-time","1800","--output",str(p),state["source"]]);state["curl_exit"]=r.returncode;assert r.returncode==0
 assert p.stat().st_size==state["expected_bytes"]
 h=hashlib.md5()
 with p.open("rb") as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b""):h.update(chunk)
 state["actual_md5"]=h.hexdigest();assert h.hexdigest()==state["expected_md5"];p.rename(D/"data.zip");state["stage"]="extracting_single_cell";save()
 with zipfile.ZipFile(D/"data.zip") as z:
  member="CRLM_SingleCell/adata_CRLM_NMP.h5ad";assert member in z.namelist();z.extract(member,D)
 state.update(stage="verified_single_cell_ready_for_qualification",matrix_path=str(D/member),zip_md5_pass=True,single_cell_zip_crc_pass=True);save();print(json.dumps(state),flush=True)
except Exception as e:
 state.update(stage="failed",error=repr(e));save();raise
