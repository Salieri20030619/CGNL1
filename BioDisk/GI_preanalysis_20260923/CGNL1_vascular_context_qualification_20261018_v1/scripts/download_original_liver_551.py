from pathlib import Path
import subprocess,json,hashlib,datetime,time,os
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1");D=Path("/mnt/BioData/bio/wangmingyuan/CGNL1_vascular_context_20261018/original_liver_551");D.mkdir(exist_ok=True);out=B/"results/original_liver_551";out.mkdir(exist_ok=True)
manifest=json.loads((B/"metadata/original_liver_549/mendeley.response").read_text());entry=next(x for x in manifest if x["filename"]=="endothelial_cells_annot.rds");info=entry["content_details"];p=D/(entry["filename"]+".part");dest=D/entry["filename"];state=dict(pid=os.getpid(),source=info["download_url"],expected_sha256=info["sha256_hash"],expected_bytes=info["size"],stage="proxy_download",attempts=[])
def save():
 state["utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat();(B/"metadata/ORIGINAL551_STATE.json").write_text(json.dumps(state,indent=2),encoding="utf-8")
save()
try:
 assert not dest.exists()
 for attempt in range(1,4):
  command=["curl","-q","--proxy","http://127.0.0.1:7897","--noproxy","","-L","--fail","-sS","--connect-timeout","15","--max-time","3600","-o",str(p)]
  if p.exists() and p.stat().st_size:command += ["-C","-"]
  result=subprocess.run(command+[info["download_url"]]);state["attempts"].append(dict(number=attempt,exit=result.returncode,bytes=p.stat().st_size if p.exists() else 0));save()
  if result.returncode==0:break
  if result.returncode not in [18,28,56]:raise RuntimeError("Non-recoverable curl exit "+str(result.returncode))
  if attempt<3:time.sleep(15)
 assert result.returncode==0 and p.stat().st_size==info["size"]
 h=hashlib.sha256()
 with p.open("rb") as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b""):h.update(chunk)
 state["actual_sha256"]=h.hexdigest();assert h.hexdigest()==info["sha256_hash"];p.rename(dest);state.update(stage="verified_reading_metadata",verified_path=str(dest));save()
 with (out/"metadata_read.log").open("w",encoding="utf-8") as log:
  result=subprocess.run(["/mnt/BioDisk/wangmingyuan/miniconda3/envs/r_seurat/bin/Rscript","--vanilla",str(B/"scripts/inspect_original_liver_551.R"),str(dest),str(out/"METADATA551.json"),str(out/"cell_metadata551.tsv")],stdout=log,stderr=subprocess.STDOUT)
 state["metadata_read_exit"]=result.returncode;assert result.returncode==0 and (out/"METADATA551.json").exists();state["stage"]="ready_for_source_comparison";save();print(json.dumps(state),flush=True)
except Exception as e:state.update(stage="failed",error=repr(e));save();raise
