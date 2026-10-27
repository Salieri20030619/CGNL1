from pathlib import Path
import json,subprocess,hashlib,datetime
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1");D=Path("/mnt/BioData/bio/wangmingyuan/CGNL1_vascular_context_20261018/normal_vascular_atlas_543");p=D/"final_ec_oct2024_cxg.h5ad.part";s=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),prior_exit=18,preserved_prefix_bytes=p.stat().st_size,stage="proxy_resume")
def save():(B/"metadata/ATLAS546_STATE.json").write_text(json.dumps(s,indent=2),encoding="utf-8")
save()
try:
 r=subprocess.run(["curl","-q","--proxy","http://127.0.0.1:7897","--noproxy","","-L","--fail","-sS","--connect-timeout","15","--max-time","1800","-C","-","-o",str(p),"https://cellgeni.cog.sanger.ac.uk/vasculature/final_ec_oct2024_cxg.h5ad"]);s["curl_exit"]=r.returncode;assert r.returncode==0;assert p.stat().st_size==358632551
 h=hashlib.md5()
 with p.open("rb") as f:
  for c in iter(lambda:f.read(8*1024*1024),b""):h.update(c)
 s["actual_md5"]=h.hexdigest();assert h.hexdigest()=="ec0a1e7f31921def7304019c37e03a63";dest=D/"final_ec_oct2024_cxg.h5ad";assert not dest.exists();p.rename(dest);s.update(stage="verified",path=str(dest));save();old=json.loads((B/"metadata/ATLAS543_STATE.json").read_text());old["files"].append(dict(kind="ec",path=str(dest),bytes=358632551,actual_md5=h.hexdigest(),expected_md5=h.hexdigest(),checksum_provenance="Publisher x-amz-meta-s3cmd-attrs MD5"));old["complete"]=True;(B/"metadata/ATLAS543_STATE.json").write_text(json.dumps(old,indent=2),encoding="utf-8");print(json.dumps(s))
except Exception as e:s.update(stage="failed",error=repr(e));save();raise
