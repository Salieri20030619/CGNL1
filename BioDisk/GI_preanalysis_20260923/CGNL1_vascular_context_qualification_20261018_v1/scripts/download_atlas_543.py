from pathlib import Path
import subprocess,json,hashlib,datetime
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_vascular_context_qualification_20261018_v1");D=Path("/mnt/BioData/bio/wangmingyuan/CGNL1_vascular_context_20261018/normal_vascular_atlas_543");D.mkdir(exist_ok=True)
files=[("lec",47622685,"f4a17b4e6e4ce1b10e3a437d8a0093b8"),("ec",358632551,"ec0a1e7f31921def7304019c37e03a63")];records=[]
for kind,size,md5 in files:
 name="final_"+kind+"_oct2024_cxg.h5ad";url="https://cellgeni.cog.sanger.ac.uk/vasculature/"+name;p=D/(name+".part");dest=D/name
 assert not dest.exists()
 r=subprocess.run(["curl","-q","--proxy","http://127.0.0.1:7897","--noproxy","","-L","--fail","-sS","--connect-timeout","15","--max-time","900","-o",str(p),url]);assert r.returncode==0
 assert p.stat().st_size==size;h=hashlib.md5()
 with p.open("rb") as f:
  for c in iter(lambda:f.read(8*1024*1024),b""):h.update(c)
 assert h.hexdigest()==md5;p.rename(dest);records.append(dict(kind=kind,path=str(dest),bytes=size,actual_md5=h.hexdigest(),expected_md5=md5,checksum_provenance="Publisher x-amz-meta-s3cmd-attrs MD5; NOT multipart ETag"));(B/"metadata/ATLAS543_STATE.json").write_text(json.dumps(dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=records,complete=len(records)==2),indent=2),encoding="utf-8");print(json.dumps(records[-1]),flush=True)
