from pathlib import Path
import importlib.util,concurrent.futures,json
B=Path("/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1")
spec=importlib.util.spec_from_file_location("bulk",B/"scripts/bulk_download_analyze_576.py");m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
m.M=B/"metadata/missing579";m.M.mkdir(exist_ok=True)
jobs=[dict(id="B01",cancer="CCA",cohort="Deng217",items=[dict(name=f"Supplementary_Table_{i}.xlsx",url=f"https://zenodo.org/records/7708351/files/Supplementary_Table_{i}.xlsx?download=1",archive=True) for i in range(1,8)]),dict(id="C03",cancer="CRC",cohort="Guo144",items=[dict(name="PMC12686132_supplementary.zip",url="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12686132/supplementaryFiles",archive=True)])]
m.save(m.M/"PLAN.json",jobs)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:result=list(p.map(m.run,jobs))
m.save(m.M/"EXIT.json",result)
