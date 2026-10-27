import os,sys,json,time,hashlib,pathlib,subprocess,concurrent.futures,urllib.request,urllib.error
ROOT=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923')
def fetch(item):
    name,url=item; dest=ROOT/'01_metadata'/name; t=time.time()
    r={'name':name,'url':url,'started':time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    try:
        cmd=['curl','-L','--fail','--connect-timeout','15','--max-time','100','--retry','1','-o',str(dest)+'.part',url]
        p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        r.update(exit_code=p.returncode,stderr=p.stderr[-600:])
        if p.returncode: raise RuntimeError('curl '+str(p.returncode))
        os.replace(str(dest)+'.part',dest)
        r.update(bytes=dest.stat().st_size,sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),status='downloaded')
    except Exception as e:r.update(status='failed',error=str(e))
    r['elapsed_seconds']=round(time.time()-t,2)
    (ROOT/'logs'/(name+'.download.json')).write_text(json.dumps(r,indent=2))
    print(json.dumps(r),flush=True)
    return r
if __name__=='__main__':
    items=[]
    for acc in ['GSE145356','GSE294559']:
        items.append((acc+'.soft.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=all&form=text&view=full'))
        items.append((acc+'.files.html',f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/'))
    items += [
      ('GC_ascites_2023.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9929081/fullTextXML'),
      ('CRC_multiome_2025.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12541563/fullTextXML'),
      ('zenodo_11058999.json','https://zenodo.org/api/records/11058999'),
      ('PXD064711.json','https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD064711'),
      ('PXD064711_files.json','https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession=PXD064711&pageSize=100'),
      ('CRC_persister_search.json','https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TITLE:%22Colorectal%20cancer%20cells%20possess%22&format=json'),
      ('GC_niche_search.json','https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:10.1053/j.gastro.2024.08.007&format=json')]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:res=list(ex.map(fetch,items))
    (ROOT/'logs'/'initial_acquisition.json').write_text(json.dumps(res,indent=2))
    (ROOT/'logs'/'initial_acquisition.complete').write_text('completed; inspect per-file statuses\n')
