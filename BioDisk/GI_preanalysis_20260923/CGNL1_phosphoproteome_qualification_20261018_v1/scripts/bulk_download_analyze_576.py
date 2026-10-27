from pathlib import Path
import os, sys, json, time, re, hashlib, zipfile, shutil, subprocess
import urllib.parse, xml.etree.ElementTree as ET
import concurrent.futures, datetime, posixpath, csv, io, socket

B = Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1')
D = Path('/mnt/BioData/bio/wangmingyuan/CGNL1_phosphosite_bulk_576')
M = B / 'metadata/bulk_576'
O = B / 'results/bulk_576'
for p in (D, M, O): p.mkdir(exist_ok=True, parents=True)

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def save(p, value):
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')
    tmp.replace(p)

def curl(url, path, timeout=3600, resume=False):
    # Data and signed URLs are never interpolated into shell commands.
    cmd = ['curl', '-q', '--location', '--fail', '--silent', '--show-error',
           '--connect-timeout', '60', '--max-time', str(timeout), '--speed-limit',
           '1024', '--speed-time', '300', '--max-filesize', '2147483648',
           '--write-out', '%{http_code}', '--output', str(path)]
    route = 'direct'
    try:
        with socket.create_connection(('127.0.0.1', 7897), timeout=1): pass
        cmd += ['--proxy', 'http://127.0.0.1:7897', '--noproxy', '']
        route = 'existing_user_proxy'
    except OSError: cmd += ['--noproxy', '*']
    if resume and path.exists() and path.stat().st_size: cmd += ['--continue-at', '-']
    r = subprocess.run(cmd + [url], capture_output=True, text=True)
    return dict(exit=r.returncode, http=r.stdout.strip(), error=r.stderr[-1800:], route=route)

def graphql(query):
    cmd = ['curl', '-q', '--noproxy', '*', '--fail', '--silent', '--show-error',
           '--connect-timeout', '30', '--max-time', '120', '-H', 'Content-Type: application/json',
           '--data-binary', '@-', 'https://pdc.cancer.gov/graphql']
    r = subprocess.run(cmd, input=json.dumps({'query': query}), capture_output=True, text=True)
    if r.returncode: raise RuntimeError('PDC metadata: ' + r.stderr[-300:])
    j = json.loads(r.stdout)
    if j.get('errors'): raise RuntimeError(str(j['errors'])[:600])
    return j['data']

def analyze_file(path, cohort):
    """Read-only target extraction, preserving headers and all patient values.
    Statistical inference waits for an explicit verified schema; unknown formats
    are reported, never silently interpreted as a valid biological comparison.
    """
    result = {'file': str(path), 'cohort': cohort, 'utc': now(), 'tables': []}
    if path.suffix.lower() == '.xlsx':
        ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        with zipfile.ZipFile(path) as z:
            strings = []
            if 'xl/sharedStrings.xml' in z.namelist():
                with z.open('xl/sharedStrings.xml') as f:
                    for _, e in ET.iterparse(f, events=('end',)):
                        if e.tag.endswith('}si'):
                            strings.append(''.join(t.text or '' for t in e.iter() if t.tag.endswith('}t')))
                            e.clear()
            def val(c):
                v = c.find('m:v', ns)
                if c.get('t') == 'inlineStr': return ''.join(t.text or '' for t in c.iter() if t.tag.endswith('}t'))
                if v is None: return None
                if c.get('t') == 's': return strings[int(v.text)]
                if c.get('t') in ('str', 'e', 'b'): return v.text
                try: return float(v.text)
                except (ValueError, TypeError): return v.text
            rels = {x.get('Id'): x.get('Target') for x in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
            for sheet in ET.fromstring(z.read('xl/workbook.xml')).findall('m:sheets/m:sheet', ns):
                target = rels[sheet.get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id')]
                member = target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/' + target)
                count = 0; heads = []; hits = []
                with z.open(member) as f:
                    for _, e in ET.iterparse(f, events=('end',)):
                        if not e.tag.endswith('}row'): continue
                        count += 1; cells = e.findall('m:c', ns)
                        first = [val(c) for c in cells[:8]]
                        hit = any(isinstance(v,str) and any(k in v.upper() for k in ('CGNL1','Q0VF96','ENSG00000128849')) for v in first)
                        if count <= 5 or hit:
                            row = {'row': int(e.get('r', count)), 'cells': {c.get('r'):val(c) for c in cells if val(c) is not None}}
                            if count <= 5: heads.append(row)
                            if hit: hits.append(row)
                        e.clear()
                result['tables'].append(dict(sheet=sheet.get('name'), rows=count, headers=heads, target_rows=hits))
    elif path.suffix.lower() in ('.tsv','.csv','.txt'):
        hits=[]; heads=[]; count=0
        with path.open(encoding='utf-8-sig') as f:
            for count, line in enumerate(f,1):
                if count<=5: heads.append(line.rstrip('\n'))
                if any(k in line.upper() for k in ('CGNL1','Q0VF96','ENSG00000128849')): hits.append({'line':count,'text':line.rstrip('\n')})
        result['tables'].append(dict(rows=count,headers=heads,target_rows=hits))
    else:
        result['status']='source_document_or_unsupported_table_format'
    if result['tables']:
        result['target_rows_found']=sum(len(t['target_rows']) for t in result['tables'])
        result['status']='target_extracted_pairing_and_scale_review_required' if result['target_rows_found'] else 'no_target_in_this_file_not_proof_of_unmeasured_cohort'
    dest=O/cohort;dest.mkdir(exist_ok=True)
    save(dest/(hashlib.sha256(str(path).encode()).hexdigest()[:12]+'.json'),result)
    return {'file':str(path),'status':result['status'],'target_rows':result.get('target_rows_found')}

def unpack(path, cohort, depth=0):
    if depth>3: raise RuntimeError('Nested archive depth exceeded')
    if not zipfile.is_zipfile(path): return [path]
    with zipfile.ZipFile(path) as z:
        if 'xl/workbook.xml' in z.namelist(): return [path]
        if sum(x.file_size for x in z.infolist())>8*1024**3: raise RuntimeError('Archive expanded size needs review')
        assert z.testzip() is None, 'Archive CRC failed'
        dest=path.parent/(path.name+'_contents');dest.mkdir(exist_ok=True)
        files=[]
        for item in z.infolist():
            if item.is_dir():continue
            n=Path(item.filename)
            if n.suffix.lower() not in ('.xlsx','.xls','.csv','.tsv','.txt','.pdf','.doc','.docx','.zip'):continue
            p=(dest/n).resolve()
            if dest.resolve() not in p.parents:raise RuntimeError('Unsafe archive path')
            if (item.external_attr>>16)&0o170000==0o120000:raise RuntimeError('Archive symlink not accepted')
            p.parent.mkdir(parents=True,exist_ok=True)
            if not p.exists():
                with z.open(item) as src,p.open('wb') as dst:shutil.copyfileobj(src,dst)
            if p.suffix.lower()=='.zip':files+=unpack(p,cohort,depth+1)
            else:files.append(p)
        return files

def download(item, cohort, state, statepath):
    folder=D/cohort;folder.mkdir(exist_ok=True)
    name=re.sub(r'[^A-Za-z0-9._-]','_',item['name']);final=folder/name
    if final.exists() and final.with_suffix(final.suffix+'.receipt.json').exists():return final
    attempts=[]
    for attempt in range(1,4):
        url=item.get('url')
        if item.get('study'):
            q='{ filesPerStudy(pdc_study_id:'+json.dumps(item['study'])+',file_name:'+json.dumps(item['name'])+',offset:0,limit:10){file_id signedUrl{url}} }'
            rows=graphql(q)['filesPerStudy'];row=next(x for x in rows if x['file_id']==item['file_id']);v=row['signedUrl'];url=v['url'] if isinstance(v,dict) else v[0]['url']
        part=folder/(name+('.part' if item.get('md5') else '.attempt'+str(attempt)+'.part'))
        state.update(stage='downloading',current_file=name,attempt=attempt,utc=now());save(statepath,state)
        r=curl(url,part,resume=bool(item.get('md5')));r.update(attempt=attempt,utc=now(),bytes=part.stat().st_size if part.exists() else 0);attempts.append(r)
        save(folder/(name+'.attempts.json'),attempts)
        if r['exit']==0:
            if item.get('size') is not None:assert part.stat().st_size==int(item['size']), 'Size mismatch'
            sha=hashlib.sha256();md=hashlib.md5()
            with part.open('rb') as f:
                for chunk in iter(lambda:f.read(8*1024*1024),b''):sha.update(chunk);md.update(chunk)
            if item.get('md5'):assert md.hexdigest()==item['md5'],'Publisher MD5 mismatch'
            if item.get('archive') or name.endswith(('.zip','.xlsx')):
                assert zipfile.is_zipfile(part),'Expected ZIP/XLSX; received other content'
                with zipfile.ZipFile(part) as z:assert z.testzip() is None,'ZIP CRC failed'
            else:
                with part.open('rb') as f:prefix=f.read(512).lower()
                assert b'<html' not in prefix and b'<!doctype html' not in prefix,'HTML received instead of table'
            assert not final.exists(),'Refuse overwrite';part.rename(final)
            save(final.with_suffix(final.suffix+'.receipt.json'),dict(utc=now(),source=item,bytes=final.stat().st_size,sha256=sha.hexdigest(),md5=md.hexdigest(),attempts=attempts))
            return final
        if r['http']=='404' or r['exit'] in (33,63):break
        if r['http'] in ('401','403') and attempt>=2:break
        if attempt<3:
            state.update(stage='waiting_retry',next_retry_utc=datetime.datetime.fromtimestamp(time.time()+900,datetime.timezone.utc).isoformat(),utc=now());save(statepath,state);time.sleep(900)
    raise RuntimeError('Download failed: '+name+' '+json.dumps(attempts[-1]))

def discover(job):
    items=list(job.get('items',[]));notes=[]
    for accession in job.get('px',[]):
        f=M/(accession+'_manifest.xml')
        result=curl('https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID='+accession+'&outputMode=XML',f,timeout=180)
        if result['exit']:
            notes.append(dict(accession=accession,result=result));continue
        try: tree=ET.parse(f)
        except ET.ParseError as err:
            notes.append(dict(accession=accession,error=str(err)));continue
        found=[];raw_count=0
        for node in tree.iter():
            name=node.get('name','');u=node.get('value','')
            if 'raw file URI' in name:raw_count+=1
            if not u.startswith(('https://','http://','ftp://')):continue
            fn=Path(urllib.parse.urlsplit(u).path).name
            if name in ('Search engine output file URI','Associated raw file URI','Other file URI') and fn.lower().endswith(('.tsv','.csv','.xlsx','.txt','.zip')) and not any(k in fn.lower() for k in ('.raw.','raw_files','rawdata')):
                found.append(dict(name=fn,url=u.replace('ftp://','https://',1),archive=fn.lower().endswith(('.zip','.xlsx'))))
        items.extend(found);notes.append(dict(accession=accession,processed_file_candidates=len(found),raw_file_count=raw_count,raw_requantification_may_be_required=not found and raw_count>0))
    for study in job.get('pdc',[]):
        allrows=[]
        for offset in range(0,20000,1000):
            q='{filesPerStudy(pdc_study_id:'+json.dumps(study)+',offset:'+str(offset)+',limit:1000){file_id file_name file_size md5sum}}'
            rows=graphql(q)['filesPerStudy'];allrows+=rows
            if len(rows)<1000:break
        save(M/(study+'_files.json'),allrows)
        try:
            q='{biospecimenPerStudy(pdc_study_id:'+json.dumps(study)+'){case_submitter_id sample_submitter_id sample_type aliquot_submitter_id}}'
            save(M/(study+'_biospecimens.json'),graphql(q))
        except Exception as err:
            notes.append(dict(study=study,biospecimen_error=repr(err)))
        for row in allrows:
            name=row['file_name']
            if name.endswith(('.sample.txt','.summary.tsv','.tmt10.tsv','.tmt11.tsv','.itraq.tsv','.xlsx')) and not any(k in name for k in ('.peptides.','.phosphopeptide.','.qcmetrics.')):
                items.append(dict(study=study,file_id=row['file_id'],name=name,size=row['file_size'],md5=row['md5sum']))
        notes.append(dict(study=study,metadata_rows=len(allrows),pagination_complete=len(rows)<1000))
    for page in job.get('landing_pages',[]):
        f=D/job['id']/('landing_'+hashlib.sha256(page.encode()).hexdigest()[:8]+'.html');f.parent.mkdir(exist_ok=True)
        result=curl(page,f,timeout=120);notes.append(dict(page=page,result=result))
        if result['exit']!=0:continue
        text=f.read_text(encoding='utf-8')
        urls=re.findall(r'(?:href|src)=["\x27]([^"\x27]+)["\x27]',text)
        for u in urls:
            u=urllib.parse.urljoin(page,u.replace('&amp;','&'))
            if urllib.parse.urlsplit(u).scheme not in ('http','https'):continue
            if re.search(r'\.(xlsx|tsv|csv|zip)(?:\?|$)',u,re.I) and any(k in u.lower() for k in ('supp','moesm','mmc','download','data')):
                items.append(dict(name=Path(urllib.parse.urlsplit(u).path).name,url=u,archive=u.lower().split('?')[0].endswith(('.zip','.xlsx'))))
    unique={json.dumps(x,sort_keys=True):x for x in items}
    return list(unique.values()),notes

def run(job):
    path=M/(job['id']+'.json');state=dict(id=job['id'],cancer=job['cancer'],cohort=job['cohort'],stage='discovering',utc=now(),files=[],analyses=[],errors=[]);save(path,state)
    try:
        reused=job.get('reuse',[])
        for x in reused:
            if Path(x).exists():state['files'].append(dict(path=x,reused=True))
            else:state['errors'].append('Missing reused evidence '+x)
        if job.get('completed_evidence'):
            for x in job['completed_evidence']:assert Path(x).exists(),x
            state.update(stage='existing_analysis_retained',analysis_evidence=job['completed_evidence'],limitation=job.get('limitation'));save(path,state);return state
        items,notes=discover(job);state['discovery_notes']=notes;state['planned_files']=len(items);save(path,state)
        for item in items:
            try:
                f=download(item,job['id'],state,path);state['files'].append(dict(path=str(f),verified=True));state.update(stage='analyzing_downloaded_tables',utc=now());save(path,state)
                for table in unpack(f,job['id']):
                    try:state['analyses'].append(analyze_file(table,job['id']))
                    except Exception as err:state['errors'].append('Table '+str(table)+': '+repr(err))
                save(path,state)
            except Exception as err:state['errors'].append(repr(err));save(path,state)
        if not items:state['stage']='source_discovery_requires_followup'
        elif state['errors']:state['stage']='partial_completion_followup_required'
        else:state['stage']='download_and_target_extraction_complete'
        state['next_action']='Resolve source scale, peptide localization, patient mapping, and count independent donors; run frozen paired analysis where schema is qualified. Never infer differential P from unknown or aggregate tables.'
    except Exception as err:state.update(stage='source_discovery_requires_followup');state['errors'].append(repr(err))
    state['utc']=now();save(path,state);return state

if __name__=='__main__':
    plan=json.loads((M/'PLAN.json').read_text())
    save(M/'RUN.json',dict(pid=os.getpid(),started=now(),cohorts=len(plan['jobs']),concurrency=3))
    # Preserve completed prior analyses while newly obtained files immediately
    # undergo source/target analysis instead of waiting for all downloads.
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        results=list(pool.map(run,plan['jobs']))
    save(M/'EXIT.json',dict(utc=now(),results=[dict(id=x['id'],stage=x['stage'],errors=x['errors']) for x in results],scientific_project_complete=False))

