from pathlib import Path
import json,gzip,hashlib,datetime,os,fcntl,collections,xml.etree.ElementTree as E,base64,zlib,struct,traceback,sys
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1');F=P/'C04_full_fraction_qualification';W=F/'reporter_normalization_20261026'
ns={'m':'http://psi.hupo.org/ms/mzml'};tag='{'+ns['m']+'}spectrum'
targets=[127.124761,127.131081,128.128116,128.134436,129.131471,129.137790,130.134825,130.141145]
lock=(W/'extract.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as stream:
  for b in iter(lambda:stream.read(8388608),b''):h.update(b)
 return h.hexdigest()
def save(name,x):
 tmp=W/(name+'.tmp');tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8');tmp.replace(W/name)
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def arrays(node):
 out={}
 for b in node.findall('.//m:binaryDataArray',ns):
  acc={v.attrib['accession'] for v in b.findall('m:cvParam',ns)}
  if not acc&{'MS:1000514','MS:1000515'}:continue
  assert ('MS:1000523' in acc)!=('MS:1000521' in acc)
  raw=base64.b64decode(b.findtext('m:binary','',ns),validate=True);raw=zlib.decompress(raw) if 'MS:1000574' in acc else raw;fmt='d' if 'MS:1000523' in acc else 'f';size=8 if fmt=='d' else 4;assert len(raw)%size==0
  out['mz' if 'MS:1000514' in acc else 'intensity']=struct.unpack('<'+fmt*(len(raw)//size),raw)
 assert len(out['mz'])==len(out['intensity'])
 return out
state=dict(pid=os.getpid(),started=stamp(),files={},complete=False,whole_project_complete=False)
try:
 spec_sha=sha(W/'SPEC_BEFORE_REPORTER_EXTRACTION.json');inventory_sha=sha(W/'ELIGIBLE_PSM_INVENTORY.json')
 if '--pilot-check' not in sys.argv:assert json.loads((W/'EXTRACTION_PILOT_ACCEPTANCE.json').read_text())['pass_all']
 inventory=json.loads((W/'ELIGIBLE_PSM_INVENTORY.json').read_text())['rows'];assert len(inventory)==18805
 geneids=set()
 for line in (P/'C04_search_pilot/mouse_plus_contaminants.fasta').open():
  if line.startswith('>') and 'GN=Cgnl1 ' in line:geneids.add(line[1:].split()[0])
 assert 'sp|Q6AW69|CGNL1_MOUSE' in geneids
 names=sorted({x['file'] for x in inventory});assert len(names)==12
 save('EXTRACTION_INPUTS.json',dict(utc=stamp(),spec_sha256=spec_sha,inventory_sha256=inventory_sha,files=names,Cgnl1_protein_identifiers=sorted(geneids),scope_includes_pY_with_zero_eligible=True))
 allstats=[]
 for name in (names[:1] if '--pilot-check' in sys.argv else names):
  stem=Path(name).stem;destination=W/(stem+'.reporters.json.gz');done=W/(stem+'.receipt.json')
  if done.exists():
   old=json.loads(done.read_text());assert old['spec_sha256']==spec_sha and old['inventory_sha256']==inventory_sha and sha(destination)==old['output_sha256'];allstats.append(old);state['files'][name]=old;continue
  assert not destination.exists(), 'Incomplete prior output requires audit'
  state.update(current_file=name,stage='reading_metadata');save('EXTRACTION_STATUS.json',state)
  path=P/'C04_raw_pilot/phospho_fraction_1.mzML' if stem=='phospho_fraction_1' else F/'converted'/stem/(stem+'.mzML')
  expected='adab96e2d5db2be34eedc1dc6529c1799a65a5b72314693a7b65b44c1b609449' if stem=='phospho_fraction_1' else json.loads((F/'converted'/stem/'ACCEPTANCE.json').read_text())['mzML_sha256']
  assert sha(path)==expected
  meta={};index=collections.defaultdict(list);third=[]
  for _,node in E.iterparse(path,events=['end']):
   if node.tag!=tag:continue
   cv={v.attrib['accession']:v.attrib.get('value') for v in node.findall('m:cvParam',ns)};level=int(cv['MS:1000511'])
   if level in [2,3]:
    scan=int(node.attrib['id'].rsplit('scan=',1)[1]);refs=tuple(sorted({v.attrib['spectrumRef'] for v in node.findall('m:precursorList/m:precursor',ns) if 'spectrumRef' in v.attrib}))
    if level==3:third.append(dict(scan=scan,refs=refs))
    else:
     rt=next(v for v in node.findall('.//m:cvParam',ns) if v.attrib['accession']=='MS:1000016');ion=node.find('m:precursorList/m:precursor/m:selectedIonList/m:selectedIon',ns);vals={} if ion is None else {v.attrib['accession']:v.attrib.get('value') for v in ion.findall('m:cvParam',ns)}
     modes=sorted({v.attrib['name'] for v in node.findall('.//m:activation/m:cvParam',ns) if 'dissociation' in v.attrib.get('name','')});mode='HCD' if modes==['beam-type collision-induced dissociation'] else 'CID' if modes==['collision-induced dissociation'] else None
     d=dict(id=node.attrib['id'],scan=scan,refs=refs,rt=float(rt.attrib['value'])*(60 if rt.attrib.get('unitAccession')=='UO:0000031' else 1),mz=float(vals.get('MS:1000744','0')),charge=vals.get('MS:1000041'),mode=mode);meta[scan]=d;index[(refs,d['charge'],mode)].append(d)
   node.clear()
  seeds=[x for x in inventory if x['file']==name];psms={(x['mode'],int(x['scan'])):x for x in seeds};assert len(psms)==len(seeds)
  groups={};unlinked=[]
  for x in seeds:
   s=meta[int(x['scan'])];assert s['mode']==x['mode'] and int(s['charge'])==int(x['charge']);other='CID' if s['mode']=='HCD' else 'HCD'
   matches=[m for m in index[(s['refs'],s['charge'],other)] if s['refs'] and abs(m['mz']-s['mz'])<=s['mz']*1e-5 and abs(m['rt']-s['rt'])<=2]
   if len(matches)!=1:unlinked.append(dict(mode=x['mode'],scan=int(x['scan']),reason='opposite_activation_not_unique',matches=len(matches)));continue
   h=s if s['mode']=='HCD' else matches[0];c=s if s['mode']=='CID' else matches[0];key=(h['scan'],c['scan'])
   if key in groups:continue
   assigned=[psms[k] for k in [('HCD',h['scan']),('CID',c['scan'])] if k in psms];features={(v['plain_peptide'],v['phospho_count']) for v in assigned};conflict=len(features)!=1
   groups[key]=dict(HCD=h['scan'],CID=c['scan'],CID_id=c['id'],feature=list(next(iter(features))) if not conflict else None,assignment_conflict=conflict,Cgnl1_associated=any(geneids.intersection(v['protein'].split(',')) for v in assigned),passing_PSMs=[dict(mode=v['mode'],scan=int(v['scan']),charge=int(v['charge']),peptide=v['plain_peptide'],phospho_count=v['phospho_count'],q_global=float(v['q_global']),q_phospho=float(v['q_phospho'])) for v in assigned])
  bycid=collections.defaultdict(list)
  for g in groups.values():bycid[g['CID_id']].append(g)
  selected={}
  for x in third:
   found=[g for ref in x['refs'] for g in bycid.get(ref,[])]
   if found:selected[x['scan']]=found
  state['stage']='extracting_MS3_reporters';save('EXTRACTION_STATUS.json',state);records=[];seen=set()
  for _,node in E.iterparse(path,events=['end']):
   if node.tag!=tag:continue
   scan=int(node.attrib['id'].rsplit('scan=',1)[1])
   if scan in selected:
    v=arrays(node);peaks=[]
    for target in targets:
     hit=[(m,i) for m,i in zip(v['mz'],v['intensity']) if abs(m-target)<=target*20e-6]
     peaks.append(dict(target_mz=target,matches=len(hit),observed_mz=hit[0][0] if len(hit)==1 else None,intensity=hit[0][1] if len(hit)==1 and hit[0][1]>0 else None))
    gg=selected[scan];complete=all(x['intensity'] is not None for x in peaks);discordant=len(gg)!=1 or any(g['assignment_conflict'] for g in gg);cgnl=any(g['Cgnl1_associated'] for g in gg)
    records.append(dict(file=name,MS3=scan,groups=gg,reporters=peaks,all8_positive_unique=complete,assignment_ambiguous=discordant,Cgnl1_associated=cgnl,reference_eligible=complete and not discordant and not cgnl,feature=gg[0]['feature'] if not discordant else None));seen.add(scan)
   node.clear()
  assert seen==set(selected) and len(records)==len(seen)
  linkedkeys={(g['HCD'],g['CID']) for rr in records for g in rr['groups']}
  orphan=[dict(HCD=g['HCD'],CID=g['CID'],reason='no_explicit_MS3_reference') for key,g in groups.items() if key not in linkedkeys]
  with gzip.open(destination,'wt',encoding='utf-8') as stream:json.dump(dict(records=records,unlinked_PSMs=unlinked,orphan_pairs=orphan),stream,allow_nan=False)
  st=dict(file=name,utc=stamp(),spec_sha256=spec_sha,inventory_sha256=inventory_sha,mzML_sha256=expected,output_sha256=sha(destination),eligible_PSMs=len(seeds),unique_precursor_pairs=len(groups),unique_MS3=len(records),unlinked_PSMs=len(unlinked),orphan_pairs=len(orphan),reference_eligible_MS3=sum(x['reference_eligible'] for x in records),incomplete_MS3=sum(not x['all8_positive_unique'] for x in records),ambiguous_MS3=sum(x['assignment_ambiguous'] for x in records),Cgnl1_MS3=sum(x['Cgnl1_associated'] for x in records));save(stem+'.receipt.json',st);allstats.append(st);state['files'][name]=st;save('EXTRACTION_STATUS.json',state);print(json.dumps(st),flush=True)
 if '--pilot-check' in sys.argv:
  save('EXTRACTION_PILOT_EXIT.json',dict(exit=0,files=allstats,full_scope_complete=False));state.update(stage='pilot_complete_full_scope_pending');save('EXTRACTION_STATUS.json',state);raise SystemExit(0)
 assert len(allstats)==12 and sum(x['eligible_PSMs'] for x in allstats)==18805
 state.update(complete=True,stage='extraction_complete_normalization_pending',finished=stamp());save('EXTRACTION_STATUS.json',state);save('EXTRACTION_EXIT.json',dict(exit=0,files=allstats,new_biological_tests=0,whole_project_complete=False))
except Exception:
 state.update(stage='failed',error=traceback.format_exc());save('EXTRACTION_STATUS.json',state);save('EXTRACTION_FAILURE.json',state);raise
