from pathlib import Path
import json,html,hashlib
b=Path(r'C:\Users\31937\CGNL1_method_scout_20261006')
d=json.loads((b/'method_matrix.json').read_text(encoding='utf-8'))
intro=(b/'review_intro.txt').read_text(encoding='utf-8')
lines=[intro]
cards=[]
labels={'evidence_depth':'核查深度','observed_methods':'原研究方法','proposed_transfer':'本项目迁移建议','required_inputs':'所需输入','priority':'取舍','availability_and_caveats':'可得性与边界'}
for r in d['records']:
 title=f"{r['id']} | {r['publication_date']} | {r['venue']} | {r['title']}"
 lines.extend([title,r['url']]+[labels[k]+'：'+r[k] for k in labels]+[''])
 body=''.join('<p><strong>'+html.escape(labels[k])+'：</strong>'+html.escape(r[k])+'</p>' for k in labels)
 cards.append('<article id="'+r['id']+'"><h2>'+html.escape(r['id']+' · '+r['title'])+'</h2><p class="meta">'+html.escape(r['publication_date']+' | '+r['venue'])+'</p>'+body+'<p><a href="'+html.escape(r['url'],quote=True)+'">原始论文/官方来源</a></p></article>')
(b/'CGNL1_method_review_20261006.txt').write_text('\n'.join(lines),encoding='utf-8')
style='body{font-family:Segoe UI,Microsoft YaHei,sans-serif;max-width:1100px;margin:40px auto;padding:0 24px;color:#192b3c;background:#f7f9fb;line-height:1.8}header,article,.intro{background:white;border:1px solid #dbe4ed;border-radius:8px;padding:22px 28px;margin:18px 0}h1{font-size:27px}h2{font-size:19px;line-height:1.5}.intro{white-space:pre-wrap}.meta{color:#526b80}a{color:#00688b}p{margin:8px 0}@media print{body{margin:0;background:white}article{break-inside:avoid}}'
doc='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CGNL1 方法调研与迁移判断 · 2026-10-06</title><style>'+style+'</style></head><body><header><h1>CGNL1 方法调研与迁移判断</h1><p>截至2026-10-06 · 20项重点研究 · 定向方法调研</p></header><section class="intro">'+html.escape(intro)+'</section>'+''.join(cards)+'</body></html>'
(b/'CGNL1_method_review_20261006.html').write_text(doc,encoding='utf-8')
assert len(d['records'])==20
assert len({r['id'] for r in d['records']})==20
assert all(r['url'].startswith('https://') for r in d['records'])
assert all(all(r.get(k) for k in labels) for r in d['records'])
assert '\ufffd' not in doc
names=['CGNL1_method_review_20261006.html','CGNL1_method_review_20261006.txt','method_matrix.json','search_scope.json']
manifest={'date':'2026-10-06','status':'method_scout_complete','new_patient_analysis_started':False,'checks':['20 unique records','source links present','transfer/input/caveat fields present','UTF-8 checked'],'files':{name:{'bytes':(b/name).stat().st_size,'sha256':hashlib.sha256((b/name).read_bytes()).hexdigest()} for name in names}}
(b/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'records':20,'files':names,'checks':'passed'},ensure_ascii=False))
