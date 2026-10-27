from common import *
D=Q/'ASSEMBLE';D.mkdir(exist_ok=True);sections=['# 后台队列证据汇总','运行版本：'+C['id'],'本文件汇总实际计算状态。细胞身份、临床来源及外部图件仍需复核；探索性RNA组成模型不能解释为测得真实细胞比例或因果证明。'];results={}
for task in C['tasks']:
 name=task['id']
 if name=='ASSEMBLE':continue
 p=Q/name/'RESULT.json';r=json.loads(p.read_text()) if p.exists() else {'status':'no_result','reason':'check task EXIT and dependency status'};results[name]=r;sections.extend(['## '+name,'```json',json.dumps(r,ensure_ascii=False,indent=2),'```'])
(D/'后台队列证据汇总.md').write_text('\n\n'.join(sections),encoding='utf-8-sig');dump(D/'RESULT.json',{'status':'assembled_pending_scientific_review','tasks':results})
