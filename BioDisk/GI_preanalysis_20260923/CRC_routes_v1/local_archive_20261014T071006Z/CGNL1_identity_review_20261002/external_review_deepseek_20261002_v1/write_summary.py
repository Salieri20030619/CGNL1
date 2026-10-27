from pathlib import Path
import json,hashlib,tarfile,datetime
R=Path(r'C:\Users\31937\CGNL1_identity_review_20261002\external_review_deepseek_20261002_v1')
report='''DeepSeek匿名审阅接收与三方比较说明

一、接收与版本
已保存DeepSeek原始Markdown并计算SHA256，另存44行结构化表。编号、顺序与原匿名模板完全一致。
DeepSeek明确声明未接触原标签、自动参考预测或既往结论。此独立性按审阅者声明记录；所交Markdown没有输入文件哈希或外部冻结清单，无法单靠文件独立验证其会话环境。原网页版本已披露历史上下文暴露，不能将两份结果统称为“两份严格独立盲审”。
本轮没有覆盖任一审阅原稿、已有簇或细胞注释，也没有重跑目标基因关联。

二、主要结论
统一同义词和主要谱系粒度后，43/44簇判断相符。这个数包括D01/C006一致缺乏确定谱系，且不表示置信度、亚型或细胞纯度全部一致。
D01/C018：两份审阅均为高置信度血管内皮，对应GSE285991原簇14。
D02/C028：两份审阅均为高置信度血管内皮，对应GSE245552原簇21。
这支持两个簇的主要内皮身份，不等于每个细胞均纯净，也不证明CGNL1相关效应、内皮亚型或功能。

三、对原自动参考的重要修正线索
D01/C013：两份审阅均支持上皮身份，而原双参考均给T/NK。双参考一致不能代替标记和审阅证据。
D01/C011：两份审阅支持pDC，原参考标签未单列这一类而给B；应按pDC候选保留。
D01/C012及D02/C041：两份审阅支持肥大细胞，原粗参考没有充分覆盖，不应继承T/NK等最近类标签。
D02/C043：基质/成纤维身份一致，但与壁细胞的边界及肌成纤维亚型仍需额外表达证据。

四、需要保守处理的审阅措辞
D02/C003：DeepSeek主要标签为“增殖细胞”，网页为T细胞；给定标记中CD3D/E提供直接T谱系证据。建议身份与状态分开，记录“T细胞；增殖状态”，并保留原始分歧。
D01/C006：两份资料都没有足够谱系证据，继续保留“谱系未定；增殖”，不能自动归T/NK。
D01/C039：未进入前30的CD3不能当作阴性，不应据此认定纯NK。保留“细胞毒淋巴细胞，NK倾向”。
D01/C035：FOXP3?不能算已观测证据；T身份明确，Treg亚型保留倾向。
D02/C038：上皮来源与完整细胞资格分开，需补原始条形码、背景及细胞级QC。
D02/C042：粒细胞可呈低RNA特征；低UMI本身不足以判为质量差、环境RNA或双细胞。
供者集中本身也不是双细胞或混合簇证据，不按某个审阅者的“高风险”描述自动删除。

五、证据逐项核对
支持基因简写展开后共核对551项。D01/C024的CD68、LYZ不在提供的该簇前30标记中，不能记作该材料已提供的支持；这不等于基因实际不表达。FOXP3?单独标为不确定项，不计实测证据。两个内皮簇的所列支持基因均能对应到提供材料。
DeepSeek文件未明确声明收到匿名QC附注。已确认的D02线粒体字段问题仍有效：当前特征匹配无法计算该指标，原0值不能解释为实测低线粒体比例。

六、后续使用范围
“簇级共识建议_未回填.tsv”是新的建议版本，保留了原始标签、分歧、信心和依据。待必要的细胞级身份及QC核查后，再用于固定的CGNL1关联复测；仍按目标原始P<0.05筛选正式阳性，完整保留FDR和关键阴性。
保持较宽的候选与支持充分的集合两个层次，不根据目标显著性改注释，不机械沿用单标记/双标记数量规则。
'''
(R/'联合审阅比较说明.txt').write_text(report,encoding='utf-8')
manifest=[]
for p in R.rglob('*'):
 if p.is_file() and p.name!='DELIVERY_MANIFEST.json' and not p.name.endswith('.tar.gz'):manifest.append({'path':p.relative_to(R).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(R/'DELIVERY_MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
archive=R/'deepseek_review_comparison_20261002_v1.tar.gz'
with tarfile.open(archive,'w:gz') as t:
 for entry in manifest:t.add(R/entry['path'],arcname=entry['path'])
 t.add(R/'DELIVERY_MANIFEST.json',arcname='DELIVERY_MANIFEST.json')
print('REPORT',R/'联合审阅比较说明.txt');print('ARCHIVE',archive)