import pathlib,sys,json,importlib.metadata,hashlib
from PIL import Image,ImageDraw
import pandas as pd,numpy as np
O=pathlib.Path(sys.argv[1]);R=O/'results/screen';jobs=json.loads((O/'metadata/LIVER_DOWNLOAD_PLAN.json').read_text());rows=[]
for j in jobs:
 p=O/'raw/liver'/(j['sample']+'.tar.gz');m=json.loads((O/'raw/liver'/(j['sample']+'.complete.json')).read_text());assert m['ok'] and p.stat().st_size==m['bytes'];rows.append(m)
(O/'results/LIVER_INPUTS_VERIFIED.json').write_text(json.dumps({'files':rows,'n':len(rows),'total_bytes':sum(x['bytes'] for x in rows),'status':'all24 individually validated'},indent=2))
s=pd.read_csv(R/'summary.tsv',sep='\t');expected={'secondary_genes3':3,'baseline_distance_sensitivity':8,'all_timepoints_distance_sensitivity':12,'baseline_source_comparisons':12,'all_timepoints_source_comparisons':12,'baseline_score_sensitivity':24,'all_timepoints_score_sensitivity':24};qc=[]
for fam,z in s.groupby('family'):
 if fam=='single_primary':assert z.FDR.isna().all();continue
 n=expected[fam];v=np.r_[z.P.fillna(1).values,np.ones(n-len(z))];order=np.argsort(v);q=np.minimum.accumulate((v[order]*n/np.arange(1,n+1))[::-1])[::-1];out=np.zeros(n);out[order]=np.minimum(q,1);delta=np.abs(out[:len(z)]-z.FDR.values);err=float(np.nanmax(delta)) if np.isfinite(delta).any() else 0;assert err<1e-10;qc.append({'family':fam,'fixed_hypotheses':n,'max_error':err})
(R/'FDR_RECHECK.json').write_text(json.dumps({'pass':True,'families':qc},indent=2));(R/'SOFTWARE.json').write_text(json.dumps({'python':sys.version,'packages':{x:importlib.metadata.version(x) for x in ['numpy','pandas','scipy','matplotlib','Pillow']}},indent=2))
files=[R/'figures'/(j['sample']+'.png') for j in jobs]
for start in range(0,len(files),4):
 im=Image.new('RGB',(1960,720),(255,255,255));draw=ImageDraw.Draw(im)
 for k,p in enumerate(files[start:start+4]):
  x=Image.open(p).convert('RGB');x.thumbnail((980,345));im.paste(x,((k%2)*980,(k//2)*360))
 im.save(R/'figures'/f'contact_{start//4+1}.jpg',quality=92)
# Add decision-critical findings up front in the concise report, without removing full results.
f=O/'空间初筛结果.md';body=f.read_text(encoding='utf-8-sig');needle='## 预设主要终点';intro='''## 当前立项判断

空间初筛通过“区域关联”这一关，尚未通过“癌细胞内独立变化”这一关。建议保留CGNL1作为局部肝环境相关候选，下一步优先做细胞分辨率来源与响应验证；当前不应直接确定肝因素—癌细胞CGNL1—连接骨架的因果轴。

- 主分析：11名患者9升2降，差值+0.541，95%CI 0.120至0.963，P=0.0169。逐患者剔除后均值仍为正（+0.429至+0.632）。该主终点单独预设，不作FDR校正。
- 200/400微米敏感性：11名患者10升1降，差值+0.529，P=0.00219，距离敏感性家族FDR=0.00876。400/800微米仅9名有效患者，8升1降，差值+0.574，保持探索性。
- 关键来源结果：肝小叶及转移基质区域均在11/11名患者中高于肿瘤区；平均差分别+1.361（P=0.000350，FDR=0.00105）和+1.113（P=0.000109，FDR=0.000434）。FDR在4候选×3来源比较的12项家族内计算。
- 关键敏感性：以spot等权均值代替区域伪bulk后，主距离差值+0.584，P=0.0643、敏感性FDR=0.386；加入组成评分及测序深度后系数+0.182，4升7降，P=0.289、敏感性FDR=0.861。因此效应方向存在区域层面支持，但对汇总方式和组成评分敏感。
- 同片微/宏转移可比较5名患者，4名微转移更高，平均差+0.357；只作为探索线索，不作推断性检验。

空间结果增加的是“信号与局部肝/基质邻域相关”的证据。它不能仅凭校正前后P值决定混杂或肝诱导机制，更不能将高表达spot直接归为纯癌细胞。现阶段最有价值的后续判别是：排除混合边缘spot后，具有独立恶性细胞身份的CGNL1信号是否仍随邻近肝/基质环境改变。

''';assert needle in body;f.write_text(body.replace(needle,intro+needle,1),encoding='utf-8-sig');print(json.dumps({'input_files':len(rows),'bytes':sum(x['bytes'] for x in rows),'FDR_families_rechecked':len(qc),'contact_sheets':6},indent=2))
