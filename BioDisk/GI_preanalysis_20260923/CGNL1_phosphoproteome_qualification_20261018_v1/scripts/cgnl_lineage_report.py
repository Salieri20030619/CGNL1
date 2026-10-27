from pathlib import Path
import json,statistics,hashlib,datetime
P=Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1');O=P/'continuous_exploration_20261023_v1';D=O/'material_search_restart_20261024';A=D/'lineage'
def read(f):return json.loads(f.read_text())
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def save(f,z):f.write_text(json.dumps(z,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');assert read(f)==z
r=read(A/'RESULTS.json');cs=read(A/'CONTRASTS.json');obs=read(A/'OBSERVATIONS.json');mapping=read(A/'SAMPLE_MAPPING.json')
assert len(mapping)==149 and len(set((x['case'],x['tissue']) for x in mapping))==149
assert len(obs)==3129 and sum(x['bcm_missing']!=x['umich_missing'] for x in obs)==0
baseline=next(x for x in read(O/'wave1/RESULTS.json') if x['test']=='C01_A')
mean=statistics.mean(x['bcm'] for x in cs);other=statistics.mean(x['umich'] for x in cs)
assert len(cs)==baseline['n']==28 and abs(mean-baseline['mean'])<1e-12
agree=sum((x['bcm']>0)==(x['umich']>0) for x in cs);assert agree==28
text=f'''CGNL1 PDAC官方处理版本对应核查（2026-10-24）
结论：当前C01的28对患者模式在UMich与BCM两个官方版本中方向一致，但处理流程不完全等价，且S198仍须按其源肽段解释。

身份与缺失模式
使用已归档的PDC生物样本表，以病例、组织类型和aliquot/sample标识符连接；未使用数值相似度猜测身份。149个组织样本均有唯一对应，没有跳过或有歧义的病例-组织组合。两个导出共有21个CGNL1位点，合计3129个样本-位点组合，其中1214个双方均有数值，缺失模式不一致数为0。
这直接支持两份导出在当前CGNL1子集具有相同覆盖。它不能证明每个非缺失值都来自直接观测，也不能替代逐谱图定位概率或上游插补掩码。

数值对应
所有可比较值的BCM−UMich中位差为{r['global_offset']:.12f}。减去该值后，最大残差为{r['max_global_offset_residual']:.12f}；所以不能将两者描述为只有一个严格相同的加法偏移。
C01 S198−S297有28对完整患者，BCM均值{mean:.12f}，UMich均值{other:.12f}；28/28患者方向一致，最大患者对比差为{r['S198_S297_max_residual']:.12f}。BCM计算与此前主分析均值一致至1e-12以内。本次是来源一致性核查，没有重新计算或挑选P值，也没有把同队列重处理计为独立生物学验证。

解释与后续重点
新增证据缩小了“处理版本是否改变C01方向”的不确定性。当前更关键的是测量对象：UMich S198行的PSNsQPTsPSLEDPAK同时标注S198与S202，不能将位点汇总结果自动解释为孤立S198磷酸化占有率。S297行源肽段与S198行也不具有完全相同的磷酸化状态构成。
仍需核查TMT-Integrator对多磷酸化肽段进入single-site汇总的规则、位点/肽段层面的贡献，以及逐观测定位证据。两个包README软件版本不同，数值接近不消除版本差异。MaxPepProb不能当作位点定位概率。
C01暂保留为“含S198的磷酸肽/位点汇总信号相对于S297的变化”候选；不更改旧数据或原始统计。取得汇总与定位依据后，再决定是否能缩窄到孤立位点假设。
C02/C03既有结论与关键阴性保持；本轮未取得新的C04材料。Goal保持active，下一步继续处理S198/S202肽段状态与源汇总规则。

复核位置
服务器material_search_restart_20261024/lineage下保存预先限定的比较方案、149例样本映射、3129条比较记录、28例对比、结果和验收。源数据和此前分析均保留。
'''
f=D/'CGNL1_PDAC_LINEAGE_20261024.txt';f.write_text(text,encoding='utf-8');assert f.read_text()==text
save(A/'ACCEPTANCE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),unique_mapping_count=149,missing_mask_agreement=True,baseline_C01_mean_reproduced=True,contrast_sign_agreement=agree,not_exact_additive_transform=True,new_tests=0,source_scope='CGNL1 shared rows only',artifacts=[dict(path=str(x),sha256=sha(x)) for x in [A/'SPEC_BEFORE_COMPARISON.json',A/'RESULTS.json',A/'OBSERVATIONS.json',A/'SAMPLE_MAPPING.json',A/'CONTRASTS.json']],report_sha256=sha(f)))
s=read(O/'STATUS.json');s.update(C01='UMich/BCM149samplemapping and3129missingmasks agree;28contrastsignsagree;S198/S202compositepeptide andperPSMlocalization unresolved',latest_lineage_report=str(f),next_task='TMTIntegrator single-site multi-phosphopeptide rollup and S198/S202 source contributions');save(O/'STATUS.json',s)
print(json.dumps(dict(report=str(f),sha256=sha(f),n=28,mean_bcm=mean,mean_umich=other)))
