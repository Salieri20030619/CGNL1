from pathlib import Path
L=Path(__file__).resolve().parent;p=L/'07_report.py';s=p.read_text(encoding='utf-8')
s=s.replace("transfer=pd.read_csv(O/'N3_transfer_performance.tsv',sep='\\t');liver=", "transfer=pd.read_csv(O/'N3_transfer_performance.tsv',sep='\\t');decomp=pd.read_csv(O/'N3_transfer_error_decomposition.tsv',sep='\\t');liver=")
s=s.replace("body=f'''<header>", """transfer_note=''
for _,z in decomp[(decomp.train!='old_discovery23')&decomp.outcome.eq('cycling')].iterrows():
 transfer_note+=f'<p>{z.train} → {z.test}：原始MSE改善{100*z.improvement:+.1f}%；去均值误差改善{100*z.centered_improvement:+.1f}%。平方偏差减少量={z.bias_squared_reduction:.5g}，去均值误差减少量={z.centered_error_reduction:.5g}。</p>'
body=f'''<header>""")
s=s.replace('<h1>新增两个原发CRC队列，检验上皮状态增量能否推广</h1>','<h1>两个新队列未复现校正后的内部增量，部分跨队列预测改善仍保留</h1>')
s=s.replace('结果用于评价一般上皮状态信息，不能直接替代肝定植特异性验证。','两队列在身份与成熟分化控制后均未改善留一患者预测；两者之间的增殖模型迁移出现正向误差改善，完整保留并分解其来源。结果用于评价一般上皮状态信息。')
s=s.replace('<section><h2>条件模型和患者外预测共同界定增量范围</h2>','<section><h2>两队列均未改善校正后的内部预测，跨队列增殖预测出现正向结果</h2>')
s=s.replace('<p>跨队列预测完全使用训练来源拟合；', '<p>提高上皮覆盖门槛到50细胞后，GSE200997保留11位患者，M1增殖系数转为正（β=+0.455，P=0.666、FDR=0.666），负向关系对覆盖门槛敏感；GSE132465在两种门槛下均保留22位患者，主结论一致。身份敏感性与全部门槛结果保留在完整统计表。</p><p>跨队列预测完全使用训练来源拟合；')
s=s.replace("'''+table(transfer)+'''</section>", "'''+table(transfer)+'''<p>新队列之间增殖预测的正向改善不能遗漏。以下为已冻结预测的事后描述性误差分解，恒等式为MSE＝平方平均误差＋去均值误差；没有用验证集重新拟合或校准，也不产生新的P值。</p>'''+transfer_note+table(decomp[decomp.outcome.eq('cycling')])+'''</section>")
s=s.replace("('跨队列预测','results/N3_transfer_performance.tsv')","('跨队列预测','results/N3_transfer_performance.tsv'),('迁移误差分解','results/N3_transfer_error_decomposition.tsv')")
s=s.replace("liver_small_cohort_descriptive_only=True,next_priority=","new_primary_to_primary_cycling_transfer=decomp[(decomp.train!='old_discovery23')&decomp.outcome.eq('cycling')].to_dict('records'),liver_small_cohort_descriptive_only=True,next_priority=")
p.write_text(s,encoding='utf-8');compile(s,str(p),'exec');print('REPORT_REVISED')
