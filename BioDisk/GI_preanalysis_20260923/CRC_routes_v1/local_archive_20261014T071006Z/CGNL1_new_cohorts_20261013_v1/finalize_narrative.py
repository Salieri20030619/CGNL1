from pathlib import Path
p=Path(__file__).resolve().parent/'07_report.py';s=p.read_text(encoding='utf-8')
s=s.replace('<h1>两个新队列未复现校正后的内部增量，部分跨队列预测改善仍保留</h1>','<h1>增殖跨队列预测出现双向增量，队列内部独立性尚未稳定复现</h1>')
s=s.replace('两队列在身份与成熟分化控制后均未改善留一患者预测；两者之间的增殖模型迁移出现正向误差改善，完整保留并分解其来源。','两者之间的增殖模型迁移出现双向正向误差改善，且去均值误差也下降，形成值得保留的候选预测线索。两队列在身份与成熟分化控制后均未改善内部留一患者预测，目标效应也未达到名义显著；两类证据并列呈现。')
s=s.replace('<p>新队列之间增殖预测的正向改善不能遗漏。以下为已冻结预测的事后描述性误差分解，','<p>两新队列之间，增殖预测的原始MSE分别改善5.4%与32.2%，去均值误差也分别降低20.7%与17.0%。因此，正向预测结果不全由整体均值偏移解释。GSE200997到GSE132465方向的总MSE减少主要来自偏差降低，同时也包含患者间误差下降；反向迁移中平均偏差变差，但去均值误差下降抵消了这一变化。</p><p>这些结果支持保留跨队列预测的候选线索；其证据属性是预测误差的描述性改善，不称为名义显著或稳定独立效应。应与两队列内部M1目标效应及留一患者预测结果一起解释。以下为已冻结预测的事后描述性误差分解，')
s=s.replace('Positive values indicate additional out-of-patient information','Positive values indicate lower held-out prediction error')
s=s.replace("new_primary_to_primary_cycling_transfer=decomp", "interpretation='Positive bidirectional cycling transfer includes centered-error improvement; within-cohort conditional inference and LOPO do not establish stable independence. Retain exploratory predictive clue.',new_primary_to_primary_cycling_transfer=decomp")
p.write_text(s,encoding='utf-8');compile(s,str(p),'exec');print('BALANCED_REPORT_READY')
