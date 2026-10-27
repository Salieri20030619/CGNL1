from pathlib import Path
import json,subprocess,datetime,shutil
B=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3\background_validation');info=json.loads((B/'LATEST_RUN.json').read_text());R=Path(info['local_run']);remote=info['remote_run'];py='/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python'
r=subprocess.run(['ssh','sgs10',py,remote+'/scripts/status.py',remote],capture_output=True,text=True,encoding='utf-8',check=True);s=json.loads(r.stdout);d=s['jobs']['DATA_PREP'];assert d['EXIT.json']['exit_code']==0 and d['download_progress.json']['successful']==110
(R/'LATEST_STATUS.json').write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8')
j=json.loads((B/'tasks.json').read_text());j['status']='first_wave_completed_pending_review';states={'V01':'computed_pending_review','V02':'inputs_QC_ready_pending_identity_review','V03':'waiting_reference_eligibility','V04':'inputs_ready_pending_external_pathology_review','V05':'waiting_evidence_gate'}
for x in j['jobs']:x['status']=states[x['id']]
(B/'tasks.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8')
text=f'''# 后台验证运行状态

最新核验：{datetime.datetime.now(datetime.timezone.utc).isoformat()}（UTC）。

运行版本：{info['run_id']}

- 下载已完成：110/110个文件，1,854,917,901字节（约1.85 GB），全部通过文件完整性检查。
- 36份人类单细胞样本完成矩阵基础检查；细胞身份及临床资格尚待审核。
- 外部Visium文件已获取并检查压缩包完整性；尚未完成外部空间复现。
- 数据准备任务于2026-09-24 13:34:36（台北时间）正常结束，退出码0，COMPLETE.json存在。
- V01严格边缘验证已完成计算，结果待统一复核。
- V03解卷积和V05机制响应验证尚未启动。

该批次没有后台计算仍在运行，SSH下载通道已随任务结束退出。

服务器目录：{remote}

成功依据为退出码、完成标记及结果文件，而不是进程消失。
'''
out=Path(r'C:\Users\31937\Downloads\CGNL1空间初筛_20260924');(out/'后台运行状态.md').write_text(text,encoding='utf-8-sig');(R/'运行状态.md').write_text(text,encoding='utf-8-sig')
subprocess.run(['scp',str(B/'tasks.json'),'sgs10:'+remote.rsplit('/',1)[0]+'/'],check=True)
subprocess.run(['scp',str(R/'运行状态.md'),'sgs10:'+remote+'/'],check=True)
print('Download completion verified and status updated')
