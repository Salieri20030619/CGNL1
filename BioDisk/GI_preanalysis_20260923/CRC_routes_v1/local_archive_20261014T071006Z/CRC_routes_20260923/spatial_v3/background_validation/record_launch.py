from pathlib import Path
import json,subprocess,datetime,shutil,sys
sys.stdout.reconfigure(encoding='utf-8')
B=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3\background_validation');info=json.loads((B/'LATEST_RUN.json').read_text());R=Path(info['local_run']);remote=info['remote_run'];python='/mnt/BioDisk/wangmingyuan/miniconda3/envs/pdac_spatial/bin/python'
r=subprocess.run(['ssh','sgs10',python,remote+'/scripts/status.py',remote],capture_output=True,text=True,encoding='utf-8',check=True);status=json.loads(r.stdout);now=datetime.datetime.now(datetime.timezone.utc).isoformat();(R/'STARTUP_CHECK.json').write_text(json.dumps({'checked_at_utc':now,'status':status},ensure_ascii=False,indent=2),encoding='utf-8')
j=json.loads((B/'tasks.json').read_text());j['status']='first_wave_launched';j['active_run']=info['run_id'];state={'V01':'computed_pending_review','V02':'input_preparation_running','V03':'waiting_reference_eligibility','V04':'input_preparation_running','V05':'waiting_evidence_gate'}
for x in j['jobs']:x['status']=state[x['id']]
(B/'tasks.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8')
out=Path(r'C:\Users\31937\Downloads\CGNL1空间初筛_20260924');out.mkdir(exist_ok=True);progress=status['jobs']['DATA_PREP'].get('download_progress.json',{});pid=status['jobs']['DATA_PREP']['STARTED.json']['pid'];v01=status['jobs']['V01'];mark=v01.get('EXIT.json',{});process=json.loads((R/'LOCAL_BACKGROUND_PROCESS.json').read_text(encoding='utf-8-sig'))
text=f'''# 后台验证运行状态

核验时间（UTC）：{now}。这是启动时快照，实时状态以服务器文件为准。

运行版本：{info['run_id']}

- V01严格边缘验证：退出码{mark.get('exit_code')}，完成标记已生成；结果为computed_pending_review，尚未作为正式验证结论发布。
- V02/V04数据准备：正在运行；服务器任务PID {pid}。启动核验时{progress.get('successful',0)}/{progress.get('expected',110)}个文件通过完整性校验，已验证{progress.get('validated_bytes',0):,}字节。
- V03解卷积：等待参考细胞身份和临床资格确认。
- V05公开干预机制响应：等待前述证据关卡，不预先认定某条轴成立。

后台网络任务由隐藏SSH进程维持，本机PID {process['pid']}，无需保持当前对话打开；数据下载仍依赖本机联网及本机代理。关闭本机或代理会影响下载，不会把中断误记为完成。

服务器目录：{remote}

查看入口：V01/run.log、V01/summary.tsv、V01/EXIT.json、DATA_PREP/run.log、DATA_PREP/download_progress.json、DATA_PREP/EXIT.json。任务成功时另有COMPLETE.json；进程消失本身不视为成功。

本地运行记录：{R}

当前完成的是首批任务启动。单细胞数据下载/矩阵QC不是恶性细胞身份确认，外部数据取得也不等于空间复现成功。后续阶段通过数据资格与来源审核后再推进。
'''
(R/'运行状态.md').write_text(text,encoding='utf-8-sig');(out/'后台运行状态.md').write_text(text,encoding='utf-8-sig');shutil.copy2(B/'LATEST_RUN.json',out/'后台最新运行.json');shutil.copy2(R/'STARTUP_CHECK.json',out/'后台启动核验.json')
subprocess.run(['scp','-r','sgs10:'+remote+'/V01',str(R)],check=True)
subprocess.run(['scp',str(B/'tasks.json'),'sgs10:'+remote.rsplit('/',1)[0]+'/'],check=True)
subprocess.run(['scp',str(R/'STARTUP_CHECK.json'),str(R/'运行状态.md'),'sgs10:'+remote+'/'],check=True)
print(json.dumps({'run':info['run_id'],'download_progress':progress,'status_file':str(out/'后台运行状态.md')},ensure_ascii=False,indent=2))
