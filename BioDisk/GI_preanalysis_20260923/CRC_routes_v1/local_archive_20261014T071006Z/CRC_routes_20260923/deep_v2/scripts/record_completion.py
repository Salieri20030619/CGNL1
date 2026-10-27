from pathlib import Path
import json,shutil
O=Path(r'C:\Users\31937\CRC_routes_20260923\deep_v2');D=O/'delivery';out=Path(r'C:\Users\31937\Downloads\CGNL1立项建议_20260924')
p=D/'交付完成.json';j=json.loads(p.read_text(encoding='utf-8'));j['remote_sync']='verified_sha256_pdf_and_zip';j['remote_delivery']='/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CRC_routes_v1/deep_v2/delivery';p.write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8');shutil.copy2(p,out/p.name)
(O/'RESUME.md').write_text('''# 项目恢复点（2026-09-24）

深化预分析及客户版立项建议已完成，交付目录delivery；Downloads/CGNL1立项建议_20260924已复制；sgs10同项目deep_v2/delivery已同步并核对PDF与ZIP的SHA256。

本轮完成：两独立bulk队列验证、单细胞RNA推断CNV敏感性、P001 RNA+ATAC检测证据、文献定位、机制相关性、8页客户PDF、完整审计、229文件证据包。18项统计复核通过，PDF逐页检查通过。

主候选CGNL1：GSE50760 18对中16升，效应+1.031，P=1.48e-4，全基因组FDR=0.00532；GSE14297严格结肠13对中12升，效应+1.251（VST单位），P=3.12e-6，FDR=0.000997。CNV95% 3/4升、99% 4/4升；小样本只作探索性。

关键未决：正常肝表达高，肝代理评分校正后证据减弱；独立癌细胞内来源验证尚未完成。机制程序无FDR阳性；不能称CGNL1–肌球蛋白/CAMSAP3轴已验证。原发侵袭前沿无一致增高。

研究定位：有条件立项，目标约10分，不保证发表。下一实质工作应优先独立单细胞/肿瘤上皮数据中的来源确认及治疗分层，再推进具体机制分支；不得重跑已完成bulk分析。湿实验尚未开展。
''',encoding='utf-8')
