from pathlib import Path
import json
B=Path(r'C:\Users\31937\CRC_routes_20260923\spatial_v3\background_validation');info=json.loads((B/'LATEST_QUEUE.json').read_text());Q=Path(info['local']);mapping={'common':'common','eligibility':'eligibility','source_reference':'source_reference','external_qc':'external_qc','mixlib':'mixlib','mixture':'mixture','external_spatial':'external_spatial','mechanism_gate':'mechanism_gate','assemble':'assemble'}
for src,dst in mapping.items():
 p=B/(src+'_w2.py');s=p.read_text(encoding='utf-8-sig')
 if src=='mixlib':s=s.replace("ix=pd.Index(sy).get_indexer(genes[gi])","index={g:i for i,g in enumerate(sy) if unique[i]};ix=np.array([index[g] for g in genes[gi]])")
 compile(s,dst,'exec');(Q/'scripts'/(dst+'.py')).write_text(s,encoding='utf-8')
c=json.loads((Q/'CONFIG.json').read_text());c['mixture']['primary_region_mask']='V01 5-of-6 mixed-edge-excluded mask; original mask sensitivity';c['mixture']['anchor_quality_gate']='>=50 distant normal hepatic spots and median fitted hepatic RNA weight>=0.4';c['mixture']['target_capture_calibration']='per-sample target scaling using only hepatic spots>=600um from macro tumor, >=10 observed target molecules';c['external_model_regions']={'hepatic_RNA_weight_min':.7,'CNV_epithelial_RNA_weight_min':.6,'hepatic_weight_max_in_tumor_like':.1,'fibroblast_weight_max_in_tumor_like':.2,'requires_morphology_review':True};c['mechanism_expansion_gate']='strict residual n>=5, >=70% positive, positive mean; FDR<0.05 if n>=10; both external model regions positive. Inventory stage only; treatment design review required.';(Q/'CONFIG.json').write_text(json.dumps(c,indent=2))
print('9 worker modules staged and syntax checked')
