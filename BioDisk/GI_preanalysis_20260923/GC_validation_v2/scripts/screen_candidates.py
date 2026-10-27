import pathlib,json,numpy as np,pandas as pd
B=pathlib.Path('/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923');V=B/'GC_validation_v2';O=V/'candidate_discovery';O.mkdir(exist_ok=True)
(O/'PLAN.md').write_text('''# Candidate discovery after FABP4 pause\nScope: gastric cancer peritoneal metastasis, no fixed cell lineage or mandatory lipid mechanism. FABP4 paused. CRC direction remains paused.\n1. All-gene patient-paired screen: genome-wide FDR <0.05, abs(logFC)>=1, AveExpr>=1, >=13/17 patients in effect direction, median paired direction concordant.\n2. Preserve all genes and models; assess GC05 exclusion, adipocyte adjustment, adipocyte+myeloid and broad stroma+immune models. Composition adjustment ranks evidence, not a causal test or mandatory elimination of microenvironment processes.\n3. Cell-source and within-cell-type patient pseudobulk checks using existing data. Missing cell types and n<10 remain exploratory; no cell-level p-values.\n4. Primary-paper novelty audit, including original dataset publication and aliases. Understudied must be context-specific, not inferred from no single search hit.\n5. Shortlist 2-3 only if evidence supports it; otherwise record rejections and unresolved data needs. No claim of functional causality before experiments.\n''',encoding='utf-8')
a=pd.read_csv(B/'04_results/GC_bulk/integrated_sensitivity_gene_table.tsv',sep='\t')
for k in ['adipocyte','adipocyte_myeloid']:
 z=pd.read_csv(V/'04_results/FABP4_audit'/f'voom_{k}_all_genes.tsv.gz',sep='\t')[['gene_id','logFC','P.Value','adj.P.Val']].rename(columns={c:c+'_'+k for c in ['logFC','P.Value','adj.P.Val']});a=a.merge(z,on='gene_id',validate='one_to_one')
a['consistent_patients']=np.where(a.logFC>0,a.positive_pairs,a.negative_pairs)
a['screen_pass']=(a['adj.P.Val']<.05)&(a.logFC.abs()>=1)&(a.AveExpr>=1)&(a.consistent_patients>=13)&(a.logFC*a.median_observed_paired_log2CPM_difference>0)
for k in ['composition','exclude05','adipocyte','adipocyte_myeloid']:
 fc='logFC_'+k;fdr='adj.P.Val_'+k
 a['same_direction_'+k]=a.logFC*a[fc]>0
 a['supported_'+k]=a['same_direction_'+k]&(a[fdr]<.05)
a['n_supported_models']=a[['supported_'+k for k in ['composition','exclude05','adipocyte','adipocyte_myeloid']]].sum(axis=1)
a['absolute_effect']=a.logFC.abs();a['paused_FABP4']=a.symbol.eq('FABP4')
a['structural_flag']=a.symbol.fillna('').str.match(r'^(IG[HKL][VDJC]|MT-|RP[SL][0-9])')
a.to_csv(O/'all_gene_screen.tsv.gz',sep='\t',index=False)
s=a[a.screen_pass].sort_values(['n_supported_models','consistent_patients','adj.P.Val'],ascending=[False,False,True]);s.to_csv(O/'screen_pass_all.tsv',sep='\t',index=False)
cols=['symbol','logFC','adj.P.Val','consistent_patients','n_supported_models','adj.P.Val_adipocyte','adj.P.Val_composition','structural_flag']
print('SCREEN_PASS',len(s),'UP',sum(s.logFC>0),'DOWN',sum(s.logFC<0));print(s.loc[~s.structural_flag,cols].head(100).to_string(index=False))
(O/'SCREEN_COMPLETE.json').write_text(json.dumps({'exit_code':0,'tested':len(a),'pass':len(s)}))
