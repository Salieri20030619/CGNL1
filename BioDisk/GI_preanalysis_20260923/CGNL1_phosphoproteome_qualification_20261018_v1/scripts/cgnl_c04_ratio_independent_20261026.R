library(jsonlite)
p<-'/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1'
w<-file.path(p,'C04_full_fraction_qualification/relative_peptide_pattern_20261026')
raw<-read.delim(file.path(p,'C04_full_fraction_qualification/reporter_normalization_20261026/CGNL1_RAW_REPORTERS.tsv'),check.names=FALSE,na.strings='',stringsAsFactors=FALSE)
map<-fromJSON(file.path(p,'C04_sources/PAK_BIOLOGICAL_SAMPLE_MAPPING.json'))$records
expected<-fromJSON(file.path(w,'RESULTS.json'));animals<-fromJSON(file.path(w,'ANIMAL_VALUES.json'));loo<-fromJSON(file.path(w,'LOO.json'))
v<-vapply(map$phospho_TMT_channel,function(c){a<-raw$intensity[raw$MS3 %in% c(26697,26749)&raw$file=='phospho_fraction_1.raw'&raw$channel==c];b<-raw$intensity[raw$MS3==19264&raw$file=='phospho_fraction_11.raw'&raw$channel==c];stopifnot(length(a)==2,length(b)==1,!anyNA(c(a,b)));log2(sum(a)/b)},0)
err<-abs(v-animals$log2_ratio[match(map$sample,animals$sample)])
a<-v[map$group=='inflamed'];b<-v[map$group=='control'];delta<-mean(a)-mean(b)
left<-vapply(seq_along(v),function(i)mean(v[map$group=='inflamed'&seq_along(v)!=i])-mean(v[map$group=='control'&seq_along(v)!=i]),0)
err<-c(err,abs(delta-expected$mean_log2_ratio_difference),abs(2^delta-expected$relative_ratio_fold_change),abs(median(a)-median(b)-expected$median_log2_ratio_difference),abs(left-loo$difference[match(map$sample,loo$left_out)]))
stopifnot(max(err)<1e-10,sum(a<median(b))==expected$inflamed_below_control_median)
q<-list(pass=TRUE,max_absolute_error=max(err),animals=8,posthoc=TRUE,not_independent_validation=TRUE,P_values_computed=FALSE)
write_json(q,file.path(w,'INDEPENDENT_R_QA.json'),auto_unbox=TRUE,pretty=TRUE,digits=17);cat(toJSON(q,auto_unbox=TRUE,digits=17),'\n')
