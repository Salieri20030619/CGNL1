library(jsonlite)
p <- '/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1'
n<-file.path(p,'C04_full_fraction_qualification/reporter_normalization_20261026');w<-file.path(p,'C04_full_fraction_qualification/disease_descriptive_20261026')
read<-function(path)read.delim(path,check.names=FALSE,stringsAsFactors=FALSE,na.strings='')
spec<-fromJSON(file.path(w,'SPEC_BEFORE_GROUP_COMPARISON.json'),simplifyVector=FALSE)
map<-fromJSON(file.path(p,'C04_sources/PAK_BIOLOGICAL_SAMPLE_MAPPING.json'))$records
raw<-read(file.path(n,'CGNL1_RAW_REPORTERS.tsv'));fac<-read(file.path(n,'NORMALIZATION_FACTORS.tsv'));observed<-read(file.path(w,'ANIMAL_VALUES.tsv'));loo_source<-read(file.path(w,'LOO_RESULTS.tsv'))
results<-fromJSON(file.path(w,'RESULTS.json'),simplifyVector=FALSE)$results
stopifnot(nrow(map)==8,!anyDuplicated(map$sample),sum(map$group=='control')==5,sum(map$group=='inflamed')==3)
errors<-c();comparisons<-list();count<-0
for(ep in spec$endpoints)for(model in names(spec$models)){
 vals<-vapply(seq_len(nrow(map)),function(i){
  z<-raw[raw$file==ep$file & raw$MS3 %in% unlist(ep$scans) & raw$channel==map$phospho_TMT_channel[i],]
  stopifnot(nrow(z)==length(ep$scans));v<-if(anyNA(z$intensity))NA_real_ else log2(sum(z$intensity))-fac$log2_offset[fac$model==model & fac$channel==map$phospho_TMT_channel[i]]
  row<-observed[observed$endpoint==ep$id & observed$model==model & observed$sample==map$sample[i],];stopifnot(nrow(row)==1,is.na(v)==is.na(row$log2_value))
  if(!is.na(v))errors<<-c(errors,abs(v-row$log2_value));v
 },0)
 a<-vals[map$group=='inflamed' & !is.na(vals)];b<-vals[map$group=='control' & !is.na(vals)];r<-Filter(function(x)x$endpoint==ep$id && x$model==model,results)[[1]]
 ok<-length(a)>=2 && length(b)>=2;stopifnot(r$evaluable==ok,r$n_control==length(b),r$n_inflamed==length(a));count<-count+1
 if(ok){
  delta<-mean(a)-mean(b);med<-median(a)-median(b);left<-vapply(seq_len(nrow(map)),function(i)mean(vals[map$group=='inflamed' & seq_along(vals)!=i])-mean(vals[map$group=='control' & seq_along(vals)!=i]),0)
  errors<-c(errors,abs(delta-r$mean_log2_difference),abs(2^delta-r$geometric_mean_fold_change),abs(med-r$median_log2_difference),abs(min(left)-r$LOO_min),abs(max(left)-r$LOO_max))
  stopifnot(sum(a>median(b))==r$inflamed_above_control_median,sum(a<median(b))==r$inflamed_below_control_median)
  for(i in seq_along(left)){z<-loo_source[loo_source$endpoint==ep$id & loo_source$model==model & loo_source$left_out_sample==map$sample[i],];stopifnot(nrow(z)==1);errors<-c(errors,abs(left[i]-z$mean_log2_difference))}
  comparisons[[length(comparisons)+1]]<-data.frame(endpoint=ep$id,model=model,mean_log2_difference=delta,fold_change=2^delta,LOO_min=min(left),LOO_max=max(left))
 } else stopifnot(is.null(r$mean_log2_difference))
}
stopifnot(count==9,max(errors)<1e-10,nrow(observed)==72,nrow(loo_source)==48)
write.csv(do.call(rbind,comparisons),file.path(w,'INDEPENDENT_R_RESULTS.csv'),row.names=FALSE,fileEncoding='UTF-8')
q<-list(pass=TRUE,endpoint_models=9,animal_values=72,LOO_values=48,max_absolute_error=max(errors),P_values_computed=FALSE)
write_json(q,file.path(w,'INDEPENDENT_R_QA.json'),auto_unbox=TRUE,pretty=TRUE,digits=17);cat(toJSON(q,auto_unbox=TRUE,digits=17),'\n')
