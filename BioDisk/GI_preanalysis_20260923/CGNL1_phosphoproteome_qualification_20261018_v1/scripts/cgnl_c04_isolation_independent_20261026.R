library(jsonlite)
w <- '/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1/C04_full_fraction_qualification/joint_spectral_qualification_20261026'
r <- fromJSON(file.path(w,'ISOLATION_RESULTS.json'),simplifyVector=FALSE)$records
errors <- c();out <- list()
for(k in seq_along(r)) {
 x<-r[[k]];mz<-vapply(x$window_centroids,function(z)z$mz,0);it<-vapply(x$window_centroids,function(z)z$intensity,0)
 stopifnot(all(mz>=x$isolation_center-x$lower_offset),all(mz<=x$isolation_center+x$upper_offset),all(it>=0))
 selected<-abs(mz-x$selected_mz)<=x$selected_mz*1e-5
 center<-abs(mz-x$isolation_center)<=x$isolation_center*1e-5
 grid<-x$selected_mz+(-2:6)*1.0033548378/x$charge
 lattice<-vapply(mz,function(m)any(abs(m-grid)<=grid*1e-5),TRUE)
 stopifnot(identical(selected,vapply(x$window_centroids,function(z)z$selected_match,TRUE)),identical(center,vapply(x$window_centroids,function(z)z$center_match,TRUE)),identical(lattice,vapply(x$window_centroids,function(z)z$lattice_match,TRUE)))
 vals<-c(sum(it[selected]),sum(it[center]),sum(it[lattice]))/sum(it)
 ref<-c(x$selected_mz_peak_fraction,x$isolation_center_peak_fraction,x$charge_spaced_lattice_fraction)
 errors<-c(errors,abs(vals-ref),abs(sum(it)-x$total_intensity))
 out[[k]]<-data.frame(file=x$file,HCD=x$HCD,selected=vals[1],center=vals[2],lattice=vals[3])
}
stopifnot(max(errors)<1e-9,length(r)==4)
write.csv(do.call(rbind,out),file.path(w,'INDEPENDENT_R_ISOLATION_QA.csv'),row.names=FALSE,fileEncoding='UTF-8')
write_json(list(pass=TRUE,records=4,max_absolute_error=max(errors),original_PIF_equivalence=FALSE),file.path(w,'INDEPENDENT_R_ISOLATION_QA.json'),auto_unbox=TRUE,pretty=TRUE)
cat('PASS 4 windows; max error',max(errors),'\n')
