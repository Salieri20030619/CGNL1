library(jsonlite)
w <- '/mnt/BioDisk/wangmingyuan/GI_preanalysis_20260923/CGNL1_phosphoproteome_qualification_20261018_v1/continuous_exploration_20261023_v1/C04_full_fraction_qualification/joint_spectral_qualification_20261026'
raw <- fromJSON(file.path(w,'TARGET_NOISE_ARRAYS.json'), simplifyVector=FALSE)
res <- fromJSON(file.path(w,'NOISE_RESULTS.json'), simplifyVector=FALSE)$records
errors <- c(); checks <- list()
for (k in seq_along(raw)) {
 a <- raw[[k]]$arrays; r <- res[[k]]
 stopifnot(raw[[k]]$scan==r$scan, raw[[k]]$file==r$file)
 mz <- unlist(a$mz); intensity <- unlist(a$intensity)
 grid <- unlist(a$noise_mz); noise <- unlist(a$noise); baseline <- unlist(a$baseline)
 ratios <- c()
 for (p in r$reporters) {
  hit <- which(abs(mz-p$target_mz)<=p$target_mz*20e-6)
  stopifnot(length(hit)==p$matches)
  if (length(hit)!=1) {stopifnot(is.null(p$I_over_noise));next}
  n <- approx(grid, noise, xout=mz[hit], rule=1)$y
  b <- approx(grid, baseline, xout=mz[hit], rule=1)$y
  stopifnot(is.finite(n),n>0,is.finite(b))
  ratio <- intensity[hit]/n; adjusted <- (intensity[hit]-b)/n
  errors <- c(errors, abs(ratio-p$I_over_noise), abs(adjusted-p$baseline_adjusted_proxy))
  ratios <- c(ratios, ratio)
 }
 stopifnot(length(ratios)==r$available_reporter_ratios)
 errors <- c(errors,abs(sum(ratios)-r$sum_available_I_over_noise))
 checks[[k]] <- data.frame(file=r$file,scan=r$scan,available_reporters=length(ratios),sum_I_over_noise=sum(ratios))
}
stopifnot(max(errors)<1e-9)
write.csv(do.call(rbind,checks),file.path(w,'INDEPENDENT_R_NOISE_QA.csv'),row.names=FALSE,fileEncoding='UTF-8')
write_json(list(pass=TRUE,max_absolute_error=max(errors),records=length(raw),original_SNR_equivalence=FALSE),file.path(w,'INDEPENDENT_R_NOISE_QA.json'),auto_unbox=TRUE,pretty=TRUE)
cat('PASS',length(raw),'spectra; max error',format(max(errors),digits=17),'\n')
