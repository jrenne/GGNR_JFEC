# Check archived simulation design, joint covariance, endpoints and independent replays.
source("code/inference/parameterization.R")
source("code/inference/premium_functions.R")
source("code/inference/uncertainty_helpers.R")
source("code/validation/cache_provenance.R")
z<-readRDS(file.path(Sys.getenv("GGNR_UNCERTAINTY_DIR","outputs/diagnostics/uncertainty"),"nonlinear_uncertainty.rds"))
ggnr_check_provenance(z,"Uncertainty simulation")
stopifnot(identical(z$generator_md5,unname(tools::md5sum(c("code/inference/nonlinear_uncertainty.R","code/inference/uncertainty_helpers.R")))))
pr<-make_inference_problem();base<-pr$evaluate_fit(pr$theta)
stopifnot(isTRUE(all.equal(z$theta,pr$theta,tolerance=1e-12)),max(abs(z$point-ggnr_uncertainty_series(base)))<1e-8,
 z$diagnostics$accepted_draws==sum(vapply(z$batches,`[[`,integer(1),"accepted")),
 z$diagnostics$proposals==sum(vapply(z$batches,`[[`,integer(1),"proposals")))
for(b in z$batches)stopifnot(b$accepted+sum(b$rejection_reasons)==b$proposals)
e<-eigen((z$parameter_covariance+t(z$parameter_covariance))/2,symmetric=TRUE);stopifnot(min(e$values)>0)
root<-sweep(e$vectors,2,sqrt(e$values),"*")
errors<-vapply(z$reference_draws,function(r){a<-ggnr_uncertainty_draw(r$input,pr,z$theta,root);stopifnot(a$ok);max(abs(a$paths-r$paths))},numeric(1))
stopifnot(max(errors)<1e-7)
for(item in c("state_bands","premium_bands")){
 bands<-z[[item]]
 for(n in sub("_estimate$","",names(bands)[grepl("_estimate$",names(bands))])){
  q<-as.matrix(bands[paste0(n,"_",c("p05","p16","median","p84","p95"))]);stopifnot(all(is.finite(q)),all(q[,1]<=q[,2]),all(q[,2]<=q[,3]),all(q[,3]<=q[,4]),all(q[,4]<=q[,5]),max(abs(bands[[paste0(n,"_estimate")]]-z$point[,n]))<1e-12)
 }
}
message("Nonlinear bands: joint covariance, admissibility counts, quantiles and three independent draw replays agree.")
