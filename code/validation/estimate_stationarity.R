# Independent finite-difference checks of the original, unfixed IEKF criterion.
source("code/inference/parameterization.R")
p <- make_inference_problem(); u <- p$theta
f <- p$evaluate_fit(u,record=TRUE)
checks <- lapply(c(1e-6,2e-7),function(h){
 parts <- parallel::mclapply(seq_along(u),function(j){
  hh<-h*max(1,abs(u[j]));aa<-bb<-u;aa[j]<-aa[j]+hh;bb[j]<-bb[j]-hh
  a<-p$evaluate_fit(aa,record=TRUE);b<-p$evaluate_fit(bb,record=TRUE)
  list(score=(a$fval_t-b$fval_t)/(2*hh),switches=sum(a$nominal_branches!=b$nominal_branches))
 },mc.cores=as.integer(Sys.getenv("GGNR_CHECK_CORES","4")))
 S<-do.call(cbind,lapply(parts,`[[`,"score"))
 data.frame(step=h,max_normalized_score=max(abs(colSums(S))/sqrt(colSums(S*S))),formula_switches=sum(vapply(parts,`[[`,numeric(1),"switches")))
})
checks<-do.call(rbind,checks)
if(max(checks$max_normalized_score)>1e-4||any(checks$formula_switches>0))
 stop("Local stationarity has not been verified at the selected point. Inspect the saved candidate before using inference.")
write.csv(checks,"outputs/diagnostics/inference/stationarity_checks.csv",row.names=FALSE)
message("Original-criterion local stationarity verified at two steps, with no lower-bound formula switches.")
