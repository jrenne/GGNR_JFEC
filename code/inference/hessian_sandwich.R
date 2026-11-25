# Full Hessian--HAC covariance, with a nonsingular linear numerical rescaling.
source("code/validation/cache_provenance.R")
source("code/inference/parameterization.R")
source("code/inference/sandwich_helpers.R")
problem <- make_inference_problem(); theta <- problem$theta
base <- problem$evaluate_fit(theta, record=TRUE)
stopifnot(abs(base$fval-problem$estimate$final_loglik)<1e-7)
B <- problem$estimate$differentiation_basis
if(is.null(B))B<-diag(length(theta))
B<-B[names(theta),names(theta),drop=FALSE]
zero<-setNames(rep(0,length(theta)),names(theta))
expand<-function(z)setNames(theta+as.vector(B%*%z),names(theta))
step<-as.numeric(Sys.getenv("GGNR_HESSIAN_STEP","0.001"))
cores<-as.integer(Sys.getenv("GGNR_HESSIAN_CORES","4"))
started<-proc.time()["elapsed"]
calculate<-function(h){
 d<-numerical_local_derivatives(zero,function(z)problem$evaluate(expand(z),base$nominal_branches),h,cores)
 if(min(eigen(d$H,symmetric=TRUE,only.values=TRUE)$values)<=0)
   stop("The local negative Hessian is not positive definite. Check convergence and differentiation; no eigenvalues are modified.")
 Bi<-solve(B);H<-t(Bi)%*%d$H%*%Bi;S<-d$scores%*%Bi
 dimnames(H)<-list(names(theta),names(theta));colnames(S)<-names(theta)
 Jz<-transformation_jacobian(zero,function(z)problem$natural(expand(z)),step_scale=1e-4)
 J<-Jz%*%Bi;colnames(J)<-names(theta)
 r<-assemble_hessian_sandwich(theta,H,S,J,problem$natural(theta),as.integer(Sys.getenv("GGNR_HAC_LAGS","12")))
 r$balanced_H<-d$H;r$linear_basis<-B;r$derivative_step_scale<-h
 r$derivative_step<-d$step;r$branches<-base$nominal_branches
 r$estimated_loglik<-base$fval;r$exact_targets<-problem$targets
 r$estimate_md5<-unname(tools::md5sum(problem$estimate_file));r$provenance<-ggnr_cache_provenance()
 r$elapsed_seconds<-unname(proc.time()["elapsed"]-started);r
}
result<-calculate(step);half<-calculate(step/2)
directory<-Sys.getenv("GGNR_INFERENCE_DIR","outputs/diagnostics/inference")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
write_hessian_inference(result,directory)
saveRDS(half,file.path(directory,"halfstep_check.rds"))
free<-result$summary$hac_sandwich_se>0
change<-abs(half$summary$hac_sandwich_se[free]/result$summary$hac_sandwich_se[free]-1)
write.csv(data.frame(median_relative_change=median(change),max_relative_change=max(change)),file.path(directory,"fresh_step_check.csv"),row.names=FALSE)
message(sprintf("Positive-definite %d-parameter Hessian; median/max SE step changes %.3g / %.3g; %.1fs.",length(theta),median(change),max(change),proc.time()["elapsed"]-started))
