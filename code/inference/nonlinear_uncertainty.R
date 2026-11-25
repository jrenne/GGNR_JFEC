# Pointwise bands from joint parameter draws, refiltering and nonlinear repricing.
source("code/inference/parameterization.R")
source("code/inference/premium_functions.R")
source("code/inference/uncertainty_helpers.R")
source("code/validation/cache_provenance.R")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")
problem <- make_inference_problem()
theta <- problem$theta
inference <- readRDS(file.path(Sys.getenv("GGNR_INFERENCE_DIR","outputs/diagnostics/inference"),
                              "hac_sandwich_covariance.rds"))
ggnr_check_provenance(inference,"Parameter covariance")
stopifnot(isTRUE(all.equal(theta,inference$theta,tolerance=1e-12)))
V <- inference$hac_working[names(theta),names(theta)]
e <- eigen((V+t(V))/2,symmetric=TRUE)
stopifnot(min(e$values)>0)
root <- sweep(e$vectors,2,sqrt(e$values),"*")
base <- problem$evaluate_fit(theta)
point <- ggnr_uncertainty_series(base)
number <- as.integer(Sys.getenv("GGNR_UNCERTAINTY_DRAWS","5000"))
cores <- as.integer(Sys.getenv("GGNR_UNCERTAINTY_CORES","4"))
first_seed <- as.integer(Sys.getenv("GGNR_UNCERTAINTY_SEED","20261126"))
stopifnot(is.finite(number),number>=100L,number%%5L==0L,is.finite(cores),cores>=1)
target <- number%/%5L
output_directory <- Sys.getenv("GGNR_UNCERTAINTY_DIR","outputs/diagnostics/uncertainty")
dir.create(output_directory,recursive=TRUE,showWarnings=FALSE)
started <- proc.time()[["elapsed"]]
all_paths <- list(); batch_records <- list(); references <- list(); batches <- list()
for (batch in 1:5) {
  seed <- first_seed+batch-1L; set.seed(seed)
  inputs <- lapply(seq_len(max(1500L,2L*target)),function(i)
    list(noise=rnorm(length(theta)),state_seed=sample.int(.Machine$integer.max,1)))
  results <- list(); used <- 0L; accepted <- 0L
  while (accepted<target && used<length(inputs)) {
    indices <- seq.int(used+1L,min(used+250L,length(inputs)))
    worker <- function(i) ggnr_uncertainty_draw(inputs[[i]],problem,theta,root)
    draw <- if (.Platform$OS.type=="unix" && cores>1L)
      parallel::mclapply(indices,worker,mc.cores=cores) else lapply(indices,worker)
    if (any(vapply(draw,inherits,logical(1),"try-error"))) stop("An uncertainty worker failed.")
    results <- c(results,draw);used <- length(results)
    accepted <- sum(vapply(results,function(x)isTRUE(x$ok),logical(1)))
  }
  if (accepted<target) stop("Insufficient admissible parameter draws; inspect this specification.")
  good <- vapply(results,function(x)isTRUE(x$ok),logical(1))
  ending <- which(good)[target]; results<-results[seq_len(ending)];good<-good[seq_len(ending)]
  accepted_results <- results[good]
  paths <- simplify2array(lapply(accepted_results,`[[`,"paths"))
  all_paths[[batch]] <- paths
  q <- ggnr_uncertainty_quantiles(paths,point,problem$data$dates,seq_len(ncol(point)))
  batches[[batch]] <- data.frame(series=colnames(point),
    mean_width90=vapply(colnames(point),function(n)
      mean(q[[paste0(n,"_p95")]]-q[[paste0(n,"_p05")]]),numeric(1)))
  reasons <- table(vapply(results[!good],`[[`,character(1),"reason"))
  batch_records[[batch]] <- list(seed=seed,accepted=target,proposals=ending,rejection_reasons=reasons)
  if (batch<=3) {
    i <- which(good)[1]
    references[[batch]] <- list(input=inputs[[i]],paths=results[[i]]$paths)
  }
}
paths <- array(unlist(all_paths),c(nrow(point),ncol(point),number))
state_bands <- ggnr_uncertainty_quantiles(paths,point,problem$data$dates,1:7)
premium_bands <- ggnr_uncertainty_quantiles(paths,point,problem$data$dates,8:13)
summary <- data.frame(series=colnames(point),
  mean_sd=apply(apply(paths,c(1,2),sd),2,mean),
  mean_width90=vapply(colnames(point),function(n) {
    z<-if(paste0(n,"_estimate")%in%names(state_bands))state_bands else premium_bands
    mean(z[[paste0(n,"_p95")]]-z[[paste0(n,"_p05")]])},numeric(1)),
  relative_width_MC_SE=vapply(seq_len(ncol(point)),function(j) {
    x<-vapply(batches,function(b)b$mean_width90[j],numeric(1));sd(x)/sqrt(5)/mean(x)},numeric(1)))
proposals <- sum(vapply(batch_records,`[[`,integer(1),"proposals"))
diagnostics <- data.frame(method="nonlinear_covariance_propagation",parameter_count=length(theta),
 sample_months=nrow(point),interval_levels="68%,90%",accepted_draws=number,proposals=proposals,
 admissibility_rate=number/proposals,cores=cores,
 elapsed_seconds=unname(proc.time()[["elapsed"]]-started),
 reconstructed_loglik=base$fval,conditional_on_exact_moment_targets=TRUE,
 conditional_on_admissibility=TRUE,finite_sample_coverage_established=FALSE)
result<-list(theta=theta,parameter_covariance=V,point=point,state_bands=state_bands,
 premium_bands=premium_bands,summary=summary,diagnostics=diagnostics,batches=batch_records,
 reference_draws=references,provenance=ggnr_cache_provenance(),
 generator_md5=unname(tools::md5sum(c("code/inference/nonlinear_uncertainty.R","code/inference/uncertainty_helpers.R"))))
write.csv(state_bands,file.path(output_directory,"state_bands.csv"),row.names=FALSE)
write.csv(premium_bands,file.path(output_directory,"premium_bands.csv"),row.names=FALSE)
write.csv(summary,file.path(output_directory,"uncertainty_summary.csv"),row.names=FALSE)
write.csv(diagnostics,file.path(output_directory,"simulation_diagnostics.csv"),row.names=FALSE)
saveRDS(result,file.path(output_directory,"nonlinear_uncertainty.rds"))
message(sprintf("Nonlinear uncertainty: %d admissible draws / %d proposals; %.1f seconds.",number,proposals,diagnostics$elapsed_seconds))
