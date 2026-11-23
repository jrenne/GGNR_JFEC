source("code/validation/cache_provenance.R")
# Recompute local empirical curvature and scores; four cores by default.
source("code/inference/parameterization.R")
source("code/inference/sandwich_helpers.R")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")

problem <- make_inference_problem()
theta <- problem$theta
base <- problem$evaluate_fit(theta, record = TRUE)
branches <- base$nominal_branches
stopifnot(abs(base$fval - problem$estimate$final_loglik) < 1e-7,
          length(branches) > 0L)
evaluate_local <- function(value) problem$evaluate(value, branches)
stopifnot(max(abs(evaluate_local(theta) - base$fval_t)) < 1e-10)
started <- proc.time()["elapsed"]
step_scale <- as.numeric(Sys.getenv("GGNR_HESSIAN_STEP", "1e-4"))
cores <- as.integer(Sys.getenv("GGNR_HESSIAN_CORES", "4"))
stopifnot(is.finite(step_scale), step_scale > 0, is.finite(cores), cores >= 1)
message("Computing the local Hessian and monthly scores at the R2 estimate.")
derivatives <- numerical_local_derivatives(theta, evaluate_local, step_scale, cores)
J <- transformation_jacobian(theta, problem$natural)
result <- assemble_hessian_sandwich(theta, derivatives$H, derivatives$scores,
  J, problem$natural(theta), as.integer(Sys.getenv("GGNR_HAC_LAGS", "12")))
result$derivative_step <- derivatives$step
result$derivative_step_scale <- step_scale
result$branches <- branches
result$estimated_loglik <- base$fval
result$elapsed_seconds <- unname(proc.time()["elapsed"] - started)
result$estimate_md5 <- unname(tools::md5sum(problem$estimate_file))
result$exact_targets <- problem$targets
output_directory <- Sys.getenv("GGNR_INFERENCE_DIR", "outputs/diagnostics/inference")
result$provenance <- ggnr_cache_provenance()
write_hessian_inference(result, output_directory)
message(sprintf("Local Hessian: %d positive, %d negative eigenvalues; %.1f seconds.",
  sum(eigen(result$H, symmetric = TRUE)$values > 0),
  sum(eigen(result$H, symmetric = TRUE)$values < 0), result$elapsed_seconds))
