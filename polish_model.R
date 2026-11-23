# Optional original-criterion polishing of all 32 free coordinates.
# The archived R1/R2 estimates are never overwritten.
source("code/inference/parameterization.R")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")
problem <- make_inference_problem()
output_file <- Sys.getenv("GGNR_POLISH_OUTPUT", "estimates/polished.rds")
protected <- normalizePath(c("estimates/baseline.rds", "estimates/baseline_r2.rds"),
                          mustWork = FALSE)
if (normalizePath(output_file, mustWork = FALSE) %in% protected) {
  stop("Choose a new output file; archived estimates are read-only.")
}
seconds <- as.numeric(Sys.getenv("GGNR_POLISH_SECONDS", "180"))
max_evaluations <- as.integer(Sys.getenv("GGNR_POLISH_MAXEVAL", "1000"))
stopifnot(is.finite(seconds), seconds > 0, is.finite(max_evaluations), max_evaluations > 0)
inference <- readRDS("outputs/diagnostics/inference/hac_sandwich_covariance.rds")
theta <- problem$theta
stopifnot(identical(names(theta), colnames(inference$H)))
scales <- pmin(0.05, pmax(1e-5, 1 / sqrt(pmax(abs(diag(inference$H)), 1))))
best_theta <- theta
best_value <- sum(problem$evaluate(theta))
initial_value <- best_value
count <- 0L
started <- proc.time()["elapsed"]
trajectory <- list()
objective <- function(z) {
  count <<- count + 1L
  if (proc.time()["elapsed"] - started > seconds) stop("Polishing time limit reached.")
  candidate <- theta + scales * z
  if (any(!is.finite(candidate)) || any(abs(candidate) > 50)) return(1e12)
  value <- tryCatch(sum(problem$evaluate(candidate)), error = function(e) -Inf)
  if (!is.finite(value)) return(1e12)
  if (value > best_value) {
    best_theta <<- candidate; best_value <<- value
    trajectory[[length(trajectory) + 1L]] <<-
      data.frame(evaluation = count, loglik = value)
  }
  -value
}
optimizer <- tryCatch(nloptr::nloptr(x0 = rep(0, length(theta)), eval_f = objective,
  opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = max_evaluations,
              maxtime = seconds, xtol_rel = 1e-10, ftol_abs = 1e-8, print_level = 0)),
  error = function(e) list(message = conditionMessage(e)))
fit <- problem$evaluate_fit(best_theta)
stopifnot(fit$fval >= initial_value - 1e-7, abs(fit$fval - best_value) < 1e-7)
result <- problem$estimate
result$parameters <- fit$objects$pars
result$working_parameters <- best_theta
result$initial_loglik <- initial_value
result$final_loglik <- fit$fval
result$exact_targets <- fit$objects$targets
result$exact_moments <- fit$objects$exact_moments
result$optimizer <- optimizer
result$evaluations <- count
result$stage <- "time-limited 32-coordinate original-criterion polish"
result$input_file <- problem$estimate_file
result$timestamp <- Sys.time()
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, output_file)
if (length(trajectory)) write.csv(do.call(rbind, trajectory),
  paste0(output_file, ".trajectory.csv"), row.names = FALSE)
message(sprintf("Polish: %.9f -> %.9f; saved separately in %s.",
                initial_value, fit$fval, output_file))
