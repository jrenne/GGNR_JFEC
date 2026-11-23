source("code/validation/cache_provenance.R")
# Archived nonlinear Gaussian-draw diagnostic; primary bands use delta_uncertainty.R.
# Each replication draws working parameters from the HAC approximation,
# reruns the estimation filter, draws from each filtered state distribution, and reprices
# bonds. The information set matches the paper, but these are not its primary bands.

source("code/model/iekf.R")
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")
Rcpp::sourceCpp("code/inference/hamilton_state_sampler.cpp")

estimate_file <- Sys.getenv(
  "GGNR_START_FILE",
  "estimates/baseline_r2.rds"
)
inference_file <- Sys.getenv(
  "GGNR_HAC_FILE",
  file.path(Sys.getenv("GGNR_INFERENCE_DIR", "outputs/diagnostics/inference"),
            "hac_sandwich_covariance.rds")
)
output_directory <- Sys.getenv(
  "GGNR_HAMILTON_DIR",
  "outputs/diagnostics/hamilton_reference"
)
draw_count <- as.integer(Sys.getenv("GGNR_HAMILTON_DRAWS", "1000"))
core_count <- as.integer(Sys.getenv(
  "GGNR_HAMILTON_CORES", as.character(min(4L, parallel::detectCores()))
))
random_seed <- as.integer(Sys.getenv("GGNR_HAMILTON_SEED", "1986"))
max_attempts <- as.integer(Sys.getenv("GGNR_HAMILTON_MAX_ATTEMPTS", "50"))
if (!is.finite(draw_count) || draw_count < 2L) stop("At least two draws are required.")
if (!is.finite(core_count) || core_count < 1L) core_count <- 1L
if (!is.finite(max_attempts) || max_attempts < 1L) max_attempts <- 50L

estimate <- readRDS(estimate_file)
inference <- readRDS(inference_file)
filter_method <- if (is.null(estimate$filter_method)) "EKF" else
  estimate$filter_method
iekf_trigger_probability <- if (is.null(estimate$iekf_trigger_probability))
  0 else estimate$iekf_trigger_probability
iekf_steps <- if (is.null(estimate$iekf_steps)) 2L else estimate$iekf_steps
theta_hat <- estimate$working_parameters
working_covariance <- inference$hac_working
if (is.null(theta_hat) || is.null(working_covariance)) {
  stop("Working estimates or their HAC covariance are unavailable.")
}
working_covariance <- working_covariance[names(theta_hat), names(theta_hat)]

source("code/inference/parameterization.R")
data <- load_paper_data()
short_rate <- load_observed_real_short_rate_inputs(data$dates)
hfi_none <- no_hfi_data(data$dates)
targets <- reduced_hfi_exact_targets(data, short_rate)
from_working <- make_parameterization(estimate$parameters, targets)$from_working
if (is.null(inference$theta) ||
    !isTRUE(all.equal(theta_hat, inference$theta, tolerance = 1e-12))) {
  stop("The parameter covariance belongs to another empirical estimate.")
}

positive_sqrt <- function(value, tolerance = 1e-12) {
  decomposition <- eigen((value + t(value)) / 2, symmetric = TRUE)
  threshold <- max(1, max(abs(decomposition$values))) * tolerance
  if (min(decomposition$values) < -100 * threshold) {
    stop("The HAC working covariance is materially indefinite.")
  }
  keep <- decomposition$values > threshold
  list(
    root = decomposition$vectors[, keep, drop = FALSE] %*%
      diag(sqrt(decomposition$values[keep]), nrow = sum(keep)),
    rank = sum(keep), eigenvalues = decomposition$values,
    threshold = threshold
  )
}
parameter_decomposition <- positive_sqrt(working_covariance)
parameter_root <- parameter_decomposition$root

source("code/inference/premium_functions.R")

baseline_fit <- run_kf_reduced_hfi_no_output(
  data, from_working(theta_hat), hfi_none, short_rate,
  filter_method = filter_method,
  iekf_trigger_probability = iekf_trigger_probability,
  iekf_steps = iekf_steps
)
if (!is.finite(baseline_fit$fval) ||
    abs(baseline_fit$fval - estimate$final_loglik) > 1e-6) {
  stop("The working-parameter mapping does not reproduce the saved estimate.")
}
selected_states <- c("r_star", "pi_star", "w")
# Trends are reported in annualized percentage points; w is dimensionless.
state_scales <- c(r_star = 1200, pi_star = 1200, w = 1)
baseline_states <- sweep(
  baseline_fit$x_upd[, selected_states, drop = FALSE], 2, state_scales, "*"
)
baseline_premia <- calculate_premia(baseline_fit, baseline_fit$x_upd)

evaluate_draw <- function(seed) {
  set.seed(seed)
  last_error <- "maximum attempts reached"
  for (attempt in seq_len(max_attempts)) {
    theta <- theta_hat + as.numeric(parameter_root %*% rnorm(ncol(parameter_root)))
    names(theta) <- names(theta_hat)
    result <- tryCatch({
      fit <- run_kf_reduced_hfi_no_output(
        data, from_working(theta), hfi_none, short_rate,
        filter_method = filter_method,
        iekf_trigger_probability = iekf_trigger_probability,
        iekf_steps = iekf_steps
      )
      state_path <- draw_filtered_states_cpp(
        fit$x_upd, fit$P_upd,
        matrix(rnorm(nrow(fit$x_upd) * ncol(fit$x_upd)),
               nrow(fit$x_upd), ncol(fit$x_upd))
      )
      premia <- calculate_premia(fit, state_path)
      if (any(!is.finite(c(state_path, premia)))) {
        stop("non-finite simulated output")
      }
      conditional_variance <- t(vapply(
        seq_len(nrow(state_path)),
        function(index) diag(fit$P_upd[, , index])[selected_states],
        numeric(length(selected_states))
      ))
      list(
        ok = TRUE, attempt = attempt,
        state_path = sweep(state_path[, selected_states, drop = FALSE],
                           2, state_scales, "*"),
        state_mean = sweep(fit$x_upd[, selected_states, drop = FALSE],
                           2, state_scales, "*"),
        state_filter_variance = sweep(conditional_variance,
                                     2, state_scales^2, "*"),
        premia = premia
      )
    }, error = function(error) {
      last_error <<- conditionMessage(error)
      NULL
    })
    if (!is.null(result)) return(result)
  }
  list(ok = FALSE, attempt = max_attempts, error = last_error)
}

set.seed(random_seed)
draw_seeds <- sample.int(.Machine$integer.max, draw_count)
start_time <- proc.time()["elapsed"]
if (.Platform$OS.type == "unix" && core_count > 1L) {
  simulations <- parallel::mclapply(
    draw_seeds, evaluate_draw, mc.cores = core_count, mc.preschedule = TRUE
  )
} else {
  simulations <- lapply(draw_seeds, evaluate_draw)
}
elapsed <- unname(proc.time()["elapsed"] - start_time)
successful <- vapply(simulations, function(value) isTRUE(value$ok), logical(1))
if (!all(successful)) {
  failures <- vapply(simulations[!successful], `[[`, character(1), "error")
  stop(sprintf("%d simulations failed; first error: %s",
               sum(!successful), failures[1]))
}

state_paths <- simplify2array(lapply(simulations, `[[`, "state_path"))
state_means <- simplify2array(lapply(simulations, `[[`, "state_mean"))
state_filter_variances <- simplify2array(
  lapply(simulations, `[[`, "state_filter_variance")
)
premium_paths <- simplify2array(lapply(simulations, `[[`, "premia"))
quantile_probabilities <- c(0.025, 0.16, 0.50, 0.84, 0.975)
quantile_names <- c("p025", "p16", "median", "p84", "p975")

array_bands <- function(values, series_names, point_estimates) {
  output <- data.frame(date = data$dates)
  for (series in seq_along(series_names)) {
    quantiles <- t(apply(values[, series, , drop = FALSE], 1, quantile,
                         probs = quantile_probabilities, na.rm = TRUE))
    colnames(quantiles) <- paste(series_names[series], quantile_names, sep = "_")
    estimate_column <- point_estimates[, series]
    output[[paste0(series_names[series], "_estimate")]] <- estimate_column
    output <- cbind(output, quantiles)
  }
  output
}
state_bands <- array_bands(
  state_paths, selected_states, baseline_states
)
premium_bands <- array_bands(
  premium_paths,
  c("nominal_5y", "nominal_10y", "real_5y", "real_10y",
    "inflation_5y", "inflation_10y"), baseline_premia
)

state_variance_decomposition <- data.frame(date = data$dates)
for (series in seq_along(selected_states)) {
  filter_component <- rowMeans(state_filter_variances[, series, ])
  parameter_component <- apply(state_means[, series, ], 1, var)
  prefix <- selected_states[series]
  state_variance_decomposition[[paste0(prefix, "_filter_variance")]] <-
    filter_component
  state_variance_decomposition[[paste0(prefix, "_parameter_variance")]] <-
    parameter_component
  state_variance_decomposition[[paste0(prefix, "_total_se")]] <-
    sqrt(filter_component + parameter_component)
}

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(state_bands, file.path(output_directory, "state_bands.csv"),
          row.names = FALSE)
write.csv(premium_bands, file.path(output_directory, "premium_bands.csv"),
          row.names = FALSE)
write.csv(state_variance_decomposition,
          file.path(output_directory, "state_variance_decomposition.csv"),
          row.names = FALSE)
attempts <- vapply(simulations, `[[`, integer(1), "attempt")
diagnostics <- data.frame(
  requested_draws = draw_count, successful_draws = sum(successful),
  total_attempts = sum(attempts), acceptance_rate = draw_count / sum(attempts),
  cores = core_count, elapsed_seconds = elapsed, seed = random_seed,
  parameter_covariance_rank = parameter_decomposition$rank,
  parameter_count = length(theta_hat),
  conditional_on_exact_moment_targets = TRUE,
  reconstructed_loglik = baseline_fit$fval,
  saved_loglik = estimate$final_loglik,
  filter_method = filter_method,
  iekf_trigger_probability = iekf_trigger_probability,
  iekf_steps = iekf_steps
)
write.csv(diagnostics, file.path(output_directory, "simulation_diagnostics.csv"),
          row.names = FALSE)
saveRDS(list(
  provenance = ggnr_cache_provenance(),
  state_bands = state_bands, premium_bands = premium_bands,
  state_variance_decomposition = state_variance_decomposition,
  diagnostics = diagnostics, parameter_seeds = draw_seeds,
  parameter_covariance = working_covariance,
  parameter_covariance_eigenvalues = parameter_decomposition$eigenvalues,
  baseline_states = baseline_states, baseline_premia = baseline_premia,
  theta = theta_hat,
  state_information_set = "filtered",
  parameter_distribution = "local Gaussian approximation using Hessian-HAC covariance",
  exact_moment_targets_treated_as_fixed = TRUE
), file.path(output_directory, "hamilton_uncertainty.rds"))

cat(sprintf(
  "Hamilton simulation: %d draws in %.1f seconds on %d core(s); acceptance %.1f%%\n",
  draw_count, elapsed, core_count, 100 * diagnostics$acceptance_rate
))
print(diagnostics)
