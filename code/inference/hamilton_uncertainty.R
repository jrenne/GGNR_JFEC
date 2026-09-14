# Hamilton-style state and risk-premium uncertainty for the preferred model.
# Each replication draws working parameters from the HAC approximation,
# reruns the estimation filter, draws from each filtered state distribution, and reprices
# bonds. Pointwise bands therefore use the same information set as the paper.

source("code/model/iekf.R")
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
Rcpp::sourceCpp("code/model/pricing.cpp")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")
Rcpp::sourceCpp("code/inference/hamilton_state_sampler.cpp")

estimate_file <- Sys.getenv(
  "GGNR_START_FILE",
  "estimates/baseline.rds"
)
inference_file <- Sys.getenv(
  "GGNR_HAC_FILE",
  "outputs/diagnostics/inference/hac_sandwich_covariance.rds"
)
output_directory <- Sys.getenv(
  "GGNR_HAMILTON_DIR",
  "outputs/diagnostics/uncertainty"
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

p_template <- modifyList(reduced_hfi_no_output_defaults(), estimate$parameters)
p_template$include_core_cpi <- FALSE
p_template$sigma_inflation_survey_long <- 0.10 / 1200
p_template$sigma_tbill_survey_long <- 0.10 / 1200
p_template$sigma_real_yield <- 0.10 / 1200
p_template$sigma_real_yield_liquidity <- 0.20 / 1200
p_template$liquidity_spread_covid <- 0

data <- load_true_release_data()
short_rate <- load_observed_real_short_rate_inputs(data$dates)
hfi_none <- list(
  standardized = rep(NaN, length(data$dates)),
  observed = rep(FALSE, length(data$dates)), source = "none",
  shock = NA_character_, event_count = 0L, month_count = 0L
)

contractive_from_raw <- function(raw) {
  decomposition <- svd(raw)
  singular <- decomposition$d / sqrt(1 + decomposition$d^2)
  decomposition$u %*% diag(singular, nrow = length(singular)) %*%
    t(decomposition$v)
}
free_measurement <- c(
  "sigma_inf", "sigma_ptr", "sigma_inflation_survey",
  "sigma_tbill_survey", "sigma_nominal_yield"
)
from_working <- function(theta) {
  p <- p_template
  p$rho_r_star <- plogis(theta["rho_r_star"])
  p$sigma_r_star <- exp(theta["sigma_r_star"]) / 1200
  p$rho_pi_star <- plogis(theta["rho_pi_star"])
  p$sigma_pi_star <- exp(theta["sigma_pi_star"]) / 1200
  raw <- matrix(theta[paste0("cyclical_raw_", seq_len(4))], 2, 2)
  p$cyclical_var <- contractive_from_raw(raw)
  dimnames(p$cyclical_var) <- list(c("m", "pi_gap"), c("m", "pi_gap"))
  p$sigma_pi_gap <- exp(theta["sigma_pi_gap"]) / 1200
  p$impact_pi_m_ratio <- unname(theta["impact_pi_m_ratio"])
  p$liquidity_spread_gfc <- exp(theta["liquidity_spread_gfc"]) / 1200
  for (name in free_measurement) {
    p[[name]] <- exp(theta[paste0("measurement_", name)])
  }
  p$lambda_0[reduced_hfi_no_output_shocks] <-
    theta[paste0("lambda0_", reduced_hfi_no_output_shocks)]
  p$lambda_w[reduced_hfi_no_output_shocks] <-
    theta[paste0("lambdaw_", reduced_hfi_no_output_shocks)]
  p$lambda_r_star_scale <- unname(theta["lambda_r_star_scale"])
  p$lambda_pi_star_scale <- unname(theta["lambda_pi_star_scale"])
  p$rho_w <- plogis(theta["rho_w"])
  p$sigma_w <- sqrt(1 - p$rho_w^2)
  p$sigma_wpi <- unname(theta["sigma_wpi"])
  p
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

inflation_compensation <- function(coefs, states, maturities) {
  vapply(maturities, function(maturity) {
    loading <- rowSums(coefs$B_X_pi[, seq_len(maturity), drop = FALSE]) /
      maturity
    intercept <- sum(coefs$A_X_exp_pi[seq_len(maturity)]) / maturity +
      0.5 * sum(coefs$Conv_pi[seq_len(maturity)]) / maturity^2
    intercept + as.numeric(states %*% loading)
  }, numeric(nrow(states)))
}

calculate_premia <- function(fit, states) {
  fitted <- function(coefs) y_fitting_r(
    t(states), coefs$A_X_for, coefs$B_X_for, coefs$A_X_exp,
    fit$objects$pars$r_lb, fit$objects$s_n, coefs$A_X_for_pi,
    coefs$B_X_pi, fit$objects$Sigma2_X, coefs$B_X_cum,
    coefs$B_X_cum_pi
  )
  physical <- fitted(fit$objects$coefs_p)
  risk_neutral <- fitted(fit$objects$coefs_q)
  maturities <- c(60L, 120L)
  nominal <- 1200 * (risk_neutral$yfit_all_n[, maturities] -
                       physical$yfit_all_n[, maturities])
  real <- 1200 * (risk_neutral$yfit_all_r[, maturities] -
                    physical$yfit_all_r[, maturities])
  inflation <- 1200 * (
    inflation_compensation(fit$objects$coefs_q, states, maturities) -
      inflation_compensation(fit$objects$coefs_p, states, maturities)
  )
  cbind(
    nominal_5y = nominal[, 1], nominal_10y = nominal[, 2],
    real_5y = real[, 1], real_10y = real[, 2],
    inflation_5y = inflation[, 1], inflation_10y = inflation[, 2]
  )
}

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
selected_states <- c("r_star", "pi_star")
baseline_states <-
  1200 * baseline_fit$x_upd[, selected_states, drop = FALSE]
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
        numeric(2)
      ))
      list(
        ok = TRUE, attempt = attempt,
        state_path = 1200 * state_path[, selected_states, drop = FALSE],
        state_mean = 1200 * fit$x_upd[, selected_states, drop = FALSE],
        state_filter_variance = 1200^2 * conditional_variance,
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
  state_paths, c("r_star", "pi_star"), baseline_states
)
premium_bands <- array_bands(
  premium_paths,
  c("nominal_5y", "nominal_10y", "real_5y", "real_10y",
    "inflation_5y", "inflation_10y"), baseline_premia
)

state_variance_decomposition <- data.frame(date = data$dates)
for (series in seq_len(2)) {
  filter_component <- rowMeans(state_filter_variances[, series, ])
  parameter_component <- apply(state_means[, series, ], 1, var)
  prefix <- c("r_star", "pi_star")[series]
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
  state_bands = state_bands, premium_bands = premium_bands,
  state_variance_decomposition = state_variance_decomposition,
  diagnostics = diagnostics, parameter_seeds = draw_seeds,
  parameter_covariance = working_covariance,
  parameter_covariance_eigenvalues = parameter_decomposition$eigenvalues,
  baseline_states = baseline_states, baseline_premia = baseline_premia,
  state_information_set = "filtered",
  parameter_distribution = "asymptotic normal using HAC sandwich covariance",
  exact_moment_targets_treated_as_fixed = TRUE
), file.path(output_directory, "hamilton_uncertainty.rds"))

cat(sprintf(
  "Hamilton simulation: %d draws in %.1f seconds on %d core(s); acceptance %.1f%%\n",
  draw_count, elapsed, core_count, 100 * diagnostics$acceptance_rate
))
print(diagnostics)
