# First-order pointwise bands combining conditional filtering variance with
# full Hessian--HAC parameter uncertainty: P_t + D_t V_theta D_t'.
source("code/inference/parameterization.R")
source("code/inference/premium_functions.R")
source("code/validation/cache_provenance.R")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")
problem <- make_inference_problem()
theta <- problem$theta
inference_file <- file.path(Sys.getenv("GGNR_INFERENCE_DIR", "outputs/diagnostics/inference"),
                            "hac_sandwich_covariance.rds")
inference <- readRDS(inference_file)
ggnr_check_provenance(inference, "Parameter covariance")
stopifnot(isTRUE(all.equal(theta, inference$theta, tolerance = 1e-12)))
V <- inference$hac_working[names(theta), names(theta)]
base <- problem$evaluate_fit(theta, record = TRUE)
branches <- base$nominal_branches
dates <- load_paper_data()$dates
scales <- c(1200, 1200, 1200, 1200, 1200, 1)
economic_series <- function(fit) {
  state <- sweep(fit$x_upd, 2, scales, "*")
  cbind(state, shadow_rate = rowSums(state[, c("r_star", "pi_star", "m")]),
        calculate_premia(fit, fit$x_upd))
}
point <- economic_series(base)
cores <- as.integer(Sys.getenv("GGNR_DELTA_CORES", "4"))
step_scale <- as.numeric(Sys.getenv("GGNR_DELTA_STEP", "1e-5"))
stopifnot(is.finite(cores), cores >= 1, is.finite(step_scale), step_scale > 0)
output_directory <- Sys.getenv("GGNR_UNCERTAINTY_DIR", "outputs/diagnostics/uncertainty")
started <- proc.time()['elapsed']
TT <- nrow(point); K <- ncol(base$x_upd); n_series <- ncol(point)

# Conditional filtering variance for states and nonlinear premium functions.
filter_variance <- matrix(0, TT, n_series, dimnames = dimnames(point))
for (t in seq_len(TT)) {
  filter_variance[t, 1:K] <- diag(base$P_upd[, , t]) * scales^2
  loading <- 1200 * c(1, 1, 1, 0, 0, 0)
  filter_variance[t, 'shadow_rate'] <- drop(t(loading) %*% base$P_upd[, , t] %*% loading)
}
price <- function(coefs) y_fitting_r(t(base$x_upd), coefs$A_X_for, coefs$B_X_for,
  coefs$A_X_exp, base$objects$pars$r_lb,
  sqrt(cumsum(diag(t(coefs$B_X_for) %*% base$objects$Sigma2_X %*% coefs$B_X_for))),
  coefs$A_X_for_pi, coefs$B_X_pi, base$objects$Sigma2_X,
  coefs$B_X_cum, coefs$B_X_cum_pi)
p <- price(base$objects$coefs_p); q <- price(base$objects$coefs_q)
G <- array(0, c(TT, 6, K))
for (m in 1:2) {
  horizon <- c(60L, 120L)[m]
  G[, m, ] <- matrix(1200 * (q$JJ_n[horizon, ] - p$JJ_n[horizon, ]), TT, K, byrow = TRUE)
  G[, 2 + m, ] <- matrix(1200 * (q$JJ_r[horizon, ] - p$JJ_r[horizon, ]), TT, K, byrow = TRUE)
  inflation_loading <- 1200 * (rowSums(base$objects$coefs_q$B_X_pi[, 1:horizon]) -
                              rowSums(base$objects$coefs_p$B_X_pi[, 1:horizon])) / horizon
  G[, 4 + m, ] <- matrix(rep(inflation_loading, each = TT), TT, K)
}
for (t in seq_len(TT)) for (j in 1:6) {
  filter_variance[t, K + 1 + j] <- drop(G[t, j, ] %*% base$P_upd[, , t] %*% G[t, j, ])
}

# Independent check of premium gradients, including the layout of the batch Jacobian.
gradient_error <- 0
for (k in seq_len(K)) {
  step <- if (k == K) 1e-5 else 1e-7
  plus <- minus <- base$x_upd
  plus[, k] <- plus[, k] + step; minus[, k] <- minus[, k] - step
  numerical <- (calculate_premia(base, plus) - calculate_premia(base, minus)) / (2 * step)
  gradient_error <- max(gradient_error, abs(numerical - G[, , k]))
}
stopifnot(gradient_error < 1e-5)

calculate <- function(step_scale, frozen = TRUE) {
  worker <- function(j) {
    h <- step_scale * max(1, abs(theta[j]))
    plus <- minus <- theta
    plus[j] <- plus[j] + h; minus[j] <- minus[j] - h
    trace <- if (frozen) branches else NULL
    (economic_series(problem$evaluate_fit(plus, trace)) -
     economic_series(problem$evaluate_fit(minus, trace))) / (2 * h)
  }
  columns <- if (.Platform$OS.type == 'unix' && cores > 1L) {
    parallel::mclapply(seq_along(theta), worker, mc.cores = cores)
  } else lapply(seq_along(theta), worker)
  D <- array(unlist(columns), c(TT, n_series, length(theta)))
  parameter_variance <- matrix(0, TT, n_series, dimnames = dimnames(point))
  for (j in seq_len(n_series)) {
    derivative <- D[, j, ]
    parameter_variance[, j] <- rowSums((derivative %*% V) * derivative)
  }
  stopifnot(all(parameter_variance >= -1e-10))
  total_sd <- sqrt(filter_variance + parameter_variance)
  list(D = D, filter_variance = filter_variance, parameter_variance = parameter_variance,
       total_sd = total_sd, point = point, step_scale = step_scale, frozen = frozen)
}
result <- calculate(step_scale)
result$theta <- theta
result$parameter_covariance <- V
result$provenance <- ggnr_cache_provenance()
result$baseline_states <- point[, c("r_star", "pi_star", "w")]
result$baseline_premia <- point[, 8:13]
result$method <- "first-order delta with local parameter derivatives and conditional filtering variance"
result$active_branches_held_fixed_for_derivatives <- TRUE
result$exact_moment_targets_treated_as_fixed <- TRUE
bands <- function(columns) {
  output <- data.frame(date = dates)
  for (name in columns) {
    output[[paste0(name, "_estimate")]] <- point[, name]
    for (suffix in c("p025", "p16", "median", "p84", "p975")) {
      probability <- c(p025 = .025, p16 = .16, median = .50, p84 = .84, p975 = .975)[suffix]
      output[[paste0(name, "_", suffix)]] <- point[, name] +
        qnorm(probability) * result$total_sd[, name]
    }
  }
  output
}
result$state_bands <- bands(colnames(point)[1:7])
result$premium_bands <- bands(colnames(point)[8:13])
result$diagnostics <- data.frame(method = "local_delta", parameter_count = length(theta),
  sample_months = nrow(point), step_scale = step_scale, cores = cores,
  elapsed_seconds = unname(proc.time()["elapsed"] - started),
  reconstructed_loglik = base$fval, saved_loglik = problem$estimate$final_loglik,
  conditional_on_exact_moment_targets = TRUE, finite_sample_coverage_established = FALSE)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
# Retain a reproducible step check, refreshed whenever the primary result changes.
half <- calculate(step_scale / 2)
half$theta <- theta
half$parameter_covariance <- V
half$provenance <- result$provenance
saveRDS(half, file.path(output_directory, "half_step.rds"))
result$diagnostics$elapsed_seconds_including_half_step <-
  unname(proc.time()["elapsed"] - started)
result$diagnostics$max_premium_state_gradient_error <- gradient_error
write.csv(result$state_bands, file.path(output_directory, "state_bands.csv"), row.names = FALSE)
write.csv(result$premium_bands, file.path(output_directory, "premium_bands.csv"), row.names = FALSE)
write.csv(result$diagnostics, file.path(output_directory, "simulation_diagnostics.csv"), row.names = FALSE)
variance_decomposition <- data.frame(date = dates)
for (name in colnames(point)) {
  variance_decomposition[[paste0(name, "_filter_variance")]] <- result$filter_variance[, name]
  variance_decomposition[[paste0(name, "_parameter_variance")]] <- result$parameter_variance[, name]
  variance_decomposition[[paste0(name, "_total_se")]] <- result$total_sd[, name]
}
write.csv(variance_decomposition, file.path(output_directory, "state_variance_decomposition.csv"), row.names = FALSE)
saveRDS(result, file.path(output_directory, "delta_uncertainty.rds"))
message(sprintf("Local uncertainty bands: %d parameters; %.1f seconds.",
  length(theta), result$diagnostics$elapsed_seconds))
