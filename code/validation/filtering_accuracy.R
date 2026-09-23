# Monte Carlo assessment of EKF linearization in the six-state model.
#
# The pricing implementation is held fixed: simulated observations and both
# filters use the same approximate bond-pricing measurement function. Thus the
# exercise isolates filtering accuracy from the separate pricing-accuracy test.

source("code/model/iekf.R")
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")

estimate_file <- Sys.getenv(
  "GGNR_FILTER_ESTIMATE",
  "estimates/baseline.rds"
)
output_directory <- Sys.getenv(
  "GGNR_FILTER_OUTPUT",
  "outputs/diagnostics/filtering"
)
replication_count <- as.integer(Sys.getenv("GGNR_FILTER_REPS", "200"))
core_count <- as.integer(Sys.getenv(
  "GGNR_FILTER_CORES", as.character(min(4L, parallel::detectCores()))
))
random_seed <- as.integer(Sys.getenv("GGNR_FILTER_SEED", "1986"))
ukf_alpha <- as.numeric(Sys.getenv("GGNR_UKF_ALPHA", "1"))
iekf_trigger_probability <- as.numeric(Sys.getenv("GGNR_IEKF_TRIGGER", "0"))
iekf_steps <- as.integer(Sys.getenv("GGNR_IEKF_STEPS", "2"))
if (!is.finite(replication_count) || replication_count < 2L) {
  stop("At least two Monte Carlo replications are required.")
}
if (!is.finite(core_count) || core_count < 1L) core_count <- 1L
if (!is.finite(ukf_alpha) || ukf_alpha <= 0) stop("UKF alpha must be positive.")
if (!is.finite(iekf_trigger_probability) || iekf_trigger_probability < 0 ||
    iekf_trigger_probability > 1) stop("Invalid iterated-EKF trigger probability.")
if (!is.finite(iekf_steps) || iekf_steps < 2L) {
  stop("The iterated EKF requires at least two measurement evaluations.")
}

estimate <- readRDS(estimate_file)
p <- modifyList(reduced_hfi_no_output_defaults(), estimate$parameters)
p$include_core_cpi <- FALSE
p$sigma_inflation_survey_long <- 0.10 / 1200
p$sigma_tbill_survey_long <- 0.10 / 1200
p$sigma_real_yield <- 0.10 / 1200
p$sigma_real_yield_liquidity <- 0.20 / 1200
p$liquidity_spread_covid <- 0

data <- load_paper_data()
short_rate <- load_observed_real_short_rate_inputs(data$dates)
hfi_none <- list(
  standardized = rep(NaN, length(data$dates)),
  observed = rep(FALSE, length(data$dates)), source = "none",
  shock = NA_character_, event_count = 0L, month_count = 0L
)
reference <- run_kf_reduced_hfi_no_output(data, p, hfi_none, short_rate)
obj <- reference$objects
p <- obj$pars
K <- obj$dims$K
TT <- obj$dims$T
state_names <- reduced_hfi_no_output_states

observed_data_matrix <- function(data) {
  cbind(
    inflation = data$macro[, 2], ptr = data$macro[, 4],
    data$surv_infexp, data$surv_tbexp, data$yields_n, data$yields_r
  )
}
Y_observed <- observed_data_matrix(data)
observation_mask <- is.finite(Y_observed)

measurement_variances <- function(tt) {
  R_s <- rep(p$sigma_inflation_survey^2, obj$dims$N_s)
  if (obj$dims$N_s >= 2L) R_s[2] <- p$sigma_inflation_survey_long^2
  R_t <- rep(p$sigma_tbill_survey^2, obj$dims$N_t)
  if (obj$dims$N_t >= 2L) R_t[2] <- p$sigma_tbill_survey_long^2
  c(
    p$sigma_inf^2, p$sigma_ptr^2, R_s, R_t,
    rep(p$sigma_nominal_yield^2, obj$dims$N_n),
    rep(obj$sigma_r_t[tt]^2, obj$dims$N_r)
  )
}

measurement_at <- function(state, tt, need_jacobian = FALSE) {
  observed_tbill <- which(observation_mask[tt, 2L + obj$dims$N_s +
                                               seq_len(obj$dims$N_t)])
  tbill <- tbill_measurement_no_output(
    state, obj, data$hstep_t[observed_tbill]
  )
  y_tbill <- rep(NaN, obj$dims$N_t)
  G_tbill <- matrix(0, obj$dims$N_t, K)
  if (length(observed_tbill)) {
    y_tbill[observed_tbill] <- tbill$y_pred
    G_tbill[observed_tbill, ] <- tbill$Gamma
  }
  yield_fit <- y_fitting_r(
    state, obj$coefs_q$A_X_for, obj$coefs_q$B_X_for,
    obj$coefs_q$A_X_exp, p$r_lb, obj$s_n,
    obj$coefs_q$A_X_for_pi, obj$coefs_q$B_X_pi, obj$Sigma2_X,
    obj$coefs_q$B_X_cum, obj$coefs_q$B_X_cum_pi
  )
  value <- c(
    obj$Gamma_macro %*% state,
    obj$Gamma_s0 + obj$Gamma_s %*% state,
    y_tbill,
    yield_fit$yfit_all_n[1, data$mats_n],
    yield_fit$yfit_all_r[1, data$mats_r] + obj$real_yield_offset_t[tt]
  )
  if (!need_jacobian) return(as.numeric(value))
  jacobian <- rbind(
    obj$Gamma_macro, obj$Gamma_s, G_tbill,
    yield_fit$JJ_n[data$mats_n, , drop = FALSE],
    yield_fit$JJ_r[data$mats_r, , drop = FALSE]
  )
  list(value = as.numeric(value), jacobian = jacobian)
}

stable_chol <- function(covariance, label) {
  covariance <- (covariance + t(covariance)) / 2
  direct <- tryCatch(chol(covariance), error = function(error) NULL)
  if (!is.null(direct)) return(direct)
  scale <- max(1, max(abs(diag(covariance))))
  for (power in c(-14, -12, -10, -8, -6)) {
    candidate <- tryCatch(
      chol(covariance + diag(scale * 10^power, nrow(covariance))),
      error = function(error) NULL
    )
    if (!is.null(candidate)) return(candidate)
  }
  stop(label, " covariance is not positive definite.")
}

gaussian_log_density <- function(innovation, covariance) {
  root <- stable_chol(covariance, "innovation")
  solved <- backsolve(root, innovation, transpose = TRUE)
  -length(innovation) * log(2 * pi) / 2 - sum(log(diag(root))) -
    sum(solved^2) / 2
}

run_fixed_ekf <- function(Y) {
  mean_state <- as.numeric(solve(diag(K) - obj$Phi, obj$Mu))
  P_pred <- stationary_covariance_no_output(obj)
  x_upd <- matrix(NA_real_, TT, K, dimnames = list(NULL, state_names))
  P_upd <- array(NA_real_, c(K, K, TT),
                 dimnames = list(state_names, state_names, NULL))
  loglik_t <- numeric(TT)
  previous_mean <- mean_state
  previous_covariance <- P_pred
  for (tt in seq_len(TT)) {
    if (tt == 1L) {
      x_pred <- mean_state
    } else {
      x_pred <- as.numeric(obj$Mu + obj$Phi %*% previous_mean)
      P_pred <- obj$Phi %*% previous_covariance %*% t(obj$Phi) +
        obj$Sigma2_X
    }
    P_pred <- (P_pred + t(P_pred)) / 2
    measurement <- measurement_at(x_pred, tt, TRUE)
    keep <- observation_mask[tt, ]
    predicted <- measurement$value[keep]
    G <- measurement$jacobian[keep, , drop = FALSE]
    innovation <- Y[tt, keep] - predicted
    R <- diag(measurement_variances(tt)[keep], nrow = sum(keep))
    S <- G %*% P_pred %*% t(G) + R
    S <- (S + t(S)) / 2
    root <- stable_chol(S, "EKF innovation")
    S_inv <- chol2inv(root)
    gain <- P_pred %*% t(G) %*% S_inv
    updated_mean <- as.numeric(x_pred + gain %*% innovation)
    updated_covariance <- P_pred - gain %*% G %*% P_pred
    updated_covariance <- (updated_covariance + t(updated_covariance)) / 2
    x_upd[tt, ] <- updated_mean
    P_upd[, , tt] <- updated_covariance
    loglik_t[tt] <- gaussian_log_density(innovation, S)
    previous_mean <- updated_mean
    previous_covariance <- updated_covariance
  }
  list(x_upd = x_upd, P_upd = P_upd,
       loglik_t = loglik_t, loglik = sum(loglik_t))
}

sigma_points <- function(mean, covariance, alpha = ukf_alpha,
                         beta = 2, kappa = 0) {
  n <- length(mean)
  lambda <- alpha^2 * (n + kappa) - n
  scale <- n + lambda
  if (!is.finite(scale) || scale <= 0) stop("Invalid UKF scaling.")
  root <- t(stable_chol(scale * covariance, "UKF state"))
  points <- matrix(mean, n, 2L * n + 1L)
  points[, 1L + seq_len(n)] <- points[, 1L + seq_len(n), drop = FALSE] + root
  points[, 1L + n + seq_len(n)] <-
    points[, 1L + n + seq_len(n), drop = FALSE] - root
  weights_mean <- c(lambda / scale, rep(1 / (2 * scale), 2L * n))
  weights_covariance <- weights_mean
  weights_covariance[1] <- weights_covariance[1] + 1 - alpha^2 + beta
  list(points = points, mean = weights_mean, covariance = weights_covariance)
}

run_fixed_ukf <- function(Y) {
  mean_state <- as.numeric(solve(diag(K) - obj$Phi, obj$Mu))
  P_pred <- stationary_covariance_no_output(obj)
  x_upd <- matrix(NA_real_, TT, K, dimnames = list(NULL, state_names))
  P_upd <- array(NA_real_, c(K, K, TT),
                 dimnames = list(state_names, state_names, NULL))
  loglik_t <- numeric(TT)
  previous_mean <- mean_state
  previous_covariance <- P_pred
  for (tt in seq_len(TT)) {
    if (tt == 1L) {
      x_pred <- mean_state
    } else {
      x_pred <- as.numeric(obj$Mu + obj$Phi %*% previous_mean)
      P_pred <- obj$Phi %*% previous_covariance %*% t(obj$Phi) +
        obj$Sigma2_X
    }
    P_pred <- (P_pred + t(P_pred)) / 2
    sigma <- sigma_points(x_pred, P_pred)
    keep <- observation_mask[tt, ]
    measurement_points <- vapply(
      seq_len(ncol(sigma$points)),
      function(index) measurement_at(sigma$points[, index], tt)[keep],
      numeric(sum(keep))
    )
    if (is.null(dim(measurement_points))) {
      measurement_points <- matrix(measurement_points, nrow = sum(keep))
    }
    predicted <- as.numeric(measurement_points %*% sigma$mean)
    centered_measurement <- sweep(measurement_points, 1, predicted)
    centered_state <- sweep(sigma$points, 1, x_pred)
    weighted_measurement <- sweep(
      centered_measurement, 2, sigma$covariance, `*`
    )
    S <- weighted_measurement %*% t(centered_measurement) +
      diag(measurement_variances(tt)[keep], nrow = sum(keep))
    S <- (S + t(S)) / 2
    cross_covariance <- sweep(
      centered_state, 2, sigma$covariance, `*`
    ) %*% t(centered_measurement)
    innovation <- Y[tt, keep] - predicted
    root <- stable_chol(S, "UKF innovation")
    S_inv <- chol2inv(root)
    gain <- cross_covariance %*% S_inv
    updated_mean <- as.numeric(x_pred + gain %*% innovation)
    updated_covariance <- P_pred - gain %*% S %*% t(gain)
    updated_covariance <- (updated_covariance + t(updated_covariance)) / 2
    stable_chol(updated_covariance, "UKF updated state")
    x_upd[tt, ] <- updated_mean
    P_upd[, , tt] <- updated_covariance
    loglik_t[tt] <- gaussian_log_density(innovation, S)
    previous_mean <- updated_mean
    previous_covariance <- updated_covariance
  }
  list(x_upd = x_upd, P_upd = P_upd,
       loglik_t = loglik_t, loglik = sum(loglik_t))
}

run_fixed_iekf <- function(Y) {
  mean_state <- as.numeric(solve(diag(K) - obj$Phi, obj$Mu))
  P_pred <- stationary_covariance_no_output(obj)
  shadow_loading <- obj$loadings$shadow_rate
  x_upd <- matrix(NA_real_, TT, K, dimnames = list(NULL, state_names))
  P_upd <- array(NA_real_, c(K, K, TT),
                 dimnames = list(state_names, state_names, NULL))
  loglik_t <- numeric(TT)
  iterated <- logical(TT)
  previous_mean <- mean_state
  previous_covariance <- P_pred
  for (tt in seq_len(TT)) {
    if (tt == 1L) {
      x_pred <- mean_state
    } else {
      x_pred <- as.numeric(obj$Mu + obj$Phi %*% previous_mean)
      P_pred <- obj$Phi %*% previous_covariance %*% t(obj$Phi) +
        obj$Sigma2_X
    }
    P_pred <- (P_pred + t(P_pred)) / 2
    shadow_mean <- as.numeric(t(shadow_loading) %*% x_pred)
    shadow_variance <- as.numeric(
      t(shadow_loading) %*% P_pred %*% shadow_loading
    )
    lower_bound_probability <- pnorm(
      (p$r_lb - shadow_mean) / sqrt(shadow_variance)
    )
    step_count <- if (lower_bound_probability >= iekf_trigger_probability) {
      iekf_steps
    } else {
      1L
    }
    iterated[tt] <- step_count > 1L
    expansion_point <- x_pred
    keep <- observation_mask[tt, ]
    R <- diag(measurement_variances(tt)[keep], nrow = sum(keep))
    for (step in seq_len(step_count)) {
      measurement <- measurement_at(expansion_point, tt, TRUE)
      predicted_at_expansion <- measurement$value[keep]
      G <- measurement$jacobian[keep, , drop = FALSE]
      adjusted_innovation <- Y[tt, keep] - predicted_at_expansion +
        as.numeric(G %*% (expansion_point - x_pred))
      S <- G %*% P_pred %*% t(G) + R
      S <- (S + t(S)) / 2
      root <- stable_chol(S, "iterated-EKF innovation")
      S_inv <- chol2inv(root)
      gain <- P_pred %*% t(G) %*% S_inv
      updated_mean <- as.numeric(x_pred + gain %*% adjusted_innovation)
      if (step == 1L) {
        ordinary_innovation <- Y[tt, keep] - predicted_at_expansion
        loglik_t[tt] <- gaussian_log_density(ordinary_innovation, S)
      }
      expansion_point <- updated_mean
    }
    updated_covariance <- P_pred - gain %*% G %*% P_pred
    updated_covariance <- (updated_covariance + t(updated_covariance)) / 2
    stable_chol(updated_covariance, "iterated-EKF updated state")
    x_upd[tt, ] <- updated_mean
    P_upd[, , tt] <- updated_covariance
    previous_mean <- updated_mean
    previous_covariance <- updated_covariance
  }
  list(x_upd = x_upd, P_upd = P_upd,
       loglik_t = loglik_t, loglik = sum(loglik_t),
       iterated = iterated, iteration_fraction = mean(iterated))
}

# Verify that the stand-alone EKF reproduces the estimation filter exactly.
observed_ekf <- run_fixed_ekf(Y_observed)
validation <- c(
  max_state_difference = max(abs(observed_ekf$x_upd - reference$x_upd)),
  loglik_difference = observed_ekf$loglik - reference$fval
)
if (validation["max_state_difference"] > 1e-10 ||
    abs(validation["loglik_difference"]) > 1e-6) {
  stop(sprintf(
    paste0("The simulation EKF does not reproduce the estimation filter: ",
           "max state difference %.3e; log-likelihood difference %.3e."),
    validation["max_state_difference"], validation["loglik_difference"]
  ))
}
observed_ukf <- run_fixed_ukf(Y_observed)
observed_iekf <- run_fixed_iekf(Y_observed)
reference_iekf <- run_kf_reduced_hfi_no_output(
  data, p, hfi_none, short_rate,
  filter_method = "IEKF",
  iekf_trigger_probability = iekf_trigger_probability,
  iekf_steps = iekf_steps
)
iekf_validation <- c(
  max_state_difference = max(
    abs(observed_iekf$x_upd - reference_iekf$x_upd)
  ),
  max_covariance_difference = max(
    abs(observed_iekf$P_upd - reference_iekf$P_upd)
  ),
  loglik_difference = observed_iekf$loglik - reference_iekf$fval,
  iteration_fraction_difference = observed_iekf$iteration_fraction -
    reference_iekf$iteration_fraction
)
if (max(abs(iekf_validation)) > 1e-8) {
  stop(
    "The stand-alone IEKF does not reproduce the estimation IEKF: ",
    paste(names(iekf_validation), signif(iekf_validation, 4), collapse = ", ")
  )
}

simulate_dataset <- function(seed) {
  set.seed(seed)
  mean_state <- as.numeric(solve(diag(K) - obj$Phi, obj$Mu))
  stationary_covariance <- stationary_covariance_no_output(obj)
  initial_root <- t(stable_chol(stationary_covariance, "initial state"))
  states <- matrix(NA_real_, TT, K, dimnames = list(NULL, state_names))
  states[1, ] <- mean_state + initial_root %*% rnorm(K)
  shocks <- matrix(rnorm(TT * K), TT, K)
  if (TT > 1L) for (tt in 2:TT) {
    states[tt, ] <- as.numeric(
      obj$Mu + obj$Phi %*% states[tt - 1L, ] +
        obj$Sigma_X %*% shocks[tt, ]
    )
  }
  Y <- matrix(NaN, TT, ncol(Y_observed))
  for (tt in seq_len(TT)) {
    keep <- observation_mask[tt, ]
    Y[tt, keep] <- measurement_at(states[tt, ], tt)[keep] +
      rnorm(sum(keep), sd = sqrt(measurement_variances(tt)[keep]))
  }
  list(states = states, observations = Y)
}

series_from_filter <- function(fit) {
  shadow_loading <- obj$loadings$shadow_rate
  means <- cbind(fit$x_upd, as.numeric(fit$x_upd %*% shadow_loading))
  colnames(means) <- c(state_names, "shadow_rate")
  variances <- matrix(NA_real_, TT, K + 1L,
                      dimnames = list(NULL, c(state_names, "shadow_rate")))
  for (tt in seq_len(TT)) {
    variances[tt, state_names] <- diag(fit$P_upd[, , tt])
    variances[tt, "shadow_rate"] <- as.numeric(
      t(shadow_loading) %*% fit$P_upd[, , tt] %*% shadow_loading
    )
  }
  list(means = means, variances = variances)
}

observed_ekf_series <- series_from_filter(observed_ekf)$means
observed_ukf_series <- series_from_filter(observed_ukf)$means
observed_iekf_series <- series_from_filter(observed_iekf)$means
observed_difference <- data.frame(
  series = colnames(observed_ekf_series),
  unit = ifelse(colnames(observed_ekf_series) == "w",
                "normalized units", "percentage points"),
  mean_difference = NA_real_, mean_absolute_difference = NA_real_,
  maximum_absolute_difference = NA_real_,
  mean_absolute_difference_after_24_months = NA_real_,
  maximum_absolute_difference_after_24_months = NA_real_
)
for (index in seq_len(nrow(observed_difference))) {
  series <- observed_difference$series[index]
  scale_factor <- if (series == "w") 1 else 1200
  difference <- scale_factor *
    (observed_ukf_series[, series] - observed_ekf_series[, series])
  observed_difference$mean_difference[index] <- mean(difference)
  observed_difference$mean_absolute_difference[index] <- mean(abs(difference))
  observed_difference$maximum_absolute_difference[index] <- max(abs(difference))
  observed_difference$mean_absolute_difference_after_24_months[index] <-
    mean(abs(difference[-seq_len(24L)]))
  observed_difference$maximum_absolute_difference_after_24_months[index] <-
    max(abs(difference[-seq_len(24L)]))
}
observed_paths <- data.frame(date = data$dates)
for (series in c("r_star", "pi_star", "m", "shadow_rate")) {
  observed_paths[[paste0(series, "_ekf")]] <-
    1200 * observed_ekf_series[, series]
  observed_paths[[paste0(series, "_ukf")]] <-
    1200 * observed_ukf_series[, series]
}

inflation_compensation_states <- function(coefs, states, maturities) {
  vapply(maturities, function(maturity) {
    loading <- rowSums(coefs$B_X_pi[, seq_len(maturity), drop = FALSE]) /
      maturity
    intercept <- sum(coefs$A_X_exp_pi[seq_len(maturity)]) / maturity +
      0.5 * sum(coefs$Conv_pi[seq_len(maturity)]) / maturity
    intercept + as.numeric(states %*% loading)
  }, numeric(nrow(states)))
}

premium_series <- function(states) {
  fitted <- function(coefs) y_fitting_r(
    t(states), coefs$A_X_for, coefs$B_X_for, coefs$A_X_exp,
    p$r_lb, sqrt(cumsum(diag(t(coefs$B_X_for) %*% obj$Sigma2_X %*% coefs$B_X_for))), coefs$A_X_for_pi, coefs$B_X_pi,
    obj$Sigma2_X, coefs$B_X_cum, coefs$B_X_cum_pi
  )
  physical <- fitted(obj$coefs_p)
  risk_neutral <- fitted(obj$coefs_q)
  maturities <- c(60L, 120L)
  nominal <- 1200 * (risk_neutral$yfit_all_n[, maturities] -
                       physical$yfit_all_n[, maturities])
  real <- 1200 * (risk_neutral$yfit_all_r[, maturities] -
                    physical$yfit_all_r[, maturities])
  inflation <- 1200 * (
    inflation_compensation_states(obj$coefs_q, states, maturities) -
      inflation_compensation_states(obj$coefs_p, states, maturities)
  )
  output <- cbind(
    nominal_5y = nominal[, 1], nominal_10y = nominal[, 2],
    real_5y = real[, 1], real_10y = real[, 2],
    inflation_5y = inflation[, 1], inflation_10y = inflation[, 2]
  )
  output
}

observed_premia_ekf <- premium_series(observed_ekf$x_upd)
observed_premia_ukf <- premium_series(observed_ukf$x_upd)
observed_premia_iekf <- premium_series(observed_iekf$x_upd)
observed_premium_difference <- do.call(rbind, lapply(
  colnames(observed_premia_ekf), function(series) {
    difference <- observed_premia_ukf[, series] - observed_premia_ekf[, series]
    data.frame(
      series = series,
      mean_difference_pp = mean(difference),
      mean_absolute_difference_pp = mean(abs(difference)),
      maximum_absolute_difference_pp = max(abs(difference)),
      mean_absolute_difference_after_24_months_pp =
        mean(abs(difference[-seq_len(24L)])),
      maximum_absolute_difference_after_24_months_pp =
        max(abs(difference[-seq_len(24L)]))
    )
  }
))

pairwise_observed_difference <- function(first, second, first_name,
                                         second_name) {
  do.call(rbind, lapply(colnames(first), function(series) {
    difference <- first[, series] - second[, series]
    data.frame(
      comparison = paste0(first_name, " minus ", second_name),
      series = series,
      mean_difference = mean(difference),
      mean_absolute_difference = mean(abs(difference)),
      median_absolute_difference = median(abs(difference)),
      maximum_absolute_difference = max(abs(difference))
    )
  }))
}
observed_iekf_state_comparison <- rbind(
  pairwise_observed_difference(
    sweep(observed_iekf_series, 2,
          ifelse(colnames(observed_iekf_series) == "w", 1, 1200), `*`),
    sweep(observed_ekf_series, 2,
          ifelse(colnames(observed_ekf_series) == "w", 1, 1200), `*`),
    "IEKF", "EKF"
  ),
  pairwise_observed_difference(
    sweep(observed_iekf_series, 2,
          ifelse(colnames(observed_iekf_series) == "w", 1, 1200), `*`),
    sweep(observed_ukf_series, 2,
          ifelse(colnames(observed_ukf_series) == "w", 1, 1200), `*`),
    "IEKF", "UKF"
  )
)
observed_iekf_premium_comparison <- rbind(
  pairwise_observed_difference(
    observed_premia_iekf, observed_premia_ekf, "IEKF", "EKF"
  ),
  pairwise_observed_difference(
    observed_premia_iekf, observed_premia_ukf, "IEKF", "UKF"
  )
)

summarize_filter <- function(fit, truth, filter_name, binding) {
  derived <- series_from_filter(fit)
  truth_all <- cbind(truth, as.numeric(truth %*% obj$loadings$shadow_rate))
  colnames(truth_all) <- c(state_names, "shadow_rate")
  regimes <- list(all = rep(TRUE, TT), binding = binding,
                  nonbinding = !binding)
  rows <- list()
  index <- 0L
  for (regime_name in names(regimes)) {
    keep <- regimes[[regime_name]]
    for (series in colnames(truth_all)) {
      index <- index + 1L
      error <- derived$means[keep, series] - truth_all[keep, series]
      standard_deviation <- sqrt(pmax(derived$variances[keep, series], 0))
      scale_factor <- if (series == "w") 1 else 1200
      unit <- if (series == "w") "normalized units" else "percentage points"
      rows[[index]] <- data.frame(
        filter = filter_name, regime = regime_name, series = series,
        unit = unit,
        observations = sum(keep),
        bias = scale_factor * mean(error),
        rmse = scale_factor * sqrt(mean(error^2)),
        mae = scale_factor * mean(abs(error)),
        average_sd = scale_factor * mean(standard_deviation),
        coverage_68 = mean(abs(error) <= qnorm(0.84) * standard_deviation),
        coverage_95 = mean(abs(error) <= qnorm(0.975) * standard_deviation)
      )
    }
  }
  do.call(rbind, rows)
}

evaluate_replication <- function(seed) {
  simulated <- simulate_dataset(seed)
  binding <- as.numeric(simulated$states %*% obj$loadings$shadow_rate) < p$r_lb
  start_ekf <- proc.time()["elapsed"]
  ekf <- run_fixed_ekf(simulated$observations)
  time_ekf <- unname(proc.time()["elapsed"] - start_ekf)
  start_ukf <- proc.time()["elapsed"]
  ukf <- run_fixed_ukf(simulated$observations)
  time_ukf <- unname(proc.time()["elapsed"] - start_ukf)
  start_iekf <- proc.time()["elapsed"]
  iekf <- run_fixed_iekf(simulated$observations)
  time_iekf <- unname(proc.time()["elapsed"] - start_iekf)
  attach_diagnostics <- function(summary, fit, runtime) {
    summary$binding_fraction <- mean(binding)
    summary$loglik_per_month <- fit$loglik / TT
    summary$runtime_seconds <- runtime
    summary$iteration_fraction <- if (is.null(fit$iteration_fraction)) 0 else
      fit$iteration_fraction
    summary
  }
  rbind(
    attach_diagnostics(
      summarize_filter(ekf, simulated$states, "EKF", binding), ekf, time_ekf
    ),
    attach_diagnostics(
      summarize_filter(iekf, simulated$states, "IEKF", binding),
      iekf, time_iekf
    ),
    attach_diagnostics(
      summarize_filter(ukf, simulated$states, "UKF", binding), ukf, time_ukf
    )
  )
}

set.seed(random_seed)
seeds <- sample.int(.Machine$integer.max, replication_count)
start_time <- proc.time()["elapsed"]
if (.Platform$OS.type == "unix" && core_count > 1L) {
  replication_results <- parallel::mclapply(
    seeds, evaluate_replication, mc.cores = core_count, mc.preschedule = TRUE
  )
} else {
  replication_results <- lapply(seeds, evaluate_replication)
}
elapsed <- unname(proc.time()["elapsed"] - start_time)
replication_summary <- do.call(rbind, Map(function(value, index) {
  value$replication <- index
  value
}, replication_results, seq_along(replication_results)))

aggregate_columns <- c(
  "bias", "rmse", "mae", "average_sd",
  "coverage_68", "coverage_95", "binding_fraction",
  "loglik_per_month", "runtime_seconds", "iteration_fraction"
)
summary_mean <- aggregate(
  replication_summary[, aggregate_columns],
  replication_summary[, c("filter", "regime", "series", "unit")],
  function(value) mean(value, na.rm = TRUE)
)
summary_sd <- aggregate(
  replication_summary[, aggregate_columns],
  replication_summary[, c("filter", "regime", "series", "unit")],
  function(value) sd(value, na.rm = TRUE)
)
names(summary_sd)[-(1:4)] <- paste0(names(summary_sd)[-(1:4)], "_mc_sd")
summary_table <- merge(summary_mean, summary_sd,
                       by = c("filter", "regime", "series", "unit"),
                       sort = FALSE)

comparison <- merge(
  subset(summary_table, filter == "EKF"),
  subset(summary_table, filter == "UKF"),
  by = c("regime", "series", "unit"), suffixes = c("_ekf", "_ukf")
)
comparison$rmse_ratio_ukf_to_ekf <- comparison$rmse_ukf /
  comparison$rmse_ekf
comparison$coverage_95_difference <- comparison$coverage_95_ukf -
  comparison$coverage_95_ekf
comparison$loglik_difference_per_month <-
  comparison$loglik_per_month_ukf - comparison$loglik_per_month_ekf
iekf_comparison <- merge(
  subset(summary_table, filter == "EKF"),
  subset(summary_table, filter == "IEKF"),
  by = c("regime", "series", "unit"), suffixes = c("_ekf", "_iekf")
)
iekf_comparison$rmse_ratio_iekf_to_ekf <- iekf_comparison$rmse_iekf /
  iekf_comparison$rmse_ekf
iekf_comparison$coverage_95_difference <- iekf_comparison$coverage_95_iekf -
  iekf_comparison$coverage_95_ekf

diagnostics <- data.frame(
  replication_count = replication_count,
  sample_months = TT,
  state_dimension = K,
  ukf_sigma_points = 2L * K + 1L,
  ukf_alpha = ukf_alpha,
  iekf_steps = iekf_steps,
  iekf_trigger_probability = iekf_trigger_probability,
  seed = random_seed,
  cores = core_count,
  elapsed_seconds = elapsed,
  observed_filter_max_state_difference = validation["max_state_difference"],
  observed_filter_loglik_difference = validation["loglik_difference"],
  observed_ukf_minus_ekf_loglik = observed_ukf$loglik - observed_ekf$loglik,
  observed_ukf_minus_ekf_loglik_per_month =
    (observed_ukf$loglik - observed_ekf$loglik) / TT,
  mean_binding_fraction = mean(replication_summary$binding_fraction)
)

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(replication_summary,
          file.path(output_directory, "replication_metrics.csv"),
          row.names = FALSE)
write.csv(summary_table, file.path(output_directory, "filtering_summary.csv"),
          row.names = FALSE)
write.csv(comparison, file.path(output_directory, "filter_comparison.csv"),
          row.names = FALSE)
write.csv(iekf_comparison,
          file.path(output_directory, "iterated_ekf_comparison.csv"),
          row.names = FALSE)
write.csv(diagnostics, file.path(output_directory, "simulation_diagnostics.csv"),
          row.names = FALSE)
write.csv(observed_difference,
          file.path(output_directory, "observed_filter_comparison.csv"),
          row.names = FALSE)
write.csv(observed_paths,
          file.path(output_directory, "observed_filter_paths.csv"),
          row.names = FALSE)
write.csv(observed_premium_difference,
          file.path(output_directory, "observed_premium_comparison.csv"),
          row.names = FALSE)
write.csv(observed_iekf_state_comparison,
          file.path(output_directory, "observed_iterated_state_comparison.csv"),
          row.names = FALSE)
write.csv(observed_iekf_premium_comparison,
          file.path(output_directory, "observed_iterated_premium_comparison.csv"),
          row.names = FALSE)
saveRDS(list(
  diagnostics = diagnostics, summary = summary_table,
  comparison = comparison, iekf_comparison = iekf_comparison,
  replication_metrics = replication_summary,
  observed_difference = observed_difference,
  observed_premium_difference = observed_premium_difference,
  observed_iekf_state_comparison = observed_iekf_state_comparison,
  observed_iekf_premium_comparison = observed_iekf_premium_comparison,
  observed_paths = observed_paths,
  parameters = p, exact_targets = obj$targets
), file.path(output_directory, "filtering_accuracy.rds"))

cat(sprintf(
  "Filtering simulation: %d replications in %.1f seconds; binding %.1f%%.\n",
  replication_count, elapsed, 100 * diagnostics$mean_binding_fraction
))
print(subset(
  comparison,
  regime == "all" & series %in% c("r_star", "pi_star", "shadow_rate"),
  select = c("series", "rmse_ekf", "rmse_ukf",
             "rmse_ratio_ukf_to_ekf", "coverage_95_ekf",
             "coverage_95_ukf", "loglik_difference_per_month")
), row.names = FALSE)
cat(sprintf("\n%d-step iterated EKF relative to the ordinary EKF:\n", iekf_steps))
print(subset(
  iekf_comparison,
  regime %in% c("all", "binding") &
    series %in% c("r_star", "pi_star", "m", "shadow_rate"),
  select = c("regime", "series", "rmse_ekf", "rmse_iekf",
             "rmse_ratio_iekf_to_ekf", "bias_iekf",
             "coverage_95_ekf", "coverage_95_iekf",
             "runtime_seconds_ekf", "runtime_seconds_iekf",
             "iteration_fraction_iekf")
), row.names = FALSE)
