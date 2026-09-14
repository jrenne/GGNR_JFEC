# Preferred parsimonious reduced-form/HFI specification.
#
# States: r_star, pi_star, short-rate component m, inflation gap, iid headline
# inflation, and the price-of-risk state w.  Output and trend growth are
# absent.  The model exactly matches the sample unconditional means and
# variances of headline inflation and the latent Gaussian shadow short rate
# inferred from the censored federal-funds-rate sample.

source("code/model/pricing_and_filter_helpers.R")

reduced_hfi_no_output_id <- "reduced_form_hfi_no_output_exact_moments_liquidity_v2"
reduced_hfi_no_output_states <- c(
  "r_star", "pi_star", "m", "pi_gap", "u", "w"
)
reduced_hfi_no_output_shocks <- c(
  "eps_r_star", "eps_pi_star", "eps_m", "eps_pi_gap", "eps_u", "eps_w"
)

reduced_hfi_no_output_defaults <- function() {
  list(
    rho_r_star = 0.995, mu_r_star = 0.001 / 12,
    sigma_r_star = 0.0005 / 12,
    rho_pi_star = 0.995, mu_pi_star = 0.020 / 12,
    sigma_pi_star = 0.001 / 12,
    cyclical_var = matrix(c(0.94, 0.01, -0.02, 0.70), 2, 2, byrow = TRUE,
      dimnames = list(c("m", "pi_gap"), c("m", "pi_gap"))),
    # sigma_m and sigma_u are starting values only: both are derived from the
    # two exact unconditional-variance restrictions at each evaluation.
    sigma_m = 0.004 / 12, sigma_pi_gap = 0.010 / 12,
    sigma_u = 0.020 / 12,
    # eps_m is an innovation to the reduced-form short-rate component. It is
    # not structurally identified as a monetary-policy shock. Expressing its
    # inflation loading relative to sigma_m keeps the normalization stable when
    # sigma_m is derived from the shadow-rate variance.
    impact_pi_m_ratio = -0.05,
    rho_w = 0.93, mu_w = 0, sigma_w = sqrt(1 - 0.93^2),
    sigma_wpi = 0,
    lambda_0 = setNames(rep(0, 6), reduced_hfi_no_output_shocks),
    lambda_w = setNames(rep(0, 6), reduced_hfi_no_output_shocks),
    # One parsimonious extension: the existing vector of market prices of
    # risk can also move with trend inflation (measured in annualized pp).
    lambda_pi_star_scale = 0,
    lambda_r_star_scale = 0,
    lambda_h = 0.30, sigma_h = 0.95, r_lb = 0,
    sigma_inf = 0.0003, sigma_ptr = 0.0003,
    sigma_inflation_survey = 0.0002,
    sigma_inflation_survey_long = NA_real_,
    sigma_tbill_survey = 0.00005,
    sigma_tbill_survey_long = NA_real_,
    sigma_nominal_yield = 0.00025, sigma_real_yield = 0.00025,
    sigma_real_yield_liquidity = 0.0008,
    # Additive real-yield spreads during the two observed TIPS-liquidity
    # episodes. They are expressed in the same monthly-rate units as yields.
    liquidity_spread_gfc = 0.20 / 1200,
    liquidity_spread_covid = 0,
    sigma_core_cpi = 0.20 / 1200,
    include_core_cpi = TRUE
  )
}

reduced_hfi_no_output_loadings <- function() {
  make <- function(entries) {
    ans <- setNames(rep(0, length(reduced_hfi_no_output_states)),
                    reduced_hfi_no_output_states)
    ans[names(entries)] <- entries
    ans
  }
  list(
    shadow_rate = make(c(r_star = 1, pi_star = 1, m = 1)),
    inflation = make(c(pi_star = 1, pi_gap = 1, u = 1)),
    core_inflation = make(c(pi_star = 1, pi_gap = 1)),
    inflation_target = make(c(pi_star = 1))
  )
}

transition_reduced_hfi_no_output_raw <- function(p) {
  states <- reduced_hfi_no_output_states
  shocks <- reduced_hfi_no_output_shocks
  K <- length(states)
  Mu <- matrix(0, K, 1, dimnames = list(states, NULL))
  Phi <- matrix(0, K, K, dimnames = list(states, states))
  Sigma <- matrix(0, K, length(shocks), dimnames = list(states, shocks))
  Phi["r_star", "r_star"] <- p$rho_r_star
  Mu["r_star", 1] <- (1 - p$rho_r_star) * p$mu_r_star
  Sigma["r_star", "eps_r_star"] <- p$sigma_r_star
  Phi["pi_star", "pi_star"] <- p$rho_pi_star
  Mu["pi_star", 1] <- (1 - p$rho_pi_star) * p$mu_pi_star
  Sigma["pi_star", "eps_pi_star"] <- p$sigma_pi_star
  cyclical <- c("m", "pi_gap")
  Phi[cyclical, cyclical] <- p$cyclical_var
  Sigma["m", "eps_m"] <- p$sigma_m
  Sigma["pi_gap", "eps_m"] <- p$impact_pi_m_ratio * p$sigma_m
  Sigma["pi_gap", "eps_pi_gap"] <- p$sigma_pi_gap
  Sigma["u", "eps_u"] <- p$sigma_u
  Phi["w", "w"] <- p$rho_w
  Mu["w", 1] <- (1 - p$rho_w) * p$mu_w
  Sigma["w", c("eps_pi_gap", "eps_w")] <- c(p$sigma_wpi, p$sigma_w)
  radius <- max(Mod(eigen(Phi, only.values = TRUE)$values))
  if (!is.finite(radius) || radius >= 1 - 1e-8) stop("Non-stationary dynamics.")
  beta_pi <- if (is.null(p$lambda_pi_star_scale)) 0 else p$lambda_pi_star_scale
  beta_r <- if (is.null(p$lambda_r_star_scale)) 0 else p$lambda_r_star_scale
  lambda_0_values <- p$lambda_0[shocks] -
    p$lambda_w[shocks] * 1200 *
      (beta_pi * p$mu_pi_star + beta_r * p$mu_r_star)
  lambda_0 <- matrix(lambda_0_values, ncol = 1,
                     dimnames = list(shocks, NULL))
  lambda_1 <- matrix(0, length(shocks), K,
                     dimnames = list(shocks, states))
  lambda_1[, "w"] <- p$lambda_w[shocks]
  lambda_1[, "pi_star"] <- p$lambda_w[shocks] * beta_pi * 1200
  lambda_1[, "r_star"] <- p$lambda_w[shocks] * beta_r * 1200
  list(
    Mu = Mu, Phi = Phi, Sigma = Sigma, Sigma2_X = Sigma %*% t(Sigma),
    MuQ = Mu + Sigma %*% lambda_0,
    PhiQ = Phi + Sigma %*% lambda_1,
    Lambda_0 = lambda_0, Lambda_1 = lambda_1,
    pars = p, spectral_radius = radius
  )
}

stationary_covariance_no_output <- function(model) {
  K <- nrow(model$Phi)
  covariance <- matrix(solve(
    diag(K * K) - kronecker(model$Phi, model$Phi),
    as.vector(model$Sigma2_X)
  ), K, K)
  (covariance + t(covariance)) / 2
}

effective_short_real_moments <- function(model, covariance) {
  L <- reduced_hfi_no_output_loadings()
  mean_state <- solve(diag(nrow(model$Phi)) - model$Phi, model$Mu)
  mu_s <- as.numeric(t(L$shadow_rate) %*% mean_state)
  mu_pi <- as.numeric(t(L$inflation) %*% mean_state)
  var_s <- as.numeric(t(L$shadow_rate) %*% covariance %*% L$shadow_rate)
  var_pi <- as.numeric(t(L$inflation) %*% covariance %*% L$inflation)
  cov_s_pi <- as.numeric(t(L$shadow_rate) %*% covariance %*% L$inflation)
  sd_s <- sqrt(var_s)
  d <- (mu_s - model$pars$r_lb) / sd_s
  probability <- pnorm(d)
  density <- dnorm(d)
  centered <- mu_s - model$pars$r_lb
  mean_effective_nominal <- model$pars$r_lb +
    centered * probability + sd_s * density
  second_centered <- (centered^2 + var_s) * probability +
    centered * sd_s * density
  variance_effective_nominal <- second_centered -
    (mean_effective_nominal - model$pars$r_lb)^2
  covariance_effective_inflation <- probability * cov_s_pi
  c(
    mean = mean_effective_nominal - mu_pi,
    variance = variance_effective_nominal + var_pi -
      2 * covariance_effective_inflation
  )
}

apply_reduced_hfi_exact_moments <- function(p, targets) {
  L <- reduced_hfi_no_output_loadings()
  # The four moment restrictions are algebraic conditional on the shadow-rate
  # targets. The latter are inferred once from the censored funds-rate sample
  # by reduced_hfi_exact_targets(), never inside the likelihood loop.
  p$mu_pi_star <- targets$mean_inflation
  p$mu_r_star <- targets$mean_shadow_rate - targets$mean_inflation
  p$sigma_u <- 0

  # Var(s) is affine in sigma_m^2 because the Lyapunov equation is linear in
  # the innovation covariance matrix. Evaluate its intercept and unit slope.
  p$sigma_m <- 0
  model_zero <- transition_reduced_hfi_no_output_raw(p)
  covariance_zero <- stationary_covariance_no_output(model_zero)
  variance_shadow_zero <- as.numeric(
    t(L$shadow_rate) %*% covariance_zero %*% L$shadow_rate
  )
  p$sigma_m <- 1
  model_unit <- transition_reduced_hfi_no_output_raw(p)
  covariance_unit <- stationary_covariance_no_output(model_unit)
  variance_shadow_slope <- as.numeric(
    t(L$shadow_rate) %*% (covariance_unit - covariance_zero) %*% L$shadow_rate
  )
  sigma_m_squared <-
    (targets$variance_shadow_rate - variance_shadow_zero) /
    variance_shadow_slope
  if (!is.finite(sigma_m_squared) || sigma_m_squared <= 0) {
    stop("The shadow-rate variance leaves no positive sigma_m squared.")
  }
  p$sigma_m <- sqrt(sigma_m_squared)

  model_without_u <- transition_reduced_hfi_no_output_raw(p)
  covariance_without_u <- stationary_covariance_no_output(model_without_u)
  persistent_inflation_variance <- as.numeric(
    t(L$inflation) %*% covariance_without_u %*% L$inflation
  )
  sigma_u_squared <- targets$variance_inflation - persistent_inflation_variance
  if (!is.finite(sigma_u_squared) || sigma_u_squared <= 0) {
    stop("The inflation variance leaves no positive sigma_u squared.")
  }
  p$sigma_u <- sqrt(sigma_u_squared)
  model <- transition_reduced_hfi_no_output_raw(p)
  covariance <- stationary_covariance_no_output(model)
  mean_state <- solve(diag(nrow(model$Phi)) - model$Phi, model$Mu)
  moments <- c(
    mean_inflation = as.numeric(t(L$inflation) %*% mean_state),
    variance_inflation = as.numeric(t(L$inflation) %*% covariance %*% L$inflation),
    mean_shadow_rate = as.numeric(t(L$shadow_rate) %*% mean_state),
    variance_shadow_rate = as.numeric(t(L$shadow_rate) %*% covariance %*% L$shadow_rate)
  )
  list(p = p, model = model, covariance = covariance,
       moments = moments, targets = targets)
}

infer_gaussian_shadow_moments <- function(observed_short_rate, lower_bound = 0) {
  observed <- observed_short_rate[is.finite(observed_short_rate)]
  target <- c(mean = mean(observed), variance = var(observed))
  censored_moments <- function(theta) {
    mu <- theta[1]
    sigma <- exp(theta[2])
    centered <- mu - lower_bound
    d <- centered / sigma
    probability <- pnorm(d)
    density <- dnorm(d)
    first_centered <- centered * probability + sigma * density
    second_centered <- (centered^2 + sigma^2) * probability +
      centered * sigma * density
    c(mean = lower_bound + first_centered,
      variance = second_centered - first_centered^2)
  }
  objective <- function(theta) {
    discrepancy <- (censored_moments(theta) - target) /
      pmax(abs(target), 1e-10)
    sum(discrepancy^2)
  }
  fit <- optim(c(target["mean"], log(sd(observed))), objective,
               method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  fitted <- censored_moments(fit$par)
  if (!is.finite(fit$value) || fit$value > 1e-6) {
    stop("Unable to infer Gaussian shadow-rate moments from censored data.")
  }
  c(mean = unname(fit$par[1]), variance = exp(fit$par[2])^2,
    fitted_effective_mean = fitted["mean"],
    fitted_effective_variance = fitted["variance"], objective = fit$value)
}

reduced_hfi_exact_targets <- function(data, observed_short_rate) {
  shadow <- infer_gaussian_shadow_moments(observed_short_rate, lower_bound = 0)
  list(
    mean_inflation = mean(data$macro[, 2], na.rm = TRUE),
    variance_inflation = var(data$macro[, 2], na.rm = TRUE),
    mean_shadow_rate = unname(shadow["mean"]),
    variance_shadow_rate = unname(shadow["variance"]),
    observed_effective_mean = mean(observed_short_rate, na.rm = TRUE),
    observed_effective_variance = var(observed_short_rate, na.rm = TRUE),
    shadow_inversion = shadow
  )
}

reduced_hfi_no_output_irf <- function(p, targets, horizon = 48L,
                                       impact_bp = 25) {
  restricted <- apply_reduced_hfi_exact_moments(p, targets)
  model <- restricted$model
  L <- reduced_hfi_no_output_loadings()
  impulse <- model$Sigma[, "eps_m", drop = FALSE]
  scale <- (impact_bp / 10000 / 12) /
    as.numeric(t(L$shadow_rate) %*% impulse)
  impulse <- impulse * scale
  response <- function(loading) {
    out <- numeric(horizon + 1L)
    power <- diag(nrow(model$Phi))
    for (h in 0:horizon) {
      out[h + 1L] <- 1200 * as.numeric(t(loading) %*% power %*% impulse)
      power <- power %*% model$Phi
    }
    out
  }
  data.frame(
    horizon = 0:horizon,
    shadow_rate_pp = response(L$shadow_rate),
    inflation_pp = response(L$inflation),
    real_rate_pp = response(L$shadow_rate - L$inflation)
  )
}
