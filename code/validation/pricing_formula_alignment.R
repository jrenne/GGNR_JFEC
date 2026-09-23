# Run from the project root in RStudio. These checks compare the pricing
# routine with the paper's formula, C++ with R, and Jacobians with finite
# differences. They do not change estimates or empirical results.
source("code/model/iekf.R")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")

check_pricing <- function(Phi, Mu, V, x, inflation_loading,
                          shadow_intercept = 0, inflation_intercept = 0) {
  K <- nrow(Phi)
  b <- c(1, 1, 1, rep(0, K - 3))
  co <- affine_coefs_r(Phi, Mu, V, shadow_intercept, b,
                       inflation_intercept, inflation_loading, maxmat = 120)
  sd <- sqrt(cumsum(diag(t(co$B_X_for) %*% V %*% co$B_X_for)))
  evaluate <- function(state, cpp = FALSE, zero_inflation = FALSE) {
    cc <- if (zero_inflation) affine_coefs_r(
      Phi, Mu, V, shadow_intercept, b, 0, rep(0, K), maxmat = 120
    ) else co
    y_fitting_r(state, cc$A_X_for, cc$B_X_for, cc$A_X_exp, 0, sd,
                cc$A_X_for_pi, cc$B_X_pi, V, cc$B_X_cum,
                cc$B_X_cum_pi, use_cpp = cpp)
  }
  r <- evaluate(x)
  cpp <- evaluate(x, TRUE)
  implementation_error <- max(abs(c(r$yfit_all_n - cpp$yfit_all_n,
                                    r$yfit_all_r - cpp$yfit_all_r,
                                    r$JJ_n - cpp$JJ_n,
                                    r$JJ_r - cpp$JJ_r)))

  # Independent direct evaluation of the two forward-rate terms. The first
  # current short rate is observed conditional on x and is handled exactly.
  g <- function(z) z * pnorm(z) + dnorm(z)
  fn <- fr <- matrix(0, ncol(x), 120)
  fn[, 1] <- pmax(shadow_intercept + as.vector(t(x) %*% b), 0)
  fr[, 1] <- fn[, 1] - as.vector(co$A_X_for_pi[1] +
                                  t(x) %*% co$B_X_pi[, 1])
  for (h in 2:120) {
    a <- as.vector(co$A_X_for[h] + t(x) %*% co$B_X_for[, h])
    abar <- as.vector(co$A_X_exp[h] + t(x) %*% co$B_X_for[, h])
    fn[, h] <- sd[h - 1] * g(a / sd[h - 1])
    interaction <- drop(t(co$B_X_cum[, h]) %*% V %*% co$B_X_cum_pi[, h])
    fr[, h] <- fn[, h] - as.vector(co$A_X_for_pi[h] +
                                    t(x) %*% co$B_X_pi[, h]) +
      pnorm(abar / sd[h - 1]) * interaction
  }
  yn <- sweep(t(apply(fn, 1, cumsum)), 2, 1:120, "/")
  yr <- sweep(t(apply(fr, 1, cumsum)), 2, 1:120, "/")
  formula_error <- max(abs(c(r$yfit_all_n - yn, r$yfit_all_r - yr)))
  zero <- evaluate(x, zero_inflation = TRUE)
  zero_error <- max(abs(c(zero$yfit_all_n - zero$yfit_all_r,
                          zero$yfit_all_n - yn)))
  derivative_error <- 0
  for (tt in seq_len(ncol(x))) {
    for (k in seq_len(K)) {
      delta <- matrix(0, K, 1)
      delta[k] <- 1e-8
      plus <- evaluate(x[, tt, drop = FALSE] + delta)
      minus <- evaluate(x[, tt, drop = FALSE] - delta)
      numerical_n <- as.vector((plus$yfit_all_n - minus$yfit_all_n) / 2e-8)
      numerical_r <- as.vector((plus$yfit_all_r - minus$yfit_all_r) / 2e-8)
      column <- (tt - 1) * K + k
      derivative_error <- max(derivative_error,
        abs(numerical_n - r$JJ_n[, column]),
        abs(numerical_r - r$JJ_r[, column]))
    }
  }
  # Independent conditional mean/variance of cumulative inflation. This
  # checks the Jensen term and the per-period normalization used for IRP.
  inflation_error <- 0
  for (h in c(12L, 60L, 120L)) {
    expected_state <- x[, 1]
    mean_sum <- variance_sum <- 0
    shock_loading <- rep(0, K)
    for (j in seq_len(h)) {
      expected_state <- Mu + Phi %*% expected_state
      mean_sum <- mean_sum + inflation_intercept +
        sum(inflation_loading * expected_state)
      shock_loading <- inflation_loading + t(Phi) %*% shock_loading
      variance_sum <- variance_sum + drop(t(shock_loading) %*% V %*% shock_loading)
    }
    exact <- (mean_sum + 0.5 * variance_sum) / h
    implemented <- sum(co$A_X_exp_pi[seq_len(h)]) / h +
      0.5 * sum(co$Conv_pi[seq_len(h)]) / h +
      sum(x[, 1] * rowSums(co$B_X_pi[, seq_len(h), drop = FALSE])) / h
    inflation_error <- max(inflation_error, abs(exact - implemented))
  }
  stopifnot(implementation_error < 1e-11, formula_error < 1e-12,
            inflation_error < 1e-12,
            zero_error < 1e-12, derivative_error < 1e-5)
  c(cpp_r = implementation_error, formula = formula_error,
    zero_inflation = zero_error, derivative = derivative_error,
    cumulative_inflation = inflation_error)
}

rho <- c(.985, .995, .90)
stationary_sd <- c(1.25, 1, 1.5) / 1200
V <- diag((stationary_sd * sqrt(1 - rho^2))^2)
x <- t(sweep(as.matrix(expand.grid(a = -1:1, b = -1:1, c = -1:1)),
             2, stationary_sd, "*"))
results <- list(three_factor = check_pricing(
  diag(rho), rep(0, 3), V, x, c(0, 1, .5), 2 / 1200, 2 / 1200
))

estimate <- readRDS("estimates/baseline.rds")
data <- load_paper_data()
short_rate <- load_observed_real_short_rate_inputs(data$dates)
fit <- run_kf_reduced_hfi_no_output(data, estimate$parameters,
  no_hfi_data(data$dates), short_rate, filter_method = "IEKF",
  iekf_trigger_probability = 0, iekf_steps = 2)
obj <- fit$objects
selected <- unique(round(seq(1, nrow(fit$x_upd), length.out = 12)))
results$estimated_physical <- check_pricing(
  obj$Phi, obj$Mu, obj$Sigma2_X, t(fit$x_upd[selected, ]),
  reduced_hfi_no_output_loadings()$inflation
)
results$estimated_risk_neutral <- check_pricing(
  obj$PhiQ, obj$MuQ, obj$Sigma2_X, t(fit$x_upd[selected, ]),
  reduced_hfi_no_output_loadings()$inflation
)
results <- do.call(rbind, results)
print(results)

# The reported premium paths must use the same prices and measure-specific
# variances, including the cumulative-inflation normalization.
price <- function(co) {
  sd <- sqrt(cumsum(diag(t(co$B_X_for) %*% obj$Sigma2_X %*% co$B_X_for)))
  y_fitting_r(t(fit$x_upd), co$A_X_for, co$B_X_for, co$A_X_exp,
    obj$pars$r_lb, sd, co$A_X_for_pi, co$B_X_pi, obj$Sigma2_X,
    co$B_X_cum, co$B_X_cum_pi)
}
yp <- price(obj$coefs_p)
yq <- price(obj$coefs_q)
published <- read.csv("outputs/diagnostics/risk_premia/risk_premia.csv")
bands <- read.csv("outputs/diagnostics/uncertainty/premium_bands.csv")
inflation_ce <- function(co, h) {
  (sum(co$A_X_for_pi[seq_len(h)]) +
     as.vector(fit$x_upd %*% rowSums(co$B_X_pi[, seq_len(h), drop = FALSE]))) / h
}
for (h in c(60L, 120L)) {
  suffix <- paste0(h / 12, "y")
  nominal <- 1200 * (yq$yfit_all_n[, h] - yp$yfit_all_n[, h])
  real <- 1200 * (yq$yfit_all_r[, h] - yp$yfit_all_r[, h])
  inflation <- 1200 * (inflation_ce(obj$coefs_q, h) - inflation_ce(obj$coefs_p, h))
  stopifnot(max(abs(published[[paste0("nominal_tp_", suffix)]] - nominal)) < 1e-9,
    max(abs(published[[paste0("real_tp_", suffix)]] - real)) < 1e-9,
    max(abs(published[[paste0("inflation_rp_", suffix)]] - inflation)) < 1e-9,
    max(abs(bands[[paste0("nominal_", suffix, "_estimate")]] - nominal)) < 1e-9,
    max(abs(bands[[paste0("real_", suffix, "_estimate")]] - real)) < 1e-9,
    max(abs(bands[[paste0("inflation_", suffix, "_estimate")]] - inflation)) < 1e-9)
}
message("Reported premium paths and Hamilton point estimates agree with the pricing formulas.")
write.csv(results, "outputs/diagnostics/pricing/formula_alignment_checks.csv")
