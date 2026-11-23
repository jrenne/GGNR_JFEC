# Prices and premia used in both local and simulation uncertainty calculations.
inflation_compensation <- function(coefs, states, maturities) {
  vapply(maturities, function(maturity) {
    loading <- rowSums(coefs$B_X_pi[, seq_len(maturity), drop = FALSE]) /
      maturity
    intercept <- sum(coefs$A_X_exp_pi[seq_len(maturity)]) / maturity +
      0.5 * sum(coefs$Conv_pi[seq_len(maturity)]) / maturity
    intercept + as.numeric(states %*% loading)
  }, numeric(nrow(states)))
}

calculate_premia <- function(fit, states) {
  fitted <- function(coefs) y_fitting_r(
    t(states), coefs$A_X_for, coefs$B_X_for, coefs$A_X_exp,
    fit$objects$pars$r_lb, sqrt(cumsum(diag(t(coefs$B_X_for) %*% fit$objects$Sigma2_X %*% coefs$B_X_for))), coefs$A_X_for_pi,
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

