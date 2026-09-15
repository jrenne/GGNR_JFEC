# Conditional real-rate/expected-inflation correlations (Mundell--Tobin effect).
# This evaluates the experiment for the preferred six-state specification.

source("code/model/iekf.R")

diagnostic_directory <- Sys.getenv(
  "GGNR_DIAGNOSTIC_DIR",
  "outputs/diagnostics/baseline"
)
output_directory <- Sys.getenv(
  "GGNR_FIGURE_DIR",
  "outputs/figures"
)
simulation_count <- as.integer(Sys.getenv("GGNR_MT_SIMULATIONS", "5000"))
forecast_horizon <- as.integer(Sys.getenv("GGNR_MT_FORECAST_HORIZON", "12"))
maturities <- c(3L, 12L, 120L)
maturity_labels <- c("3-month maturity", "1-year maturity", "10-year maturity")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

diagnostics <- readRDS(file.path(diagnostic_directory, "diagnostics.rds"))
fit <- diagnostics$fit
obj <- fit$objects
states <- t(fit$x_upd)
dates <- load_true_release_data()$dates
K <- obj$dims$K

# Distribution of X_{t+h} conditional on the filtered state X_t.
Phi_h <- diag(K)
mu_h <- rep(0, K)
cov_h <- matrix(0, K, K)
Phi_power <- diag(K)
for (h in seq_len(forecast_horizon)) {
  mu_h <- obj$Mu + obj$Phi %*% mu_h
  cov_h <- cov_h + Phi_power %*% obj$Sigma2_X %*% t(Phi_power)
  Phi_power <- obj$Phi %*% Phi_power
}
Phi_h <- Phi_power
cov_h <- (cov_h + t(cov_h)) / 2

set.seed(20260908)
standard_draws <- matrix(rnorm(K * simulation_count), K, simulation_count)
conditional_innovations <- t(chol(cov_h)) %*% standard_draws

# Physical-measure average expected inflation at each requested maturity.
expected_inflation_loading <- matrix(0, length(maturities), K)
expected_inflation_intercept <- numeric(length(maturities))
for (j in seq_along(maturities)) {
  maturity <- maturities[j]
  expected_inflation_loading[j, ] <- rowSums(
    obj$coefs_p$B_X_pi[, seq_len(maturity), drop = FALSE]
  ) / maturity
  expected_inflation_intercept[j] <- sum(
    obj$coefs_p$A_X_exp_pi[seq_len(maturity)] / maturity +
      0.5 * obj$coefs_p$Conv_pi[seq_len(maturity)] / maturity^2
  )
}

# Pricing formula without Jacobians, which makes the simulation inexpensive.
selected_real_yields <- function(X) {
  coefs <- obj$coefs_q
  max_maturity <- max(maturities)
  sample_size <- ncol(X)
  shadow_forwards <- sweep(
    t(X) %*% coefs$B_X_for[, seq_len(max_maturity), drop = FALSE],
    2, as.vector(coefs$A_X_exp[seq_len(max_maturity)]), "+"
  )
  fitted_forwards <- shadow_forwards
  fitted_forwards[, 1] <- pmax(shadow_forwards[, 1], obj$pars$r_lb)
  if (max_maturity >= 2L) {
    mu <- sweep(shadow_forwards[, 2:max_maturity, drop = FALSE],
                2, obj$pars$r_lb, "-")
    probabilities <- sweep(mu, 2, obj$s_n[seq_len(max_maturity - 1L)], "/")
    probabilities <- pnorm(probabilities)
    nominal_mu <- sweep(mu, 2, as.vector(coefs$A_X_for[2:max_maturity] -
                                          coefs$A_X_exp[2:max_maturity]), "+")
    nominal_z <- sweep(nominal_mu, 2, obj$s_n[seq_len(max_maturity - 1L)], "/")
    fitted_forwards[, 2:max_maturity] <- obj$pars$r_lb +
      nominal_mu * pnorm(nominal_z) +
      sweep(dnorm(nominal_z), 2, obj$s_n[seq_len(max_maturity - 1L)], "*")
  }
  inflation_forwards <- sweep(
    t(X) %*% coefs$B_X_pi[, seq_len(max_maturity), drop = FALSE],
    2, as.vector(coefs$A_X_for_pi[seq_len(max_maturity)]), "+"
  )
  probabilities_all <- cbind(
    as.numeric(shadow_forwards[, 1] > obj$pars$r_lb), probabilities
  )
  covariance_correction <- c(
    0,
    colSums(
      coefs$B_X_cum[, 2:max_maturity, drop = FALSE] *
        (obj$Sigma2_X %*% coefs$B_X_cum_pi[, 2:max_maturity, drop = FALSE])
    )
  )
  real_forwards <- fitted_forwards - inflation_forwards +
    sweep(probabilities_all, 2, covariance_correction, "*")
  cumulative <- t(apply(real_forwards, 1, cumsum))
  sapply(maturities, function(maturity) cumulative[, maturity] / maturity)
}

conditional_correlations <- matrix(
  NA_real_, nrow = ncol(states), ncol = length(maturities),
  dimnames = list(NULL, paste0(maturities, "m"))
)
# Check the inexpensive simulation formula against the common pricing routine.
check_states <- states[, seq_len(min(3L, ncol(states))), drop = FALSE]
check_full <- y_fitting_r(check_states, obj$coefs_q$A_X_for,
  obj$coefs_q$B_X_for, obj$coefs_q$A_X_exp, obj$pars$r_lb, obj$s_n,
  obj$coefs_q$A_X_for_pi, obj$coefs_q$B_X_pi, obj$Sigma2_X,
  obj$coefs_q$B_X_cum, obj$coefs_q$B_X_cum_pi)
stopifnot(max(abs(selected_real_yields(check_states) -
                   check_full$yfit_all_r[, maturities])) < 1e-12)
for (tt in seq_len(ncol(states))) {
  conditional_mean <- as.numeric(mu_h + Phi_h %*% states[, tt])
  simulated_states <- conditional_innovations + conditional_mean
  expected_inflation <- sweep(
    t(expected_inflation_loading %*% simulated_states),
    2, expected_inflation_intercept, "+"
  )
  real_yields <- selected_real_yields(simulated_states)
  conditional_correlations[tt, ] <- vapply(
    seq_along(maturities),
    function(j) cor(expected_inflation[, j], real_yields[, j]),
    numeric(1)
  )
}

output <- data.frame(date = dates, conditional_correlations, check.names = FALSE)
write.csv(output, file.path(output_directory, "mundell_tobin_effect.csv"),
          row.names = FALSE)

draw_figure <- function() {
  colors <- c("black", "grey35", "grey65")
  correlation_range <- range(c(0, conditional_correlations), finite = TRUE)
  correlation_padding <- 0.06 * diff(correlation_range)
  matplot(
    dates, conditional_correlations, type = "l",
    lwd = c(2.8, 2.6, 2.6), lty = c(1, 2, 3), col = colors,
    xlab = "", ylab = "Conditional correlation",
    ylim = correlation_range + c(-correlation_padding, correlation_padding)
  )
  abline(h = 0, col = "grey60", lwd = 1.2)
  grid(col = "grey88", lty = 1)
  legend(
    "bottomleft", maturity_labels,
    col = colors, lty = c(1, 2, 3), lwd = c(2.8, 2.6, 2.6),
    bg = adjustcolor("white", alpha.f = 0.9), cex = 1.08
  )
}

png(file.path(output_directory, "mundell_tobin_effect.png"),
    width = 1800, height = 1050, res = 180)
par(mar = c(4.0, 5.0, 1.0, 1.2), cex.axis = 1.18, cex.lab = 1.22,
    mgp = c(2.8, 0.85, 0), tcl = -0.35)
draw_figure()
dev.off()

pdf(file.path(output_directory, "mundell_tobin_effect.pdf"),
    width = 10, height = 5.8, useDingbats = FALSE)
par(mar = c(4.0, 5.0, 1.0, 1.2), cex.axis = 1.18, cex.lab = 1.22,
    mgp = c(2.8, 0.85, 0), tcl = -0.35)
draw_figure()
dev.off()

print(data.frame(
  maturity_months = maturities,
  mean_correlation = colMeans(conditional_correlations),
  minimum = apply(conditional_correlations, 2, min),
  maximum = apply(conditional_correlations, 2, max)
))
