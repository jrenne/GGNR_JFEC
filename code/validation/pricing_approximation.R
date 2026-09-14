# Pricing-approximation validation in a transparent, stylized three-factor model.
# The current state is observed. This script deliberately contains no filter.
suppressPackageStartupMessages({
  library(Rcpp)
  library(ggplot2)
})

source("code/model/pricing_and_filter_helpers.R")
Rcpp::sourceCpp("code/validation/three_factor_pricing_mc.cpp")

output_dir <- "outputs/diagnostics/pricing"
figure_dir <- "outputs/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

rho <- c(real = 0.985, target = 0.995, cycle = 0.90)
stationary_sd_annual_pp <- c(real = 1.25, target = 1.00, cycle = 1.50)
stationary_sd <- stationary_sd_annual_pp / 1200
innovation_sd <- stationary_sd * sqrt(1 - rho^2)
Phi_q <- diag(rho)
Sigma2 <- diag(innovation_sd^2)
mu_q <- rep(0, 3)
shadow_loading <- c(1, 1, 1)
shadow_intercept <- 2 / 1200
inflation_intercept <- 2 / 1200
lower_bound <- 0
horizons <- c(12L, 24L, 60L, 120L)

grid <- expand.grid(real_level = -1:1, target_level = -1:1, cycle_level = -1:1)
states <- sweep(as.matrix(grid), 2, stationary_sd, "*")
colnames(states) <- names(rho)
grid$state_id <- seq_len(nrow(grid))
grid$shadow_rate <- 1200 * (shadow_intercept + as.vector(states %*% shadow_loading))
grid$region <- cut(grid$shadow_rate, breaks = c(-Inf, 0, 1, Inf),
                   labels = c("Below the bound", "Near the bound", "Above the bound"),
                   right = TRUE)

unconstrained_prices <- function(inflation_loading) {
  one_price <- function(x0, h, real = FALSE) {
    state_loadings <- matrix(0, h + 1L, 3)
    state_loadings[seq_len(h), ] <- matrix(-shadow_loading, h, 3, byrow = TRUE)
    intercept <- -h * shadow_intercept
    if (real) {
      state_loadings[2:(h + 1L), ] <- state_loadings[2:(h + 1L), , drop = FALSE] +
        matrix(inflation_loading, h, 3, byrow = TRUE)
      intercept <- intercept + h * inflation_intercept
    }
    powers <- vector("list", h + 1L)
    powers[[1]] <- diag(3)
    if (h >= 1L) for (j in seq_len(h)) powers[[j + 1L]] <- powers[[j]] %*% Phi_q
    mean_log <- intercept + sum(vapply(0:h, function(j) {
      sum(state_loadings[j + 1L, ] * (powers[[j + 1L]] %*% x0))
    }, numeric(1)))
    variance_log <- 0
    for (ell in seq_len(h)) {
      shock_loading <- rep(0, 3)
      for (j in ell:h) {
        shock_loading <- shock_loading +
          t(powers[[j - ell + 1L]]) %*% state_loadings[j + 1L, ]
      }
      variance_log <- variance_log + drop(t(shock_loading) %*% Sigma2 %*% shock_loading)
    }
    exp(mean_log + 0.5 * variance_log)
  }
  nominal <- outer(seq_len(nrow(states)), seq_along(horizons), Vectorize(
    function(s, j) one_price(states[s, ], horizons[j], FALSE)))
  real <- outer(seq_len(nrow(states)), seq_along(horizons), Vectorize(
    function(s, j) one_price(states[s, ], horizons[j], TRUE)))
  list(nominal = nominal, real = real)
}

run_case <- function(eta, pairs = 150000L, seed = 271828L) {
  inflation_loading <- c(0, 1, eta)
  controls <- unconstrained_prices(inflation_loading)
  coefs <- affine_coefs_r(
    Phi_q, mu_q, Sigma2, shadow_intercept, shadow_loading,
    inflation_intercept, inflation_loading, maxmat = max(horizons)
  )
  s_n <- sqrt(cumsum(diag(t(coefs$B_X_for) %*% Sigma2 %*% coefs$B_X_for)))
  approx <- y_fitting_r(
    t(states), coefs$A_X_for, coefs$B_X_for, coefs$A_X_exp,
    lower_bound, s_n, coefs$A_X_for_pi, coefs$B_X_pi, Sigma2,
    coefs$B_X_cum, coefs$B_X_cum_pi, use_cpp = FALSE
  )
  mc <- three_factor_bond_mc(
    states, rho, innovation_sd, shadow_loading, inflation_loading,
    shadow_intercept, inflation_intercept, lower_bound, horizons,
    controls$nominal, controls$real, pairs, seed
  )

  out <- do.call(rbind, lapply(seq_along(horizons), function(j) {
    h <- horizons[j]
    do.call(rbind, lapply(c("Nominal", "Real"), function(security) {
      is_nominal <- security == "Nominal"
      p <- if (is_nominal) mc$price_nominal[, j] else mc$price_real[, j]
      sep <- if (is_nominal) mc$se_price_nominal[, j] else mc$se_price_real[, j]
      ya <- if (is_nominal) approx$yfit_all_n[, h] else approx$yfit_all_r[, h]
      ym <- -1200 / h * log(p)
      # Delta-method standard error in basis points.
      se_bp <- 100 * 1200 / h * sep / p
      data.frame(grid, eta = eta, maturity_months = h,
                 maturity = paste0(h / 12, "y"), security = security,
                 approximation = ya * 1200, monte_carlo = ym,
                 error_bp = 100 * (ya * 1200 - ym), mc_se_bp = se_bp,
                 npaths = mc$npaths)
    }))
  }))
  out
}

run_simulation <- exists("run_pricing_exercise") && isTRUE(run_pricing_exercise)
if (!run_simulation) {
  message("Reusing saved Monte Carlo results...")
  baseline <- read.csv(file.path(output_dir, "state_by_state_results.csv"))
  sensitivity <- read.csv(file.path(output_dir, "covariance_sensitivity_results.csv"))
} else {
  message("Running baseline Monte Carlo benchmark...")
  baseline <- run_case(eta = 0.5, pairs = 1000000L, seed = 271828L)
  message("Running compact covariance sensitivity...")
  sensitivity <- rbind(
    run_case(eta = -0.5, pairs = 100000L, seed = 314159L),
    baseline,
    run_case(eta = 1.5, pairs = 100000L, seed = 161803L)
  )
}

write.csv(baseline, file.path(output_dir, "state_by_state_results.csv"), row.names = FALSE)
write.csv(sensitivity, file.path(output_dir, "covariance_sensitivity_results.csv"), row.names = FALSE)

summarize_errors <- function(d, grouping) {
  pieces <- split(d, interaction(d[grouping], drop = TRUE, lex.order = TRUE))
  do.call(rbind, lapply(pieces, function(z) {
    key <- z[1, grouping, drop = FALSE]
    cbind(key, data.frame(
      mean_error_bp = mean(z$error_bp),
      mean_abs_error_bp = mean(abs(z$error_bp)),
      rmse_bp = sqrt(mean(z$error_bp^2)),
      max_abs_error_bp = max(abs(z$error_bp)),
      median_mc_se_bp = median(z$mc_se_bp),
      max_mc_se_bp = max(z$mc_se_bp)
    ))
  }))
}

summary_main <- summarize_errors(baseline, c("security", "maturity_months", "maturity"))
summary_region <- summarize_errors(baseline, c("security", "region"))
summary_sensitivity <- summarize_errors(subset(sensitivity, security == "Real"),
                                        c("eta", "security", "maturity_months", "maturity"))
write.csv(summary_main, file.path(output_dir, "error_summary_by_maturity.csv"), row.names = FALSE)
write.csv(summary_region, file.path(output_dir, "error_summary_by_bound_region.csv"), row.names = FALSE)
write.csv(summary_sensitivity, file.path(output_dir, "error_summary_covariance_sensitivity.csv"), row.names = FALSE)

parameter_table <- data.frame(
  factor = c("Real-rate factor", "Inflation-target factor", "Transitory common factor"),
  persistence = rho,
  unconditional_sd_annual_pp = stationary_sd_annual_pp,
  shadow_rate_loading = shadow_loading,
  inflation_loading_baseline = c(0, 1, 0.5)
)
write.csv(parameter_table, file.path(output_dir, "calibration.csv"), row.names = FALSE)

baseline$maturity <- factor(baseline$maturity, levels = c("1y", "2y", "5y", "10y"))
p <- ggplot(baseline, aes(shadow_rate, error_bp, shape = maturity)) +
  geom_hline(yintercept = 0, linewidth = 0.55, colour = "grey55") +
  geom_point(size = 2.25, stroke = 0.65, colour = "grey25", alpha = 0.82,
             position = position_jitter(width = 0.035, height = 0)) +
  facet_wrap(~security, ncol = 2) +
  scale_shape_manual(values = c(1, 4, 3, 8)) +
  labs(x = "Current shadow short rate (percent)",
       y = "Approximation error\n(basis points)", shape = "Maturity") +
  theme_bw(base_size = 12.5) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90"),
        strip.text = element_text(face = "bold", size = 12.5),
        legend.position = "top", legend.title = element_text(size = 11.5),
        legend.text = element_text(size = 11.5), axis.title = element_text(size = 12.5),
        axis.text = element_text(size = 11.5),
        plot.margin = margin(7, 7, 7, 18))
ggsave(file.path(figure_dir, "pricing_approximation_errors.pdf"), p,
       width = 7.1, height = 3.55, device = "pdf")
ggsave(file.path(figure_dir, "pricing_approximation_errors.png"), p,
       width = 7.1, height = 3.55, dpi = 220)

capture.output({
  cat("Three-factor pricing-approximation validation\n")
  cat("Generated:", format(Sys.time(), tz = "UTC"), "UTC\n")
  cat("Baseline Monte Carlo paths:", unique(baseline$npaths), "\n\n")
  print(summary_main, row.names = FALSE)
  cat("\nBy bound region:\n")
  print(summary_region, row.names = FALSE)
  cat("\nCovariance sensitivity:\n")
  print(summary_sensitivity, row.names = FALSE)
}, file = file.path(output_dir, "summary.txt"))

print(summary_main, row.names = FALSE)
cat("\nMaximum baseline Monte Carlo SE:", max(baseline$mc_se_bp), "bp\n")
