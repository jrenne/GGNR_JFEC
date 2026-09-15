# Robustness check that treats every observation in calendar 2020 as missing.

source("code/model/iekf.R")
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
Rcpp::sourceCpp("code/model/pricing.cpp")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")

baseline_file <- Sys.getenv(
  "GGNR_BASELINE_FILE",
  "estimates/baseline.rds"
)
robustness_file <- Sys.getenv(
  "GGNR_ROBUSTNESS_FILE",
  "estimates/omit_2020.rds"
)
output_directory <- Sys.getenv(
  "GGNR_COVID_DIR",
  "outputs/diagnostics/omit_2020"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

baseline <- readRDS(baseline_file)
robustness <- readRDS(robustness_file)
data_full <- load_true_release_data()
short_full <- load_observed_real_short_rate_inputs(data_full$dates)
data_omit <- data_full
short_omit <- short_full
omit <- data_full$dates >= as.Date("2020-01-01") &
  data_full$dates <= as.Date("2020-12-31")
for (name in c("macro", "yields_n", "yields_r", "surv_infexp",
               "surv_gdpexp", "surv_tbexp")) {
  data_omit[[name]][omit, ] <- NaN
}
data_omit$core_inflation[omit] <- NaN
data_omit$any_survey_release[omit] <- FALSE
data_omit$any_nonmonthly_release[omit] <- FALSE
short_omit[omit] <- NaN
hfi_none <- list(standardized = rep(NaN, length(data_full$dates)))

run_saved_filter <- function(data, short_rate, estimate) {
  run_kf_reduced_hfi_no_output(
    data, estimate$parameters, hfi_none, short_rate,
    filter_method = if (is.null(estimate$filter_method)) "EKF" else
      estimate$filter_method,
    iekf_trigger_probability = if (is.null(estimate$iekf_trigger_probability))
      0 else estimate$iekf_trigger_probability,
    iekf_steps = if (is.null(estimate$iekf_steps)) 2L else estimate$iekf_steps
  )
}
fit_baseline <- run_saved_filter(data_full, short_full, baseline)
fit_omit <- run_saved_filter(data_omit, short_omit, robustness)
if (abs(fit_omit$fval - robustness$final_loglik) > 1e-6) {
  stop("The omitted-2020 filter does not reproduce the saved likelihood.")
}

calculate_premia <- function(fit) {
  states <- t(fit$x_upd)
  obj <- fit$objects
  fitted <- function(coefs) y_fitting_r(
    states, coefs$A_X_for, coefs$B_X_for, coefs$A_X_exp,
    obj$pars$r_lb, sqrt(cumsum(diag(t(coefs$B_X_for) %*% obj$Sigma2_X %*% coefs$B_X_for))), coefs$A_X_for_pi, coefs$B_X_pi,
    obj$Sigma2_X, coefs$B_X_cum, coefs$B_X_cum_pi
  )
  y_q <- fitted(obj$coefs_q)
  y_p <- fitted(obj$coefs_p)
  maturity <- 120L
  inflation_compensation <- function(coefs) {
    loading <- rowSums(coefs$B_X_pi[, seq_len(maturity), drop = FALSE]) /
      maturity
    intercept <- sum(coefs$A_X_exp_pi[seq_len(maturity)]) / maturity +
      0.5 * sum(coefs$Conv_pi[seq_len(maturity)]) / maturity
    intercept + as.numeric(t(loading) %*% states)
  }
  cbind(
    nominal_tp_10y = 1200 * (y_q$yfit_all_n[, maturity] -
                               y_p$yfit_all_n[, maturity]),
    real_tp_10y = 1200 * (y_q$yfit_all_r[, maturity] -
                            y_p$yfit_all_r[, maturity]),
    inflation_rp_10y = 1200 * (
      inflation_compensation(obj$coefs_q) -
        inflation_compensation(obj$coefs_p)
    )
  )
}
premia_baseline <- calculate_premia(fit_baseline)
premia_omit <- calculate_premia(fit_omit)

state_baseline <- 1200 * fit_baseline$x_upd[, c("r_star", "pi_star")]
state_omit <- 1200 * fit_omit$x_upd[, c("r_star", "pi_star")]
dates <- data_full$dates
post2019 <- dates >= as.Date("2019-01-01")
summary <- data.frame(
  object = c("r_star", "pi_star", colnames(premia_baseline)),
  mean_absolute_difference_full = c(
    colMeans(abs(state_omit - state_baseline)),
    colMeans(abs(premia_omit - premia_baseline))
  ),
  mean_absolute_difference_since_2019 = c(
    colMeans(abs(state_omit[post2019, ] - state_baseline[post2019, ])),
    colMeans(abs(premia_omit[post2019, ] - premia_baseline[post2019, ]))
  )
)
write.csv(summary, file.path(output_directory, "robustness_summary.csv"),
          row.names = FALSE)

rstar_2020 <- data.frame(
  estimate = c("Baseline", "Omit 2020"),
  minimum = c(min(state_baseline[omit, "r_star"]),
              min(state_omit[omit, "r_star"])),
  maximum = c(max(state_baseline[omit, "r_star"]),
              max(state_omit[omit, "r_star"])),
  peak_to_trough = c(diff(range(state_baseline[omit, "r_star"])),
                     diff(range(state_omit[omit, "r_star"])))
)
write.csv(rstar_2020, file.path(output_directory, "rstar_2020_range.csv"),
          row.names = FALSE)

selected_parameters <- c("rho_r_star", "rho_pi_star", "rho_w",
                         "liquidity_spread_gfc")
parameter_comparison <- data.frame(
  parameter = selected_parameters,
  baseline = unlist(baseline$parameters[selected_parameters]),
  omit_2020 = unlist(robustness$parameters[selected_parameters])
)
parameter_comparison$baseline[parameter_comparison$parameter ==
  "liquidity_spread_gfc"] <- 1200 * parameter_comparison$baseline[
    parameter_comparison$parameter == "liquidity_spread_gfc"]
parameter_comparison$omit_2020[parameter_comparison$parameter ==
  "liquidity_spread_gfc"] <- 1200 * parameter_comparison$omit_2020[
    parameter_comparison$parameter == "liquidity_spread_gfc"]
write.csv(parameter_comparison,
          file.path(output_directory, "parameter_comparison.csv"),
          row.names = FALSE)

draw_panel <- function(mask, panel_title, legend_position = NULL) {
  limits <- range(c(state_baseline[mask, "r_star"],
                    state_omit[mask, "r_star"]), finite = TRUE)
  plot(dates[mask], state_baseline[mask, "r_star"], type = "l",
       col = "black", lwd = 2.6, xlab = "", ylab = "Percentage points",
       main = panel_title, ylim = limits)
  lines(dates[mask], state_omit[mask, "r_star"], col = "grey50",
        lwd = 2.4, lty = 2)
  abline(h = 0, col = "grey70", lwd = 1.0)
  grid(col = "grey90", lty = 1)
  if (!is.null(legend_position)) legend(
    legend_position, c("Baseline IEKF", "All 2020 observations omitted"),
    col = c("black", "grey50"), lty = c(1, 2), lwd = c(2.6, 2.4),
    cex = 1.05, bg = adjustcolor("white", alpha.f = 0.90)
  )
}
draw_figure <- function() {
  par(mfrow = c(2, 1), mar = c(3.6, 5.0, 2.9, 1),
      cex.axis = 1.15, cex.lab = 1.20, cex.main = 1.25,
      font.main = 2, mgp = c(2.8, 0.85, 0), tcl = -0.35)
  draw_panel(rep(TRUE, length(dates)), "(a) Full sample", "bottomleft")
  draw_panel(dates >= as.Date("2015-01-01"), "(b) 2015 to 2026")
}
png(file.path(output_directory, "rstar_omit2020.png"),
    width = 1800, height = 1250, res = 180)
draw_figure()
dev.off()
pdf(file.path(output_directory, "rstar_omit2020.pdf"),
    width = 10, height = 6.9, useDingbats = FALSE)
draw_figure()
dev.off()

cat("Omitted-2020 robustness summary:\n")
print(summary)
print(rstar_2020)
print(parameter_comparison)
