source("code/outputs/graphics_helpers.R")
# Robustness check that treats every observation in calendar 2020 as missing.

source("code/model/iekf.R")
source("code/model/fast_interface.R")
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")

baseline_file <- Sys.getenv(
  "GGNR_BASELINE_FILE",
  Sys.getenv("GGNR_START_FILE", "estimates/baseline_r2.rds")
)
output_directory <- Sys.getenv(
  "GGNR_COVID_DIR",
  "outputs/diagnostics/omit_2020"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

baseline <- readRDS(baseline_file)
# Omit measurements only: retain all baseline parameters and moment targets.
# This is a conditional filtering check, not a new parameter estimate.
data_full <- load_paper_data()
short_full <- load_observed_real_short_rate_inputs(data_full$dates)
omit <- format(data_full$dates, "%Y") == "2020"
fit_baseline <- evaluate_filter_fast(data_full, baseline$parameters,
  filter_method=baseline$filter_method, trigger=baseline$iekf_trigger_probability,
  steps=baseline$iekf_steps)
obj <- fit_baseline$objects
Y <- cbind(data_full$macro[,2],data_full$macro[,4],data_full$surv_infexp,
  data_full$surv_tbexp,data_full$yields_n,data_full$yields_r)
Y[omit,] <- NaN
mu <- as.vector(solve(diag(6)-obj$Phi,obj$Mu))
fit_omit <- fast_iekf_cpp(obj,Y,data_full$mats_n,data_full$mats_r,data_full$hstep_t,
  mu,stationary_covariance_no_output(obj),TRUE,baseline$iekf_trigger_probability,
  baseline$iekf_steps,NULL,FALSE)
colnames(fit_omit$x_upd) <- reduced_hfi_no_output_states
fit_omit$objects <- obj
robustness <- baseline

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

draw_panel <- function(mask, panel_title) {
  limits <- range(c(state_baseline[mask, "r_star"],
                    state_omit[mask, "r_star"]), finite = TRUE)
  ggnr_plot(dates[mask], state_baseline[mask, "r_star"], type = "l",
       col = "black", lwd = 2.6, xlab = "", ylab = "Percentage points",
       main = panel_title, ylim = limits)
  lines(dates[mask], state_omit[mask, "r_star"], col = "grey50",
        lwd = 2.4, lty = 2)
  abline(h = 0, col = "grey70", lwd = 1.0)
  grid(col = "grey90", lty = 1)
}
draw_figure <- function() {
  par(mfrow = c(2, 1), oma = c(0,0,1.8,0), mar = c(3.6, 5.0, 2.9, 1),
      cex.axis = 1.15, cex.lab = 1.20, cex.main = 1.25,
      font.main = 2, mgp = c(2.8, 0.85, 0), tcl = -0.35)
  draw_panel(rep(TRUE, length(dates)), "(a) Full sample")
  draw_panel(dates >= as.Date("2015-01-01"), "(b) 2015 to 2026")
  ggnr_figure_legend(legend = c("Baseline IEKF", "2020 omitted; parameters fixed"),
    col = c("black","grey50"), lty = c(1,2), lwd = c(2.6,2.4), ncol = 2, cex = 1.05)

}
png(file.path(output_directory, "rstar_omit2020.png"),
    width = 1800, height = 1250, res = 180)
draw_figure()
ggnr_dev_off()
figure_directory <- Sys.getenv("GGNR_FIGURE_DIR", "outputs/figures")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
ggnr_pdf(file.path(figure_directory, "rstar_omit2020.pdf"),
    width = 10, height = 6.9, useDingbats = FALSE)
draw_figure()
ggnr_dev_off()

cat("Omitted-2020 conditional filtering summary:\n")
print(summary)
print(rstar_2020)
print(parameter_comparison)
