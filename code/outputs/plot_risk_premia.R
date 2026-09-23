# Risk-premium diagnostics for the preferred reduced-form estimate.

source("code/model/iekf.R")
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")

estimate_file <- Sys.getenv(
  "GGNR_START_FILE",
  "estimates/baseline.rds"
)
output_directory <- Sys.getenv(
  "GGNR_PREMIA_DIR",
  "outputs/diagnostics/risk_premia"
)
figure_directory <- Sys.getenv("GGNR_FIGURE_DIR", "outputs/figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

estimate <- readRDS(estimate_file)
data <- load_paper_data()
short_rate <- load_observed_real_short_rate_inputs(data$dates)
hfi_none <- list(standardized = rep(NaN, length(data$dates)))
fit <- run_kf_reduced_hfi_no_output(
  data, estimate$parameters, hfi_none, short_rate,
  filter_method = if (is.null(estimate$filter_method)) "EKF" else estimate$filter_method,
  iekf_trigger_probability = if (is.null(estimate$iekf_trigger_probability)) 0 else estimate$iekf_trigger_probability,
  iekf_steps = if (is.null(estimate$iekf_steps)) 2L else estimate$iekf_steps
)
obj <- fit$objects
states <- t(fit$x_upd)

fitted_yields <- function(coefs) y_fitting_r(
  states, coefs$A_X_for, coefs$B_X_for, coefs$A_X_exp,
  obj$pars$r_lb, sqrt(cumsum(diag(t(coefs$B_X_for) %*% obj$Sigma2_X %*% coefs$B_X_for))), coefs$A_X_for_pi, coefs$B_X_pi,
  obj$Sigma2_X, coefs$B_X_cum, coefs$B_X_cum_pi
)
y_q <- fitted_yields(obj$coefs_q)
y_p <- fitted_yields(obj$coefs_p)
maturities <- c(`5y` = 60L, `10y` = 120L)
nominal_term_premium <- sapply(maturities, function(maturity) {
  y_q$yfit_all_n[, maturity] - y_p$yfit_all_n[, maturity]
})
real_term_premium <- sapply(maturities, function(maturity) {
  y_q$yfit_all_r[, maturity] - y_p$yfit_all_r[, maturity]
})

inflation_compensation <- function(coefs) {
  maxmat <- ncol(coefs$B_X_pi)
  out <- matrix(NA_real_, nrow(data$macro), length(maturities),
                dimnames = list(NULL, names(maturities)))
  for (j in seq_along(maturities)) {
    maturity <- maturities[j]
    loading <- rowSums(coefs$B_X_pi[, seq_len(maturity), drop = FALSE]) / maturity
    intercept <- sum(coefs$A_X_exp_pi[seq_len(maturity)]) / maturity +
      0.5 * sum(coefs$Conv_pi[seq_len(maturity)]) / maturity
    out[, j] <- intercept + as.numeric(t(loading) %*% states)
  }
  out
}
inflation_risk_premium <- inflation_compensation(obj$coefs_q) -
  inflation_compensation(obj$coefs_p)

benchmark <- readRDS("data/external/model_free_risk_premia.rds")
benchmark_acm <- benchmark$nominal_acm[match(data$dates, benchmark$date)]
benchmark_kw <- benchmark$nominal_kw[match(data$dates, benchmark$date)]
benchmark_real_dkw <- benchmark$real_dkw[match(data$dates, benchmark$date)]
benchmark_inflation_dkw <- benchmark$inflation_dkw[
  match(data$dates, benchmark$date)
]

survey_nominal_10y <- 1200 * (data$yields_n[, ncol(data$yields_n)] -
  data$surv_tbexp[, 2])
survey_real_10y <- 1200 * (data$yields_r[, ncol(data$yields_r)] -
  data$surv_tbexp[, 2] + data$surv_infexp[, 2])
survey_inflation_10y <- 1200 * (data$yields_n[, ncol(data$yields_n)] -
  data$yields_r[, ncol(data$yields_r)] - data$surv_infexp[, 2])

premium_data <- data.frame(
  date = data$dates,
  nominal_tp_5y = 1200 * nominal_term_premium[, "5y"],
  nominal_tp_10y = 1200 * nominal_term_premium[, "10y"],
  nominal_tp_10y_acm = benchmark_acm,
  nominal_tp_10y_kw = benchmark_kw,
  nominal_tp_10y_survey = survey_nominal_10y,
  real_tp_5y = 1200 * real_term_premium[, "5y"],
  real_tp_10y = 1200 * real_term_premium[, "10y"],
  real_tp_10y_dkw = benchmark_real_dkw,
  real_tp_10y_survey = survey_real_10y,
  inflation_rp_5y = 1200 * inflation_risk_premium[, "5y"],
  inflation_rp_10y = 1200 * inflation_risk_premium[, "10y"],
  inflation_rp_10y_dkw = benchmark_inflation_dkw,
  inflation_rp_10y_survey = survey_inflation_10y
)
write.csv(premium_data, file.path(output_directory, "risk_premia.csv"),
          row.names = FALSE)

colors <- c(model5 = "black", model10 = "grey35",
            benchmark1 = "grey55", benchmark2 = "grey70",
            survey = "grey35")
draw_panel <- function(model5, model10, benchmark1, survey, title,
                       benchmark1_label, benchmark2 = NULL,
                       benchmark2_label = NULL) {
  all_values <- c(model5, model10, benchmark1, benchmark2, survey)
  limits <- range(c(0, all_values), finite = TRUE)
  plot(data$dates, model5, type = "l", lwd = 2.5, col = colors["model5"],
       xlab = "", ylab = "percentage points", main = title, ylim = limits)
  lines(data$dates, model10, lwd = 2.5, lty = 2, col = colors["model10"])
  lines(data$dates, benchmark1, lwd = 2.2, lty = 3, col = colors["benchmark1"])
  if (!is.null(benchmark2)) {
    lines(data$dates, benchmark2, lwd = 2.2, lty = 4, col = colors["benchmark2"])
  }
  points(data$dates, survey, pch = 4, cex = 1.20, lwd = 2.5,
         col = colors["survey"])
  abline(h = 0, col = "grey55", lwd = 1.2)
  labels <- c("Model 5y", "Model 10y", benchmark1_label)
  line_types <- c(1, 2, 3)
  line_colors <- colors[c("model5", "model10", "benchmark1")]
  if (!is.null(benchmark2)) {
    labels <- c(labels, benchmark2_label)
    line_types <- c(line_types, 4)
    line_colors <- c(line_colors, colors["benchmark2"])
  }
  labels <- c(labels, "Survey model-free 10y")
  legend("topright", legend = labels, ncol = 2,
         col = c(line_colors, colors["survey"]),
         lty = c(line_types, NA), pch = c(rep(NA, length(line_types)), 4),
         lwd = c(rep(2.3, length(line_types)), 2.5), pt.cex = 1.25, cex = 1.02,
         bg = adjustcolor("white", alpha.f = 0.85))
}

png(file.path(figure_directory, "risk_premia.png"),
    width = 1800, height = 1800, res = 180)
par(mfrow = c(3, 1), mar = c(3.6, 5.0, 3.0, 1),
    cex.axis = 1.15, cex.lab = 1.20, cex.main = 1.35, font.main = 2,
    mgp = c(2.8, 0.85, 0), tcl = -0.35)
par(cex = 1) # Restore full-size text after the multi-panel layout.
draw_panel(
  premium_data$nominal_tp_5y, premium_data$nominal_tp_10y,
  premium_data$nominal_tp_10y_acm, premium_data$nominal_tp_10y_survey,
  "(a) Nominal term premium", "ACM 10y", premium_data$nominal_tp_10y_kw, "KW 10y"
)
draw_panel(
  premium_data$real_tp_5y, premium_data$real_tp_10y,
  premium_data$real_tp_10y_dkw, premium_data$real_tp_10y_survey,
  "(b) Real term premium", "DKW 10y"
)
draw_panel(
  premium_data$inflation_rp_5y, premium_data$inflation_rp_10y,
  premium_data$inflation_rp_10y_dkw, premium_data$inflation_rp_10y_survey,
  "(c) Inflation risk premium", "DKW 10y"
)
dev.off()

pdf(file.path(figure_directory, "risk_premia.pdf"), width = 10, height = 10)
par(mfrow = c(3, 1), mar = c(3.6, 5.0, 3.0, 1),
    cex.axis = 1.15, cex.lab = 1.20, cex.main = 1.35, font.main = 2,
    mgp = c(2.8, 0.85, 0), tcl = -0.35)
par(cex = 1) # Restore full-size text after the multi-panel layout.
draw_panel(
  premium_data$nominal_tp_5y, premium_data$nominal_tp_10y,
  premium_data$nominal_tp_10y_acm, premium_data$nominal_tp_10y_survey,
  "(a) Nominal term premium", "ACM 10y", premium_data$nominal_tp_10y_kw, "KW 10y"
)
draw_panel(
  premium_data$real_tp_5y, premium_data$real_tp_10y,
  premium_data$real_tp_10y_dkw, premium_data$real_tp_10y_survey,
  "(b) Real term premium", "DKW 10y"
)
draw_panel(
  premium_data$inflation_rp_5y, premium_data$inflation_rp_10y,
  premium_data$inflation_rp_10y_dkw, premium_data$inflation_rp_10y_survey,
  "(c) Inflation risk premium", "DKW 10y"
)
dev.off()

summary_table <- data.frame(
  series = names(premium_data)[-1],
  mean = vapply(premium_data[-1], mean, numeric(1), na.rm = TRUE),
  standard_deviation = vapply(premium_data[-1], sd, numeric(1), na.rm = TRUE),
  minimum = vapply(premium_data[-1], min, numeric(1), na.rm = TRUE),
  maximum = vapply(premium_data[-1], max, numeric(1), na.rm = TRUE),
  observations = vapply(premium_data[-1], function(x) sum(is.finite(x)), integer(1))
)
write.csv(summary_table, file.path(output_directory, "risk_premia_summary.csv"),
          row.names = FALSE)
print(summary_table)
