# Observable-fit diagnostics for the preferred six-state reduced-form model.

source("code/model/iekf.R")
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
Rcpp::sourceCpp("code/model/pricing.cpp")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")

estimate_file <- Sys.getenv(
  "GGNR_START_FILE",
  "estimates/baseline.rds"
)
output_directory <- Sys.getenv(
  "GGNR_DIAGNOSTIC_DIR",
  "outputs/diagnostics/baseline"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

result <- readRDS(estimate_file)
filter_method <- if (is.null(result$filter_method)) "EKF" else
  result$filter_method
iekf_trigger_probability <- if (is.null(result$iekf_trigger_probability)) 0 else
  result$iekf_trigger_probability
iekf_steps <- if (is.null(result$iekf_steps)) 2L else result$iekf_steps
data <- load_true_release_data()
short_rate <- load_observed_real_short_rate_inputs(data$dates)
hfi_none <- list(
  standardized = rep(NaN, length(data$dates)), observed = rep(FALSE, length(data$dates)),
  source = "none", shock = NA_character_, event_count = 0L, month_count = 0L
)
fit <- run_kf_reduced_hfi_no_output(
  data, result$parameters, hfi_none, short_rate,
  filter_method = filter_method,
  iekf_trigger_probability = iekf_trigger_probability,
  iekf_steps = iekf_steps
)
obj <- fit$objects
T <- nrow(fit$x_upd)

macro_fit <- t(obj$Gamma_macro %*% t(fit$x_upd))
core_fit <- as.numeric(obj$Gamma_core %*% t(fit$x_upd))
inflation_survey_fit <- matrix(NaN, T, obj$dims$N_s)
tbill_survey_fit <- matrix(NaN, T, obj$dims$N_t)
nominal_yield_fit <- matrix(NaN, T, obj$dims$N_n)
real_yield_fit <- matrix(NaN, T, obj$dims$N_r)
for (tt in seq_len(T)) {
  state <- fit$x_upd[tt, ]
  inflation_survey_fit[tt, ] <- as.numeric(obj$Gamma_s0 + obj$Gamma_s %*% state)
  tbill_survey_fit[tt, ] <- tbill_measurement_no_output(state, obj, data$hstep_t)$y_pred
  yf <- y_fitting_r(
    state, obj$coefs_q$A_X_for, obj$coefs_q$B_X_for,
    obj$coefs_q$A_X_exp, obj$pars$r_lb, obj$s_n,
    obj$coefs_q$A_X_for_pi, obj$coefs_q$B_X_pi, obj$Sigma2_X,
    obj$coefs_q$B_X_cum, obj$coefs_q$B_X_cum_pi
  )
  nominal_yield_fit[tt, ] <- yf$yfit_all_n[1, data$mats_n]
  real_yield_fit[tt, ] <-
    yf$yfit_all_r[1, data$mats_r] + obj$real_yield_offset_t[tt]
}

rolling_12m_average <- function(x) as.numeric(stats::filter(
  x, rep(1 / 12, 12), sides = 1
))
items <- list(
  headline_cpi = list(data$macro[, 2], macro_fit[, 1]),
  headline_cpi_12m = list(
    rolling_12m_average(data$macro[, 2]),
    rolling_12m_average(macro_fit[, 1])
  ),
  ptr = list(data$macro[, 4], macro_fit[, 2]),
  cpi_forecast_1y = list(data$surv_infexp[, 1], inflation_survey_fit[, 1]),
  cpi_forecast_10y = list(data$surv_infexp[, 2], inflation_survey_fit[, 2]),
  tbill_forecast_1y = list(data$surv_tbexp[, 1], tbill_survey_fit[, 1]),
  tbill_forecast_10y = list(data$surv_tbexp[, 2], tbill_survey_fit[, 2])
)
if (isTRUE(obj$pars$include_core_cpi)) {
  items$core_cpi <- list(data$core_inflation, core_fit)
}
for (j in seq_len(ncol(data$yields_n))) {
  name <- paste0("nominal_yield_", data$mats_n[j], "m")
  items[[name]] <- list(data$yields_n[, j], nominal_yield_fit[, j])
}
for (j in seq_len(ncol(data$yields_r))) {
  name <- paste0("real_yield_", data$mats_r[j], "m")
  items[[name]] <- list(data$yields_r[, j], real_yield_fit[, j])
}

summarize_item <- function(name, item) {
  keep <- is.finite(item[[1]]) & is.finite(item[[2]])
  observed <- item[[1]][keep]
  fitted <- item[[2]][keep]
  errors <- observed - fitted
  data.frame(
    observable = name, observations = sum(keep),
    rmse_annual_pp = 1200 * sqrt(mean(errors^2)),
    mae_annual_pp = 1200 * mean(abs(errors)),
    correlation = if (sum(keep) > 2) cor(observed, fitted) else NA_real_
  )
}
fit_summary <- do.call(rbind, Map(summarize_item, names(items), items))
write.csv(fit_summary, file.path(output_directory, "fit_summary.csv"), row.names = FALSE)

relevant_survey_release <- rowSums(is.finite(data$surv_infexp)) > 0 |
  rowSums(is.finite(data$surv_tbexp)) > 0
state_changes <- rbind(NA_real_, diff(fit$x_upd))
state_change_summary <- do.call(rbind, lapply(seq_len(ncol(state_changes)), function(j) {
  state_scale <- if (colnames(state_changes)[j] == "w") 1 else 1200
  do.call(rbind, lapply(c(FALSE, TRUE), function(release) {
    values <- abs(state_scale * state_changes[relevant_survey_release == release, j])
    data.frame(
      state = colnames(state_changes)[j], survey_release = release,
      observations = sum(is.finite(values)), median_abs_change = median(values, na.rm = TRUE),
      p95_abs_change = as.numeric(quantile(values, 0.95, na.rm = TRUE)),
      maximum_abs_change = max(values, na.rm = TRUE)
    )
  }))
}))
write.csv(state_change_summary, file.path(output_directory, "state_change_summary.csv"),
          row.names = FALSE)

plot_panel <- function(item, title) {
  observed <- 1200 * item[[1]]
  fitted <- 1200 * item[[2]]
  limits <- range(c(observed, fitted), finite = TRUE)
  plot(data$dates, fitted, type = "l", col = "#0072B2", lwd = 1.5,
       xlab = "", ylab = "annualized percentage points", main = title,
       ylim = limits)
  points(data$dates, observed, col = "#D55E00", pch = 16, cex = 0.35)
}

png(file.path(output_directory, "macro_and_survey_fits.png"),
    width = 1800, height = 1500, res = 180)
par(mfrow = c(3, 2), mar = c(3, 4, 2.2, 1))
for (name in c("headline_cpi_12m", "ptr", "cpi_forecast_1y",
               "cpi_forecast_10y", "tbill_forecast_1y", "tbill_forecast_10y")) {
  plot_panel(items[[name]], gsub("_", " ", name))
}
dev.off()

png(file.path(output_directory, "filtered_states.png"),
    width = 1800, height = 1500, res = 180)
par(mfrow = c(3, 2), mar = c(3, 4, 2.2, 1))
for (j in seq_len(ncol(fit$x_upd))) {
  state_scale <- if (colnames(fit$x_upd)[j] == "w") 1 else 1200
  state_unit <- if (colnames(fit$x_upd)[j] == "w") "standardized units" else
    "annualized percentage points"
  plot(data$dates, state_scale * fit$x_upd[, j], type = "l",
       col = "#0072B2", lwd = 1, xlab = "", ylab = state_unit,
       main = gsub("_", " ", colnames(fit$x_upd)[j]))
}
dev.off()

png(file.path(output_directory, "yield_fits.png"),
    width = 1800, height = 1500, res = 180)
par(mfrow = c(3, 2), mar = c(3, 4, 2.2, 1))
for (name in c("nominal_yield_3m", "nominal_yield_24m", "nominal_yield_120m",
               "real_yield_24m", "real_yield_60m", "real_yield_120m")) {
  plot_panel(items[[name]], gsub("_", " ", name))
}
dev.off()

saveRDS(list(fit = fit, fitted = list(
  macro = macro_fit, core = core_fit, inflation_survey = inflation_survey_fit,
  tbill_survey = tbill_survey_fit, nominal_yields = nominal_yield_fit,
  real_yields = real_yield_fit
), summary = fit_summary), file.path(output_directory, "diagnostics.rds"))
print(fit_summary)
print(state_change_summary)
