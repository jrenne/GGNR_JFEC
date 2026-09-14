# Robustness check using real-yield observations only in ordinary-liquidity
# post-2004 months. All other observables remain available throughout.

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
  "estimates/liquid_tips_only.rds"
)
output_directory <- Sys.getenv(
  "GGNR_LIQUID_TIPS_DIR",
  "outputs/diagnostics/liquid_tips"
)
table_directory <- Sys.getenv(
  "GGNR_TABLE_DIR", "outputs/tables"
)
figure_directory <- Sys.getenv(
  "GGNR_FIGURE_DIR", "outputs/figures"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

baseline <- readRDS(baseline_file)
robustness <- readRDS(robustness_file)
if (!isTRUE(robustness$liquid_tips_only)) {
  stop("The robustness estimate was not produced with liquid TIPS only.")
}

data_full <- load_true_release_data()
short_rate <- load_observed_real_short_rate_inputs(data_full$dates)
dates <- data_full$dates
gfc <- dates >= as.Date("2008-01-01") & dates < as.Date("2010-01-01")
covid <- dates >= as.Date("2020-03-01") & dates < as.Date("2021-03-01")
ordinary_post2004 <- dates >= as.Date("2004-01-01") & !gfc & !covid
excluded_episodes <- gfc | covid
data_liquid <- data_full
data_liquid$yields_r[!ordinary_post2004, ] <- NaN
hfi_none <- list(
  standardized = rep(NaN, length(dates)), observed = rep(FALSE, length(dates)),
  source = "none", shock = NA_character_, event_count = 0L, month_count = 0L
)

run_saved_filter <- function(data, estimate) {
  run_kf_reduced_hfi_no_output(
    data, estimate$parameters, hfi_none, short_rate,
    filter_method = if (is.null(estimate$filter_method)) "IEKF" else
      estimate$filter_method,
    iekf_trigger_probability = if (is.null(estimate$iekf_trigger_probability))
      0 else estimate$iekf_trigger_probability,
    iekf_steps = if (is.null(estimate$iekf_steps)) 2L else estimate$iekf_steps
  )
}
fit_baseline <- run_saved_filter(data_full, baseline)
fit_robust <- run_saved_filter(data_liquid, robustness)
if (abs(fit_robust$fval - robustness$final_loglik) > 1e-6) {
  stop("The liquid-TIPS filter does not reproduce the saved likelihood.")
}

calculate_premia <- function(fit) {
  states <- t(fit$x_upd)
  obj <- fit$objects
  fitted <- function(coefs) y_fitting_r(
    states, coefs$A_X_for, coefs$B_X_for, coefs$A_X_exp,
    obj$pars$r_lb, obj$s_n, coefs$A_X_for_pi, coefs$B_X_pi,
    obj$Sigma2_X, coefs$B_X_cum, coefs$B_X_cum_pi
  )
  y_q <- fitted(obj$coefs_q)
  y_p <- fitted(obj$coefs_p)
  maturity <- 120L
  inflation_compensation <- function(coefs) {
    loading <- rowSums(coefs$B_X_pi[, seq_len(maturity), drop = FALSE]) /
      maturity
    intercept <- sum(coefs$A_X_exp_pi[seq_len(maturity)]) / maturity +
      0.5 * sum(coefs$Conv_pi[seq_len(maturity)]) / maturity^2
    intercept + as.numeric(t(loading) %*% states)
  }
  cbind(
    real_yield_10y = 1200 * y_q$yfit_all_r[, maturity],
    nominal_tp = 1200 * (y_q$yfit_all_n[, maturity] -
                           y_p$yfit_all_n[, maturity]),
    real_tp_5y = 1200 * (y_q$yfit_all_r[, 60L] -
                           y_p$yfit_all_r[, 60L]),
    real_tp = 1200 * (y_q$yfit_all_r[, maturity] -
                        y_p$yfit_all_r[, maturity]),
    inflation_rp = 1200 * (
      inflation_compensation(obj$coefs_q) -
        inflation_compensation(obj$coefs_p)
    )
  )
}

baseline_objects <- cbind(
  r_star = 1200 * fit_baseline$x_upd[, "r_star"],
  pi_star = 1200 * fit_baseline$x_upd[, "pi_star"],
  calculate_premia(fit_baseline)
)
robust_objects <- cbind(
  r_star = 1200 * fit_robust$x_upd[, "r_star"],
  pi_star = 1200 * fit_robust$x_upd[, "pi_star"],
  calculate_premia(fit_robust)
)
extended_real_yield_10y <- 1200 * data_full$yields_r[, ncol(data_full$yields_r)]
model_free_real_tp_10y <- 1200 * (
  data_full$yields_r[, ncol(data_full$yields_r)] -
    data_full$surv_tbexp[, 2] + data_full$surv_infexp[, 2]
)
difference <- robust_objects - baseline_objects

masks <- list(
  full_sample = rep(TRUE, length(dates)),
  pre2004 = dates < as.Date("2004-01-01"),
  ordinary_post2004 = ordinary_post2004,
  omitted_gfc_covid = excluded_episodes
)
summary_rows <- lapply(names(masks), function(sample_name) {
  mask <- masks[[sample_name]]
  data.frame(
    sample = sample_name,
    object = colnames(difference),
    mean_difference = colMeans(difference[mask, , drop = FALSE]),
    mean_absolute_difference = colMeans(abs(difference[mask, , drop = FALSE])),
    maximum_absolute_difference = apply(
      abs(difference[mask, , drop = FALSE]), 2, max
    ),
    months = sum(mask), row.names = NULL
  )
})
summary <- do.call(rbind, summary_rows)
write.csv(summary, file.path(output_directory, "robustness_summary.csv"),
          row.names = FALSE)

series <- data.frame(
  date = dates,
  baseline_objects,
  setNames(as.data.frame(robust_objects),
           paste0(colnames(robust_objects), "_liquid_tips"))
)
write.csv(series, file.path(output_directory, "robustness_series.csv"),
          row.names = FALSE)

selected_parameters <- c("rho_r_star", "rho_pi_star", "rho_w",
                         "lambda_r_star_scale", "lambda_pi_star_scale")
parameter_comparison <- data.frame(
  parameter = selected_parameters,
  baseline = unlist(baseline$parameters[selected_parameters]),
  liquid_tips = unlist(robustness$parameters[selected_parameters]),
  row.names = NULL
)
write.csv(parameter_comparison,
          file.path(output_directory, "parameter_comparison.csv"),
          row.names = FALSE)

labels <- c(
  r_star = "$r_t^*$", pi_star = "$\\pi_t^*$",
  nominal_tp = "Nominal term premium",
  real_tp = "Real term premium",
  inflation_rp = "Inflation risk premium"
)
ordinary <- summary[summary$sample == "ordinary_post2004", ]
episodes <- summary[summary$sample == "omitted_gfc_covid", ]
table_lines <- c(
  "\\begin{table}[H]", "\\centering",
  "\\caption{Robustness to illiquid real-yield observations}",
  "\\label{tab:liquid_tips_robustness}",
  "\\begin{threeparttable}", "\\small",
  "\\begin{tabular}{lrrrr}", "\\toprule",
  " & \\multicolumn{2}{c}{Ordinary post-2004 months} & \\multicolumn{2}{c}{Omitted GFC/COVID months} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  " & Mean diff. & Mean abs. diff. & Mean diff. & Mean abs. diff. \\\\",
  "\\midrule"
)
for (object_name in names(labels)) {
  i <- match(object_name, ordinary$object)
  j <- match(object_name, episodes$object)
  table_lines <- c(table_lines, sprintf(
    "%s & %.3f & %.3f & %.3f & %.3f \\\\", labels[[object_name]],
    ordinary$mean_difference[i], ordinary$mean_absolute_difference[i],
    episodes$mean_difference[j], episodes$mean_absolute_difference[j]
  ))
}
table_lines <- c(
  table_lines, "\\bottomrule", "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]", "\\footnotesize",
  paste0(
    "\\item Notes: Entries are annualized percentage points and compare the ",
    "liquid-TIPS re-estimate with the baseline. The sensitivity estimate uses ",
    "real-yield observations only in ordinary-liquidity months from 2004 onward; ",
    "all other observables remain available throughout. The GFC and acute-COVID ",
    "windows contain 36 months. `Mean diff.' is sensitivity minus baseline."
  ),
  "\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}"
)
writeLines(table_lines,
           file.path(table_directory, "table_liquid_tips_robustness.tex"))

draw_sensitivity_panel <- function(baseline_series, sensitivity_series,
                                   panel_title, comparison_series = NULL,
                                   comparison_type = c("line", "points", "none"),
                                   comparison_label = NULL,
                                   show_model_legend = FALSE) {
  comparison_type <- match.arg(comparison_type)
  limits <- range(c(baseline_series, sensitivity_series, comparison_series),
                  finite = TRUE)
  plot(dates, baseline_series, type = "n", ylim = limits, xlab = "",
       ylab = "Annualized percentage points", main = panel_title)
  usr <- par("usr")
  rect(min(dates), usr[3], as.Date("2003-12-31"), usr[4],
       col = "grey94", border = NA)
  grid(col = "grey85", lty = 1)
  if (comparison_type == "line") {
    lines(dates, comparison_series, col = "grey70", lwd = 1.9)
  } else if (comparison_type == "points") {
    points(dates, comparison_series, col = "grey65", pch = 4,
           cex = 1.15, lwd = 2.2)
  }
  lines(dates, baseline_series, col = "black", lwd = 2.7)
  lines(dates, sensitivity_series, col = "grey35", lwd = 2.5, lty = 2)
  abline(v = as.Date("2004-01-01"), col = "grey25", lwd = 1.2, lty = 3)
  if (show_model_legend) {
    legend_labels <- c("Baseline", "Sensitivity estimate")
    legend_colors <- c("black", "grey35")
    legend_lty <- c(1, 2)
    legend_pch <- c(NA, NA)
    legend_lwd <- c(2.7, 2.5)
    if (!is.null(comparison_series)) {
      legend_labels <- c(legend_labels, comparison_label)
      legend_colors <- c(legend_colors, "grey70")
      legend_lty <- c(legend_lty, if (comparison_type == "line") 1 else NA)
      legend_pch <- c(legend_pch, if (comparison_type == "points") 4 else NA)
      legend_lwd <- c(legend_lwd, if (comparison_type == "line") 1.9 else 2.2)
    }
    legend("topright", legend_labels, col = legend_colors,
           lty = legend_lty, pch = legend_pch, lwd = legend_lwd,
           pt.cex = 1.15, cex = 1.03,
           bg = adjustcolor("white", alpha.f = 0.92))
  } else if (!is.null(comparison_series)) {
    legend("topright", comparison_label, col = "grey65",
           lty = if (comparison_type == "line") 1 else NA,
           pch = if (comparison_type == "points") 4 else NA,
           lwd = if (comparison_type == "line") 1.9 else 2.2,
           pt.cex = 1.15, cex = 1.03,
           bg = adjustcolor("white", alpha.f = 0.92))
  }
}

draw_sensitivity_figure <- function() {
  par(mfrow = c(3, 1), mar = c(3.7, 5.3, 2.9, 1.0),
      cex.axis = 1.15, cex.lab = 1.18, cex.main = 1.25,
      font.main = 2, mgp = c(3.0, 0.85, 0), tcl = -0.35)
  draw_sensitivity_panel(
    baseline_objects[, "real_yield_10y"],
    robust_objects[, "real_yield_10y"],
    "(a) Ten-year real yield", extended_real_yield_10y,
    "line", "Extended real-yield data", TRUE
  )
  draw_sensitivity_panel(
    baseline_objects[, "real_tp"], robust_objects[, "real_tp"],
    "(b) Ten-year real term premium", model_free_real_tp_10y,
    "points", "Survey model-free 10y"
  )
  draw_sensitivity_panel(
    baseline_objects[, "r_star"], robust_objects[, "r_star"],
    "(c) Real-rate trend"
  )
}

pdf(file.path(figure_directory, "pre2004_real_yield_sensitivity.pdf"),
    width = 10, height = 10.0, useDingbats = FALSE)
draw_sensitivity_figure()
dev.off()

cat("Liquid-TIPS robustness summary:\n")
print(summary)
cat("\nSelected parameters:\n")
print(parameter_comparison)
