# Publication-ready fit and state figures for the final six-state estimate.
# Figure-level titles are intentionally omitted; LaTeX captions carry them.

if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Rcpp is required.")
source("code/model/iekf.R")
Rcpp::sourceCpp("code/model/pricing.cpp")
Rcpp::sourceCpp("code/model/pricing_helpers.cpp")

diagnostic_directory <- Sys.getenv(
  "GGNR_DIAGNOSTIC_DIR",
  "outputs/diagnostics/baseline"
)
output_directory <- Sys.getenv(
  "GGNR_FIGURE_DIR",
  "outputs/figures"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
diagnostics <- readRDS(file.path(diagnostic_directory, "diagnostics.rds"))
data <- load_true_release_data()
dates <- data$dates
fit <- diagnostics$fitted
estimate_file <- Sys.getenv(
  "GGNR_START_FILE",
  "estimates/baseline.rds"
)
estimate <- readRDS(estimate_file)
short_rate <- load_observed_real_short_rate_inputs(dates)
hfi_none <- list(standardized = rep(NaN, length(dates)))
filter_with_covariance <- run_kf_reduced_hfi_no_output(
  data, estimate$parameters, hfi_none, short_rate,
  filter_method = if (is.null(estimate$filter_method)) "EKF" else estimate$filter_method,
  iekf_trigger_probability = if (is.null(estimate$iekf_trigger_probability)) 0 else estimate$iekf_trigger_probability,
  iekf_steps = if (is.null(estimate$iekf_steps)) 2L else estimate$iekf_steps
)
states <- filter_with_covariance$x_upd
state_sd <- t(vapply(seq_along(dates), function(tt) {
  sqrt(pmax(diag(filter_with_covariance$P_upd[, , tt]), 0))
}, numeric(ncol(states))))
colnames(state_sd) <- colnames(states)

colors <- c(model = "black", observed = "grey45", comparison = "grey30",
            comparison2 = "grey65")
publication_par <- function(...) par(
  cex.axis = 1.15, cex.lab = 1.20, cex.main = 1.15,
  font.main = 2, mgp = c(2.6, 0.8, 0), tcl = -0.35, ...
)
draw_fit <- function(observed, fitted, panel_title, legend_position = NULL) {
  observed <- 1200 * observed
  fitted <- 1200 * fitted
  limits <- range(c(observed, fitted), finite = TRUE)
  plot(dates, fitted, type = "l", col = colors["model"], lwd = 2.4,
       xlab = "", ylab = "Annualized percentage points",
       main = panel_title, ylim = limits)
  points(dates, observed, col = colors["observed"], pch = 16, cex = 0.68)
  grid(col = "grey88", lty = 1)
  if (!is.null(legend_position)) {
    legend(legend_position, c("Model", "Observed"),
           col = colors[c("model", "observed")], lty = c(1, NA),
           pch = c(NA, 16), pt.cex = 1.25, lwd = c(2.4, NA), cex = 1.05,
           bg = adjustcolor("white", alpha.f = 0.9))
  }
}

rolling_12m <- function(value) as.numeric(stats::filter(
  value, rep(1 / 12, 12), sides = 1
))

pdf(file.path(output_directory, "fig_macro_fit.pdf"),
    width = 10, height = 7.6, useDingbats = FALSE)
publication_par(mfrow = c(2, 2), mar = c(3.4, 4.8, 2.8, 0.8))
draw_fit(rolling_12m(data$macro[, 2]), rolling_12m(fit$macro[, 1]),
         "(a) Twelve-month CPI inflation", "topright")
draw_fit(data$macro[, 4], fit$macro[, 2], "(b) Perceived inflation target")
draw_fit(data$surv_infexp[, 2], fit$inflation_survey[, 2],
         "(c) Ten-year CPI forecast")
draw_fit(data$surv_tbexp[, 2], fit$tbill_survey[, 2],
         "(d) Ten-year T-bill forecast")
dev.off()

pdf(file.path(output_directory, "fig_survey_fit.pdf"),
    width = 10, height = 7.6, useDingbats = FALSE)
publication_par(mfrow = c(2, 2), mar = c(3.4, 4.8, 2.8, 0.8))
draw_fit(data$surv_infexp[, 1], fit$inflation_survey[, 1],
         "(a) One-year CPI forecast", "topright")
draw_fit(data$surv_infexp[, 2], fit$inflation_survey[, 2],
         "(b) Ten-year CPI forecast")
draw_fit(data$surv_tbexp[, 1], fit$tbill_survey[, 1],
         "(c) One-year T-bill forecast")
draw_fit(data$surv_tbexp[, 2], fit$tbill_survey[, 2],
         "(d) Ten-year T-bill forecast")
dev.off()

pdf(file.path(output_directory, "fig_nominal_yield_fit.pdf"),
    width = 10, height = 8.8, useDingbats = FALSE)
publication_par(mfrow = c(3, 1), mar = c(3.4, 4.8, 2.8, 0.8))
for (entry in list(c(1, 3), c(3, 24), c(7, 120))) {
  j <- entry[1]
  draw_fit(data$yields_n[, j], fit$nominal_yields[, j],
           paste0(entry[2], "-month nominal yield"),
           if (j == 1) "topright" else NULL)
}
dev.off()

pdf(file.path(output_directory, "fig_real_yield_fit.pdf"),
    width = 10, height = 8.8, useDingbats = FALSE)
publication_par(mfrow = c(3, 1), mar = c(3.4, 4.8, 2.8, 0.8))
for (entry in list(c(1, 24), c(2, 60), c(4, 120))) {
  j <- entry[1]
  draw_fit(data$yields_r[, j], fit$real_yields[, j],
           paste0(entry[2], "-month real yield"),
           if (j == 1) "topright" else NULL)
}
dev.off()

benchmark <- readRDS("data/external/natural_rate_comparisons.rds")
benchmark_lw <- benchmark$LW[match(dates, benchmark$date)]
benchmark_hlw <- benchmark$HLW[match(dates, benchmark$date)]
r_star <- 1200 * states[, "r_star"]
pi_star <- 1200 * states[, "pi_star"]
w <- states[, "w"]
r_star_band <- 1.96 * 1200 * state_sd[, "r_star"]
pi_star_band <- 1.96 * 1200 * state_sd[, "pi_star"]
w_band <- 1.96 * state_sd[, "w"]
pdf(file.path(output_directory, "fig_states.pdf"),
    width = 10, height = 8.8, useDingbats = FALSE)
publication_par(mfrow = c(3, 1), mar = c(3.4, 4.8, 2.8, 0.8))
limits <- range(c(r_star - r_star_band, r_star + r_star_band,
                  benchmark_hlw, benchmark_lw), finite = TRUE)
plot(dates, r_star, type = "n",
     xlab = "", ylab = "Percent", main = "(a) Real-rate trend", ylim = limits)
polygon(c(dates, rev(dates)),
        c(r_star - r_star_band, rev(r_star + r_star_band)),
        col = "grey85", border = NA)
lines(dates, r_star, lwd = 2.4, col = colors["model"])
lines(dates, benchmark_hlw, col = colors["comparison"], lty = 2, lwd = 2.0)
lines(dates, benchmark_lw, col = colors["comparison2"], lty = 3, lwd = 2.0)
grid(col = "grey90", lty = 1)
legend("bottomleft", c("Model", "95% filtered interval", "HLW", "LW"),
       col = c(colors["model"], "grey85", colors["comparison"], colors["comparison2"]),
       lty = c(1, NA, 2, 3), lwd = c(2.4, NA, 2.0, 2.0),
       pch = c(NA, 15, NA, NA), pt.cex = c(NA, 1.5, NA, NA), cex = 1.00,
       bg = adjustcolor("white", alpha.f = 0.9))
plot(dates, pi_star, type = "n", xlab = "", ylab = "Percent",
     main = "(b) Inflation trend",
     ylim = range(c(pi_star - pi_star_band, pi_star + pi_star_band,
                   1200 * data$macro[, 4]), finite = TRUE))
polygon(c(dates, rev(dates)),
        c(pi_star - pi_star_band, rev(pi_star + pi_star_band)),
        col = "grey85", border = NA)
lines(dates, pi_star, lwd = 2.4, col = colors["model"])
points(dates, 1200 * data$macro[, 4], col = colors["observed"], pch = 16, cex = 0.68)
grid(col = "grey90", lty = 1)
legend("topright", c("Model", "95% filtered interval", "PTR observations"),
       col = c(colors["model"], "grey85", colors["observed"]), lty = c(1, NA, NA),
       pch = c(NA, 15, 16), pt.cex = c(NA, 1.5, 1.25),
       lwd = c(2.4, NA, NA), cex = 1.00,
       bg = adjustcolor("white", alpha.f = 0.9))
plot(dates, w, type = "n", xlab = "", ylab = "Normalized units",
     main = "(c) Price-of-risk state",
     ylim = range(c(w - w_band, w + w_band), finite = TRUE))
polygon(c(dates, rev(dates)), c(w - w_band, rev(w + w_band)),
        col = "grey85", border = NA)
lines(dates, w, lwd = 2.4, col = colors["model"])
abline(h = 0, col = "grey65", lwd = 0.8)
grid(col = "grey90", lty = 1)
dev.off()

effective <- 1200 * data$effective_federal_funds_rate
shadow <- 1200 * (states[, "r_star"] + states[, "pi_star"] + states[, "m"])
pdf(file.path(output_directory, "fig_shadow_rate.pdf"),
    width = 10, height = 5.2, useDingbats = FALSE)
publication_par(mar = c(3.8, 4.8, 1, 0.8))
plot(dates, shadow, type = "l", lwd = 2.4, col = colors["model"],
     xlab = "", ylab = "Percent",
     ylim = range(c(shadow, effective), finite = TRUE))
lines(dates, effective, col = colors["observed"], lwd = 2.0, lty = 2)
abline(h = 0, col = "grey65", lwd = 0.8)
grid(col = "grey90", lty = 1)
legend("topright", c("Model shadow rate", "Effective federal funds rate"),
       col = colors[c("model", "observed")], lty = c(1, 2), lwd = c(2.4, 2.0),
       cex = 1.08, bg = adjustcolor("white", alpha.f = 0.9))
dev.off()

write.csv(data.frame(date = dates, r_star = r_star, pi_star = pi_star,
                     w = w, shadow_rate = shadow,
                     effective_federal_funds_rate = effective),
          file.path(output_directory, "states_and_shadow_rate.csv"),
          row.names = FALSE)
cat("Saved publication fit, state, and shadow-rate figures.\n")
