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
publication_par <- function(...) {
  par(
    cex.axis = 1.15, cex.lab = 1.20, cex.main = 1.15,
    font.main = 2, mgp = c(2.6, 0.8, 0), tcl = -0.35, ...
  )
  par(cex = 1) # Set after mfrow, which otherwise shrinks panel text.
}
draw_fit <- function(observed, fitted, panel_title, legend_position = NULL) {
  observed <- 1200 * observed
  fitted <- 1200 * fitted
  limits <- range(c(observed, fitted), finite = TRUE)
  plot(dates, fitted, type = "l", col = colors["model"], lwd = 2.4,
       xlab = "", ylab = "Percentage points",
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
    width = 10, height = 10.2, useDingbats = FALSE)
publication_par(mfrow = c(3, 2), mar = c(3.4, 4.8, 2.8, 0.8))
draw_fit(rolling_12m(data$macro[, 2]), rolling_12m(fit$macro[, 1]),
         "(a) Twelve-month CPI inflation", "topright")
draw_fit(data$macro[, 4], fit$macro[, 2], "(b) Perceived inflation target")
draw_fit(data$surv_infexp[, 1], fit$inflation_survey[, 1],
         "(c) One-year CPI forecast")
draw_fit(data$surv_infexp[, 2], fit$inflation_survey[, 2],
         "(d) Ten-year CPI forecast")
draw_fit(data$surv_tbexp[, 1], fit$tbill_survey[, 1],
         "(e) One-year T-bill forecast")
draw_fit(data$surv_tbexp[, 2], fit$tbill_survey[, 2],
         "(f) Ten-year T-bill forecast")
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

# One-sided LW and HLW estimates; frozen vintages are documented in data/README.md.
benchmark <- readRDS("data/external/natural_rate_comparisons.rds")
benchmark_lw <- benchmark$LW[match(dates, benchmark$date)]
benchmark_hlw <- benchmark$HLW[match(dates, benchmark$date)]
r_star <- 1200 * states[, "r_star"]
pi_star <- 1200 * states[, "pi_star"]
w <- states[, "w"]
# Pointwise Hamilton-style bands include filtering and parameter uncertainty.
uncertainty_directory <- Sys.getenv(
  "GGNR_HAMILTON_DIR", "outputs/diagnostics/uncertainty"
)
state_bands <- read.csv(file.path(uncertainty_directory, "state_bands.csv"))
stopifnot(identical(as.Date(state_bands$date), as.Date(dates)))
stopifnot(isTRUE(all.equal(state_bands$r_star_estimate, as.numeric(r_star))),
          isTRUE(all.equal(state_bands$pi_star_estimate, as.numeric(pi_star))),
          isTRUE(all.equal(state_bands$w_estimate, as.numeric(w))))

draw_state_panel <- function(prefix, estimate, panel_title, ylab,
                             comparisons = numeric()) {
  lower <- state_bands[[paste0(prefix, "_p025")]]
  upper <- state_bands[[paste0(prefix, "_p975")]]
  plot(dates, estimate, type = "n", xlab = "", ylab = ylab,
       main = panel_title, ylim = range(c(lower, upper, comparisons), finite = TRUE))
  polygon(c(dates, rev(dates)), c(lower, rev(upper)),
          col = "grey85", border = NA)
  grid(col = "grey92", lty = 1)
  lines(dates, estimate, lwd = 3, col = "black")
}

pdf(file.path(output_directory, "fig_states.pdf"),
    width = 8, height = 9, pointsize = 12.5, useDingbats = FALSE)
par(mfrow = c(3, 1), mar = c(2.6, 4.1, 2.3, 0.8))
# Reset the automatic mfrow text shrinkage for a three-panel figure.
par(cex = 1, cex.axis = 1.08, cex.lab = 1.10, cex.main = 1.14,
    font.main = 2, mgp = c(2.6, 0.75, 0), tcl = -0.3, las = 1)
draw_state_panel("r_star", r_star, expression(bold("(a) Real-rate trend ") * (r[t]^"*")), "Percent",
                 c(benchmark_hlw, benchmark_lw))
# Connect the available quarterly estimates; missing monthly rows otherwise
# break every line segment and make the comparison series disappear.
hlw_available <- is.finite(benchmark_hlw)
lw_available <- is.finite(benchmark_lw)
lines(dates[hlw_available], benchmark_hlw[hlw_available],
      col = colors["comparison"], lty = 2, lwd = 2.6)
lines(dates[lw_available], benchmark_lw[lw_available],
      col = colors["comparison2"], lty = 3, lwd = 2.6)
legend("bottomleft", c("Model", "95% interval", "HLW", "LW"), ncol = 2,
       col = c("black", "grey85", colors["comparison"], colors["comparison2"]),
       lty = c(1, NA, 2, 3), lwd = c(3, NA, 2.6, 2.6),
       pch = c(NA, 15, NA, NA), pt.cex = 1.5, cex = 1,
       bg = "white", box.col = "grey75")

draw_state_panel("pi_star", pi_star, expression(bold("(b) Inflation trend ") * (pi[t]^"*")), "Percent",
                 1200 * data$macro[, 4])
points(dates, 1200 * data$macro[, 4], col = colors["observed"],
       pch = 16, cex = 0.85)
legend("topright", c("Model", "95% interval", "PTR"), ncol = 3,
       col = c("black", "grey85", colors["observed"]), lty = c(1, NA, NA),
       pch = c(NA, 15, 16), pt.cex = c(1, 1.5, 1),
       lwd = c(3, NA, NA), cex = 1, bg = "white", box.col = "grey75")

draw_state_panel("w", w, expression(bold("(c) Price-of-risk state ") * (w[t])), "Normalized units")
abline(h = 0, col = "grey60", lwd = 1, lty = 3)
lines(dates, w, lwd = 3)
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
